"""Tests for SeqIO TwoBitIO module."""

import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest
import warnings
from concurrent.futures import ThreadPoolExecutor
from unittest import mock

import numpy as np
import support

import Bio
from Bio import BiopythonParserWarning
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.Seq import UndefinedSequenceError
from Bio.SeqIO import _twoBitIO  # type: ignore
from Bio.SeqRecord import SeqRecord


class Parsing(unittest.TestCase):
    """Test parsing 2bit files."""

    def setUp(self):
        path = support.DATA / "TwoBit" / "sequence.fa"
        with SeqIO.parse(path, "fasta") as records:
            self.records = list(records)

    def test_littleendian(self, step=5):
        path = support.DATA / "TwoBit" / "sequence.littleendian.2bit"
        with open(path, "rb") as stream:
            records = SeqIO.parse(stream, "twobit")
            self.assertEqual(records.byteorder, "little")
            self.assertEqual(len(self.records), len(records))
            for record1, record2 in zip(self.records, records):
                self.assertEqual(record1.id, record2.id)
                seq1 = record1.seq
                seq2 = record2.seq
                self.assertEqual(seq1, seq2)
                n = len(seq1)
                for i in range(0, n, step):
                    for j in range(i, n, step):
                        self.assertEqual(seq1[i:j], seq2[i:j])
                        self.assertEqual(repr(seq1[i:j]), repr(seq2[i:j]))

    def test_bigendian(self, step=5):
        path = support.DATA / "TwoBit" / "sequence.bigendian.2bit"
        with open(path, "rb") as stream:
            records = SeqIO.parse(stream, "twobit")
            self.assertEqual(len(records), 6)
            self.assertEqual(records.byteorder, "big")
            for record1, record2 in zip(self.records, records):
                self.assertEqual(record1.id, record2.id)
                seq1 = record1.seq
                seq2 = record2.seq
                self.assertEqual(seq1, seq2)
                n = len(seq1)
                for i in range(0, n, step):
                    for j in range(i, n, step):
                        self.assertEqual(seq1[i:j], seq2[i:j])
                        self.assertEqual(repr(seq1[i:j]), repr(seq2[i:j]))

    def test_extended_slices_and_indices(self):
        """Check slicing and indexing against an in-memory Seq."""
        path = support.DATA / "TwoBit" / "sequence.littleendian.2bit"
        with open(path, "rb") as stream:
            records = SeqIO.parse(stream, "twobit")
            for record1, record2 in zip(self.records, records):
                seq1 = record1.seq
                seq2 = record2.seq
                n = len(seq1)
                for i in range(-n - 3, n + 3):
                    if -n <= i < n:
                        self.assertEqual(seq1[i], seq2[i])
                    else:
                        with self.assertRaises(IndexError):
                            seq2[i]
                positions = (None, -n - 5, -n, -7, -1, 0, 1, 5, 7, n - 1, n, n + 5)
                steps = (None, 1, 2, 3, 4, 5, 7, -1, -2, -3, -4, -7)
                steps += (sys.maxsize, -sys.maxsize, -sys.maxsize - 1)
                for start in positions:
                    for stop in positions:
                        for step in steps:
                            key = slice(start, stop, step)
                            self.assertEqual(seq1[key], seq2[key], msg=key)

    def test_sparse_slice_reads_only_selected_span(self):
        """Check a sparse slice reads bytes up to its last element, not its stop."""
        path = support.DATA / "TwoBit" / "sequence.littleendian.2bit"
        with open(path, "rb") as stream:
            seq = next(SeqIO.parse(stream, "twobit")).seq
            with mock.patch.object(np, "fromfile", wraps=np.fromfile) as fromfile:
                for key, expected in (
                    (slice(None, None, sys.maxsize), 1),  # first base only
                    (slice(None, None, -sys.maxsize), 1),  # last base only
                    (slice(0, 100, 40), 21),  # bases 0, 40, 80: bytes 0-20
                    (slice(100, 0, -40), 21),  # bases 100, 60, 20: bytes 5-25
                ):
                    self.assertEqual(seq[key], self.records[0].seq[key], msg=key)
                    count = fromfile.call_args.kwargs["count"]
                    self.assertEqual(count, expected, msg=key)

    def test_convert_zero_step(self):
        """Check the C helper rejects a zero step instead of dividing by it."""
        blocks = np.empty((0, 2), dtype="uint32")
        with self.assertRaises(ValueError):
            _twoBitIO.convert(b"\x00", 0, 4, 0, blocks, blocks)

    def test_convert_block_ending_before_start(self):
        """Check the C helper skips N and mask blocks that end before they start.

        The parser rejects such blocks, but applyNs used to pass the negative
        length of one to memset, writing far past the end of the sequence.
        The calls run in a child process, so that on the old code a crash
        fails this test rather than the whole test run.
        """
        code = """\
import numpy as np
from Bio.SeqIO import _twoBitIO
data = bytes(4)  # 16 bases of T
wrapped = np.array([[10, 4]], dtype="uint32")  # starts at 10, ends at 4
empty = np.empty((0, 2), dtype="uint32")
for start, end, step in ((0, 16, 1), (0, 16, 2), (15, -1, -1)):
    print(_twoBitIO.convert(data, start, end, step, wrapped, empty))
    print(_twoBitIO.convert(data, start, end, step, empty, wrapped))
"""
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(
            filter(
                None,
                [os.path.dirname(os.path.dirname(Bio.__file__)), env.get("PYTHONPATH")],
            )
        )
        result = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            env=env,
            timeout=60,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        expected = [repr(b"T" * 16)] * 2 + [repr(b"T" * 8)] * 2 + [repr(b"T" * 16)] * 2
        self.assertEqual(result.stdout.splitlines(), expected, result.stderr)

    def test_sequence_long(self):
        path = support.DATA / "TwoBit" / "sequence.long.2bit"
        with open(path, "rb") as stream:
            with self.assertRaises(ValueError) as cm:
                SeqIO.parse(stream, "twobit")
            self.assertEqual(
                str(cm.exception),
                "version-1 twoBit files with 64-bit offsets for index are currently not supported",
            )

    def test_block_count_too_large(self):
        """Check a block count is checked against the file size before allocating."""
        path = support.DATA / "TwoBit" / "sequence.littleendian.2bit"
        data = path.read_bytes()
        fromfile = np.fromfile

        def checked_fromfile(file, dtype, count):
            # Fail before numpy reserves room for more values than the file has
            self.assertLessEqual(count, len(data) // 4)
            return fromfile(file, dtype=dtype, count=count)

        # The record of seq11111 starts at byte 81, so its nBlockCount is at
        # byte 85, and its maskBlockCount, after two N blocks, at byte 105.
        for field, position in (("nBlockCount", 85), ("maskBlockCount", 105)):
            for count in (100, 0x40000000, 0xFFFFFFFF):
                patched = bytearray(data)
                patched[position : position + 4] = count.to_bytes(4, "little")
                with self.subTest(field=field, count=count):
                    # np.fromfile needs a real file, not an io.BytesIO
                    with tempfile.TemporaryFile() as stream:
                        stream.write(patched)
                        stream.seek(0)
                        with mock.patch.object(np, "fromfile", checked_fromfile):
                            with self.assertRaises(ValueError) as cm:
                                SeqIO.parse(stream, "twobit")
                    self.assertEqual(
                        str(cm.exception),
                        "Found %s %u for sequence seq11111, "
                        "too many blocks for a file of 770 bytes" % (field, count),
                    )

    @contextlib.contextmanager
    def _patched_file(self, patches):
        """Yield the sample file, with some 32-bit fields patched, as a file."""
        path = support.DATA / "TwoBit" / "sequence.littleendian.2bit"
        patched = bytearray(path.read_bytes())
        for position, value in patches:
            patched[position : position + 4] = value.to_bytes(4, "little")
        # np.fromfile needs a real file, not an io.BytesIO
        with tempfile.TemporaryFile() as stream:
            stream.write(patched)
            stream.seek(0)
            yield stream

    def _parse_patched(self, patches):
        """Parse the patched sample file and read it all.

        Returns a dictionary mapping each name to its sequence as a string.
        """
        with self._patched_file(patches) as stream:
            return {
                record.id: str(record.seq) for record in SeqIO.parse(stream, "twobit")
            }

    def test_block_wrapping_past_2_32(self):
        """Check N and mask blocks whose end wraps round past 2**32 are rejected."""
        path = support.DATA / "TwoBit" / "sequence.littleendian.2bit"
        data = path.read_bytes()
        # The record of seq11111 (480 bases) starts at byte 81. Its first N
        # block starts at 77 and has its size at byte 97; its first mask
        # block has its start at byte 109 and its size at byte 117.
        mask_start = int.from_bytes(data[109:113], "little")
        cases = (
            # sizes that wrap the 32-bit end round to before the start
            ("N", 97, 77, 2**32 - 70),
            ("mask", 117, mask_start, 2**32 - 1),
        )
        for kind, position, start, size in cases:
            with self.subTest(kind=kind, size=size):
                # Raised on opening the file, before any sequence is read
                with self._patched_file([(position, size)]) as stream:
                    with self.assertRaises(ValueError) as cm:
                        SeqIO.parse(stream, "twobit")
                self.assertEqual(
                    str(cm.exception),
                    "Found %s block of size %u at position %u for sequence "
                    "seq11111, which has length 480" % (kind, size, start),
                )

    def test_block_past_end_of_sequence(self):
        """Check N and mask blocks running past the sequence end are clipped.

        Such a block is clipped to the end of its sequence with a warning, as
        the reader did before, and the other records of the file still read.
        """
        expected = self._parse_patched([])
        seq = expected["seq11111"]
        self.assertEqual(len(seq), 480)
        # The last N block of seq11111 starts at 305 with its size at byte
        # 101, its last mask block starts at 356 with its size at byte 121.
        # Each is followed by no block of its kind; mask blocks must not
        # overlap, as applyMask lower cases a base once for each block.
        masked_ns = "".join("n" if c.islower() else "N" for c in seq[305:])
        cases = (
            ("N", 101, 305, 480 - 305 + 1, seq[:305] + masked_ns),
            ("N", 101, 305, 2**32 - 306, seq[:305] + masked_ns),
            ("mask", 121, 356, 480 - 356 + 1, seq[:356] + seq[356:].lower()),
            ("mask", 121, 356, 2**32 - 357, seq[:356] + seq[356:].lower()),
        )
        for kind, position, start, size, seq11111 in cases:
            with self.subTest(kind=kind, size=size):
                with self._patched_file([(position, size)]) as stream:
                    with self.assertWarns(BiopythonParserWarning) as cm:
                        records = SeqIO.to_dict(SeqIO.parse(stream, "twobit"))
                    sequences = {name: str(r.seq) for name, r in records.items()}
                    data = records["seq11111"].seq._data
                self.assertEqual(
                    str(cm.warning),
                    "Found %s block of size %u at position %u for sequence "
                    "seq11111, which has length 480; clipping it to the end of "
                    "the sequence" % (kind, size, start),
                )
                self.assertEqual(sequences["seq11111"], seq11111)
                self.assertEqual(sequences, {**expected, "seq11111": seq11111})
                # The parser itself clips the block, not only the C reader
                table = data.nBlocks if kind == "N" else data.maskBlocks
                self.assertEqual(list(table[-1]), [start, 480])
        # A block may also start past the end; it is clipped to [480, 480]
        for kind, patches in (
            ("N", [(93, 500), (101, 10)]),
            ("mask", [(113, 500), (121, 10)]),
        ):
            with self.subTest(kind=kind, start=500):
                with self._patched_file(patches) as stream:
                    with self.assertWarns(BiopythonParserWarning):
                        records = SeqIO.to_dict(SeqIO.parse(stream, "twobit"))
                    data = records["seq11111"].seq._data
                    table = data.nBlocks if kind == "N" else data.maskBlocks
                    self.assertEqual(list(table[-1]), [480, 480])
                    self.assertEqual(len(str(records["seq11111"].seq)), 480)
        # A block may end exactly at the end of the sequence, without a warning
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            sequences = self._parse_patched([(97, 480 - 77)])
        self.assertEqual(sequences["seq11111"][77:].upper(), "N" * (480 - 77))

    def test_zero_size_blocks(self):
        """Check N and mask blocks of size zero leave the bases as they are."""
        expected = self._parse_patched([])
        seq = expected["seq11111"]
        # The N block at 77 (size at byte 97) is 19 bases of N, the mask
        # block at 178 (size at byte 117) is 16 bases of lower case.
        self.assertEqual(seq[77:96], "N" * 19)
        self.assertEqual(seq[178:194], seq[178:194].lower())
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            sequences = self._parse_patched([(97, 0)])
        self.assertNotIn("N", sequences["seq11111"][77:96].upper())
        self.assertRegex(sequences["seq11111"][77:96], "^[ACGT]+$")
        self.assertEqual(sequences["seq11111"][:77], seq[:77])
        self.assertEqual(sequences["seq11111"][96:], seq[96:])
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonParserWarning)
            sequences = self._parse_patched([(117, 0)])
        self.assertEqual(sequences["seq11111"][178:194], seq[178:194].upper())
        self.assertEqual(sequences["seq11111"][:178], seq[:178])
        self.assertEqual(sequences["seq11111"][194:], seq[194:])

    def test_overlapping_block_tables(self):
        """Check overlapping records cannot claim more blocks than the file holds."""
        # The records start one word apart in a run of 18 words holding 8,
        # followed by 18 zero words.  Each record therefore has 8 N blocks,
        # no mask blocks and a zero reserved field, and its block table fits
        # in the file, but the 17 tables overlap: together they need
        # 17 * 64 bytes, far more than the 279 bytes of the file.
        names = [b"%02d" % k for k in range(17)]
        start = 16 + 7 * len(names)
        data = (0x1A412743).to_bytes(4, "little")
        data += (0).to_bytes(4, "little")
        data += len(names).to_bytes(4, "little")
        data += (0).to_bytes(4, "little")
        for k, name in enumerate(names):
            data += bytes([len(name)]) + name + (start + 4 * k).to_bytes(4, "little")
        data += (8).to_bytes(4, "little") * 18 + bytes(4 * 18)
        self.assertEqual(len(data), 279)
        with tempfile.TemporaryFile() as stream:
            stream.write(data)
            stream.seek(0)
            with self.assertRaises(ValueError) as cm:
                SeqIO.parse(stream, "twobit")
        # The records before 04 take 4 * 64 = 256 bytes of block tables
        self.assertEqual(
            str(cm.exception),
            "Found nBlockCount 8 for sequence 04, too many blocks for a file of 279 bytes",
        )

    def test_shared_record(self):
        """Check a record shared by several names is read only once."""
        path = support.DATA / "TwoBit" / "sequence.littleendian.2bit"
        data = bytearray(path.read_bytes())
        # Point the index entry of seq222, whose offset is at byte 36, at the
        # record of seq11111 at byte 81.
        data[36:40] = (81).to_bytes(4, "little")
        with tempfile.TemporaryFile() as stream:
            stream.write(data)
            stream.seek(0)
            records = SeqIO.parse(stream, "twobit")
            self.assertEqual(len(records), 6)
            seq = records["seq222"].seq
            self.assertIs(seq, records["seq11111"].seq)
            self.assertEqual(seq, self.records[0].seq)

    def test_truncated_file(self):
        """Check a file that ends too early raises ValueError."""
        path = support.DATA / "TwoBit" / "sequence.littleendian.2bit"
        data = path.read_bytes()
        # The header takes bytes 0-15 and the index bytes 16-80. The record of
        # seq6, the last one, takes bytes 726-765, followed by its DNA. A file
        # cut before byte 766 is therefore missing part of the header, the
        # index or a record.
        sizes = (
            4,  # before the version
            14,  # in the reserved field of the header
            16,  # before the index
            20,  # in the name of seq11111
            38,  # in the offset of seq222
            726,  # before the record of seq6
            744,  # in the maskBlockCount of seq6
            764,  # in the reserved field of seq6
        )
        for size in range(766):
            with self.subTest(size=size):
                with tempfile.TemporaryFile() as stream:
                    stream.write(data[:size])
                    stream.seek(0)
                    with self.assertRaises(ValueError) as cm:
                        SeqIO.parse(stream, "twobit")
                # Other sizes give other messages, such as "Unknown signature"
                # or a block count too large for the file
                if size in sizes:
                    self.assertEqual(
                        str(cm.exception), "Unexpected end of file at byte %d" % size
                    )

    def test_file_truncated_after_parsing(self):
        """Check a file that shrinks after parsing raises ValueError on access."""
        path = support.DATA / "TwoBit" / "sequence.bigendian.2bit"
        with tempfile.TemporaryFile() as stream:
            stream.write(path.read_bytes())
            stream.seek(0)
            seq = next(SeqIO.parse(stream, "twobit")).seq
            # Keep only the first 10 bytes (40 bases) of the 480 bases of DNA
            stream.truncate(seq._data.offset + 10)
            for key in (slice(None), slice(None, None, 3), slice(100, 200), -1):
                with self.subTest(key=key):
                    with self.assertRaises(ValueError) as cm:
                        seq[key]
                    self.assertTrue(
                        str(cm.exception).startswith("Unexpected end of file at byte")
                    )
            # The DNA that is still in the file can still be read
            self.assertEqual(seq[:40], self.records[0].seq[:40])

    def test_dna_size_too_large(self):
        """Check a sequence length is checked against the file size."""
        # One record holding one byte of DNA (TCAG), with no blocks
        data = (0x1A412743).to_bytes(4, "little")
        data += (0).to_bytes(4, "little")
        data += (1).to_bytes(4, "little")
        data += (0).to_bytes(4, "little")
        data += b"\x04chr1" + (25).to_bytes(4, "little")
        for dnaSize in (4, 5, 4_000_000_000):
            record = dnaSize.to_bytes(4, "little") + bytes(12) + b"\x1b"
            with self.subTest(dnaSize=dnaSize):
                # np.fromfile needs a real file, not an io.BytesIO
                with tempfile.TemporaryFile() as stream:
                    stream.write(data + record)
                    stream.seek(0)
                    if dnaSize == 4:
                        seq = next(SeqIO.parse(stream, "twobit")).seq
                        self.assertEqual(seq, "TCAG")
                        continue
                    with self.assertRaises(ValueError) as cm:
                        SeqIO.parse(stream, "twobit")
                self.assertEqual(
                    str(cm.exception),
                    "Found dnaSize %u for sequence chr1, "
                    "too long for a file of 42 bytes" % dnaSize,
                )

    def test_sequence_count_too_large(self):
        """Check the index is not read past the end of the file."""
        # A header claiming 2**32 - 1 sequences, with no index after it
        data = (0x1A412743).to_bytes(4, "little")
        data += (0).to_bytes(4, "little")
        data += (0xFFFFFFFF).to_bytes(4, "little")
        data += (0).to_bytes(4, "little")

        class Stream(io.BytesIO):
            """Fail, rather than hang, if the parser keeps reading at the end."""

            empty_reads = 0

            def read(self, size=-1):
                data = super().read(size)
                if size and not data:
                    self.empty_reads += 1
                    if self.empty_reads > 100:
                        raise AssertionError("kept reading at end of file")
                return data

        with self.assertRaises(ValueError) as cm:
            SeqIO.parse(Stream(data), "twobit")
        self.assertEqual(str(cm.exception), "Unexpected end of file at byte 16")


class TestComparisons(unittest.TestCase):
    """Test comparisons of sequences read from 2bit files to Seq and other objects."""

    def setUp(self):
        path = support.DATA / "TwoBit" / "sequence.bigendian.2bit"
        self.stream = open(path, "rb")
        records = SeqIO.parse(self.stream, "twobit")
        record1 = next(records)
        record2 = next(records)
        self.seq1a = record1.seq
        self.seq2a = record2.seq
        path = support.DATA / "TwoBit" / "sequence.fa"
        with SeqIO.parse(path, "fasta") as records:
            record1 = next(records)
            record2 = next(records)
        self.seq1b = record1.seq
        self.seq2b = record2.seq

    def tearDown(self):
        self.stream.close()

    def test_eq(self):
        seq1a = self.seq1a
        seq2a = self.seq2a
        seq1b = self.seq1b
        seq2b = self.seq2b
        self.assertEqual(seq1a, seq1b)
        self.assertEqual(seq2a, seq2b)
        self.assertEqual(seq1a, seq1a)
        self.assertEqual(seq2a, seq2a)
        with self.assertRaises(UndefinedSequenceError):
            seq1a == Seq(None, len(seq1a))
        with self.assertRaises(UndefinedSequenceError):
            seq2a == Seq(None, len(seq2a))
        with self.assertRaises(UndefinedSequenceError):
            seq1a == Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2a == Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq1a)) == seq1a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq2a)) == seq2a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) == seq1a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) == seq2a

    def test_ne(self):
        seq1a = self.seq1a
        seq2a = self.seq2a
        seq1b = self.seq1b
        seq2b = self.seq2b
        self.assertNotEqual(seq1a, seq2a)
        self.assertNotEqual(seq1a, seq2b)
        with self.assertRaises(UndefinedSequenceError):
            seq1a != Seq(None, len(seq1a))
        with self.assertRaises(UndefinedSequenceError):
            seq2a != Seq(None, len(seq2a))
        with self.assertRaises(UndefinedSequenceError):
            seq1a != Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2a != Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq1a)) != seq1a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, len(seq2a)) != seq2a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) != seq1a
        with self.assertRaises(UndefinedSequenceError):
            Seq(None, 10) != seq2a

    def test_lt(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
        self.assertLess(seq1, seq2)
        self.assertLess("AA", seq1)
        self.assertLess(seq1, "TT")
        self.assertLess("AA", seq2)
        self.assertLess(seq2, "TTT")
        self.assertLess(b"AA", seq1)
        self.assertLess(seq1, b"TT")
        self.assertLess(b"AA", seq2)
        self.assertLess(seq2, b"TTT")
        self.assertLess(Seq("AA"), seq1)
        self.assertLess(seq1, Seq("TT"))
        self.assertLess(Seq("AA"), seq2)
        self.assertLess(seq2, Seq("TTT"))
        self.assertLess(MutableSeq("AA"), seq1)
        self.assertLess(seq1, MutableSeq("TT"))
        self.assertLess(MutableSeq("AA"), seq2)
        self.assertLess(seq2, MutableSeq("TTT"))
        with self.assertRaises(UndefinedSequenceError):
            seq1 < Seq(None, len(seq1))
        with self.assertRaises(UndefinedSequenceError):
            seq2 < Seq(None, len(seq2))
        with self.assertRaises(UndefinedSequenceError):
            seq1 < Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2 < Seq(None, 10)

    def test_le(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
        self.assertLessEqual(seq1, seq2)
        self.assertLessEqual(seq1, "TT")
        self.assertLessEqual("TT", seq2)
        self.assertLessEqual(seq1, b"TT")
        self.assertLessEqual("TT", seq2)
        self.assertLessEqual(seq1, Seq("TT"))
        self.assertLessEqual("TT", seq2)
        self.assertLessEqual(seq1, MutableSeq("TT"))
        self.assertLessEqual(MutableSeq("TT"), seq2)
        self.assertLessEqual("AA", seq1)
        self.assertLessEqual("AA", seq2)
        self.assertLessEqual(b"AA", seq1)
        self.assertLessEqual(b"AA", seq2)
        self.assertLessEqual(Seq("AA"), seq1)
        self.assertLessEqual(Seq("AA"), seq2)
        self.assertLessEqual(MutableSeq("AA"), seq1)
        self.assertLessEqual(MutableSeq("AA"), seq2)
        self.assertLessEqual("GC", seq1)
        self.assertLessEqual("GC", seq2)
        self.assertLessEqual(b"GC", seq1)
        self.assertLessEqual(b"GC", seq2)
        self.assertLessEqual(Seq("GC"), seq1)
        self.assertLessEqual(Seq("GC"), seq2)
        self.assertLessEqual(MutableSeq("GC"), seq1)
        self.assertLessEqual(MutableSeq("GC"), seq2)
        with self.assertRaises(UndefinedSequenceError):
            seq1 <= Seq(None, len(seq1))
        with self.assertRaises(UndefinedSequenceError):
            seq2 <= Seq(None, len(seq2))
        with self.assertRaises(UndefinedSequenceError):
            seq1 <= Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2 <= Seq(None, 10)

    def test_gt(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
        self.assertGreater(seq2, seq1)
        self.assertGreater("TT", seq1)
        self.assertGreater(seq2, "TT")
        self.assertGreater(b"TT", seq1)
        self.assertGreater(seq2, b"TT")
        self.assertGreater(Seq("TT"), seq1)
        self.assertGreater(seq2, Seq("TT"))
        self.assertGreater(MutableSeq("TT"), seq1)
        self.assertGreater(seq2, MutableSeq("TT"))
        self.assertGreater(seq1, "AA")
        self.assertGreater(seq2, "AA")
        self.assertGreater(seq1, b"AA")
        self.assertGreater(seq2, b"AA")
        self.assertGreater(seq1, Seq("AA"))
        self.assertGreater(seq2, Seq("AA"))
        self.assertGreater(seq1, MutableSeq("AA"))
        self.assertGreater(seq2, MutableSeq("AA"))
        self.assertGreater(seq1, "GC")
        self.assertGreater(seq2, "GC")
        self.assertGreater(seq1, b"GC")
        self.assertGreater(seq2, b"GC")
        self.assertGreater(seq1, Seq("GC"))
        self.assertGreater(seq2, Seq("GC"))
        self.assertGreater(seq1, MutableSeq("GC"))
        self.assertGreater(seq2, MutableSeq("GC"))
        with self.assertRaises(UndefinedSequenceError):
            seq1 > Seq(None, len(seq1))
        with self.assertRaises(UndefinedSequenceError):
            seq2 > Seq(None, len(seq2))
        with self.assertRaises(UndefinedSequenceError):
            seq1 > Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2 > Seq(None, 10)

    def test_ge(self):
        seq1 = self.seq1a
        seq2 = self.seq2a
        self.assertGreaterEqual(seq2, seq1)
        self.assertGreaterEqual("TT", seq1)
        self.assertGreaterEqual(seq2, "TT")
        self.assertGreaterEqual(b"TT", seq1)
        self.assertGreaterEqual(seq2, b"TT")
        self.assertGreaterEqual(Seq("TT"), seq1)
        self.assertGreaterEqual(seq2, Seq("TT"))
        self.assertGreaterEqual(MutableSeq("TT"), seq1)
        self.assertGreaterEqual(seq2, MutableSeq("TT"))
        self.assertGreaterEqual(seq1, "AA")
        self.assertGreaterEqual(seq2, "AA")
        self.assertGreaterEqual(seq1, b"AA")
        self.assertGreaterEqual(seq2, b"AA")
        self.assertGreaterEqual(seq1, Seq("AA"))
        self.assertGreaterEqual(seq2, Seq("AA"))
        self.assertGreaterEqual(seq1, MutableSeq("AA"))
        self.assertGreaterEqual(seq2, MutableSeq("AA"))
        self.assertGreaterEqual(seq1, "GC")
        self.assertGreaterEqual(seq2, "GC")
        self.assertGreaterEqual(seq1, b"GC")
        self.assertGreaterEqual(seq2, b"GC")
        self.assertGreaterEqual(seq1, Seq("GC"))
        self.assertGreaterEqual(seq2, Seq("GC"))
        self.assertGreaterEqual(seq1, MutableSeq("GC"))
        self.assertGreaterEqual(seq2, MutableSeq("GC"))
        with self.assertRaises(UndefinedSequenceError):
            seq1 >= Seq(None, len(seq1))
        with self.assertRaises(UndefinedSequenceError):
            seq2 >= Seq(None, len(seq2))
        with self.assertRaises(UndefinedSequenceError):
            seq1 >= Seq(None, 10)
        with self.assertRaises(UndefinedSequenceError):
            seq2 >= Seq(None, 10)


class TestBaseClassMethods(unittest.TestCase):
    """Test if methods from the base class are called correctly."""

    def setUp(self):
        path = support.DATA / "TwoBit" / "sequence.bigendian.2bit"
        self.stream = open(path, "rb")
        records = SeqIO.parse(self.stream, "twobit")
        self.record1_twobit = next(records)
        self.seq1_twobit = self.record1_twobit.seq
        self.record2_twobit = next(records)
        self.seq2_twobit = self.record2_twobit.seq
        path = support.DATA / "TwoBit" / "sequence.fa"
        with SeqIO.parse(path, "fasta") as records:
            self.record1_fasta = next(records)
            self.record2_fasta = next(records)
        self.seq1_fasta = self.record1_fasta.seq
        self.seq2_fasta = self.record2_fasta.seq

    def tearDown(self):
        self.stream.close()

    def test_getitem(self):
        self.assertEqual(self.seq1_twobit, self.seq1_fasta)
        self.assertEqual(self.seq2_twobit, self.seq2_fasta)
        self.assertEqual(self.seq1_twobit[:], self.seq1_fasta[:])
        self.assertEqual(self.seq2_twobit[:], self.seq2_fasta[:])
        self.assertEqual(self.seq1_twobit[30], self.seq1_fasta[30])
        self.assertEqual(self.seq2_twobit[30], self.seq2_fasta[30])
        self.assertEqual(self.seq1_twobit[-30], self.seq1_fasta[-30])
        self.assertEqual(self.seq2_twobit[-30], self.seq2_fasta[-30])
        self.assertEqual(self.record1_twobit.seq, self.record1_fasta.seq)
        self.assertEqual(self.record2_twobit.seq, self.record2_fasta.seq)
        self.assertEqual(self.record1_twobit[:].seq, self.record1_fasta[:].seq)
        self.assertEqual(self.record2_twobit[:].seq, self.record2_fasta[:].seq)
        self.assertEqual(self.record1_twobit[30], self.record1_fasta[30])
        self.assertEqual(self.record2_twobit[30], self.record2_fasta[30])
        self.assertEqual(self.record1_twobit[-30], self.record1_fasta[-30])
        self.assertEqual(self.record2_twobit[-30], self.record2_fasta[-30])

    def test_bytes(self):
        b = bytes(self.seq1_twobit)
        self.assertIsInstance(b, bytes)
        self.assertEqual(len(b), 480)
        self.assertEqual(b, bytes(self.seq1_fasta))
        b = bytes(self.seq1_twobit[:10])
        self.assertEqual(len(b), 10)
        self.assertIsInstance(b, bytes)
        self.assertEqual(b, b"GTATACCCCT")

    def test_hash(self):
        self.assertEqual(hash(self.seq1_twobit), hash(self.seq1_fasta))

    def test_add(self):
        self.assertIsInstance(self.seq1_twobit + "ABCD", Seq)
        self.assertEqual(self.seq1_twobit + "ABCD", self.seq1_fasta + "ABCD")
        self.assertIsInstance(self.record1_twobit + "ABCD", SeqRecord)
        record1_twobit = self.record1_twobit + "ABCD"
        record1_fasta = self.record1_fasta + "ABCD"
        self.assertEqual(self.record1_twobit.seq, self.record1_fasta.seq)

    def test_radd(self):
        self.assertIsInstance("ABCD" + self.seq1_twobit, Seq)
        self.assertEqual("ABCD" + self.seq1_twobit, "ABCD" + self.seq1_fasta)
        self.assertIsInstance("ABCD" + self.record1_twobit, SeqRecord)
        record1_twobit = "ABCD" + self.record1_twobit
        record1_fasta = "ABCD" + self.record1_fasta
        self.assertEqual(self.record1_twobit.seq, self.record1_fasta.seq)

    def test_mul(self):
        self.assertIsInstance(2 * self.seq1_twobit, Seq)
        self.assertEqual(2 * self.seq1_twobit, 2 * self.seq1_fasta)
        self.assertIsInstance(self.seq1_twobit * 2, Seq)
        self.assertEqual(self.seq1_twobit * 2, self.seq1_fasta * 2)

    def test_contains(self):
        for seq in (
            self.seq1_twobit,
            self.seq1_fasta,
            self.record1_twobit,
            self.record1_fasta,
        ):
            self.assertIn("ACCCCT", seq)
            self.assertNotIn("ACGTACGT", seq)

    def test_repr(self):
        self.assertIsInstance(repr(self.seq1_twobit), str)
        self.assertEqual(repr(self.seq1_twobit), repr(self.seq1_fasta))

    def test_str(self):
        self.assertIsInstance(str(self.seq1_twobit), str)
        self.assertEqual(str(self.seq1_twobit), str(self.seq1_fasta))

    def test_count(self):
        self.assertEqual(self.seq1_twobit.count("CT"), self.seq1_fasta.count("CT"))
        self.assertEqual(
            self.seq1_twobit.count("CT", 75), self.seq1_fasta.count("CT", 75)
        )
        self.assertEqual(
            self.seq1_twobit.count("CT", 125, 250),
            self.seq1_fasta.count("CT", 125, 250),
        )
        self.assertEqual(
            self.record1_twobit.count("CT"), self.record1_fasta.count("CT")
        )
        self.assertEqual(
            self.record1_twobit.count("CT", 75), self.record1_fasta.count("CT", 75)
        )
        self.assertEqual(
            self.record1_twobit.count("CT", 125, 250),
            self.record1_fasta.count("CT", 125, 250),
        )

    def test_find(self):
        self.assertEqual(self.seq1_twobit.find("CT"), self.seq1_fasta.find("CT"))
        self.assertEqual(
            self.seq1_twobit.find("CT", 75), self.seq1_fasta.find("CT", 75)
        )
        self.assertEqual(
            self.seq1_twobit.find("CT", 75, 100), self.seq1_fasta.find("CT", 75, 100)
        )
        self.assertEqual(
            self.seq1_twobit.find("CT", None, 100),
            self.seq1_fasta.find("CT", None, 100),
        )

    def test_rfind(self):
        self.assertEqual(self.seq1_twobit.rfind("CT"), self.seq1_fasta.rfind("CT"))
        self.assertEqual(
            self.seq1_twobit.rfind("CT", 450), self.seq1_fasta.rfind("CT", 450)
        )
        self.assertEqual(
            self.seq1_twobit.rfind("CT", None, 100),
            self.seq1_fasta.rfind("CT", None, 100),
        )
        self.assertEqual(
            self.seq1_twobit.rfind("CT", 75, 100), self.seq1_fasta.rfind("CT", 75, 100)
        )

    def test_index(self):
        self.assertEqual(self.seq1_twobit.index("CT"), self.seq1_fasta.index("CT"))
        self.assertEqual(
            self.seq1_twobit.index("CT", 75), self.seq1_fasta.index("CT", 75)
        )
        self.assertEqual(
            self.seq1_twobit.index("CT", None, 100),
            self.seq1_fasta.index("CT", None, 100),
        )
        for seq in (self.seq1_twobit, self.seq1_fasta):
            self.assertRaises(ValueError, seq.index, "CT", 75, 100)

    def test_rindex(self):
        self.assertEqual(self.seq1_twobit.rindex("CT"), self.seq1_fasta.rindex("CT"))
        self.assertEqual(
            self.seq1_twobit.rindex("CT", None, 100),
            self.seq1_fasta.rindex("CT", None, 100),
        )
        for seq in (self.seq1_twobit, self.seq1_fasta):
            self.assertRaises(ValueError, seq.rindex, "CT", 450)
            self.assertRaises(ValueError, seq.rindex, "CT", 75, 100)

    def test_startswith(self):
        for seq in (self.seq1_twobit, self.seq1_fasta):
            self.assertTrue(seq.startswith("GTAT"))
            self.assertTrue(seq.startswith("TGGG", start=10))
            self.assertTrue(seq.startswith("TGGG", start=10, end=14))
            self.assertFalse(seq.startswith("TGGG", start=10, end=12))

    def test_endswith(self):
        for seq in (self.seq1_twobit, self.seq1_fasta):
            self.assertTrue(seq.endswith("ACCG"))
            self.assertTrue(seq.endswith("ACCG", 476))
            self.assertTrue(seq.endswith("GCAC", 472, 478))
            self.assertFalse(seq.endswith("GCAC", 476, 478))

    def test_split(self):
        self.assertEqual(self.seq1_twobit.split(), self.seq1_fasta.split())
        self.assertEqual(self.seq1_twobit.split("C"), self.seq1_fasta.split("C"))
        self.assertEqual(self.seq1_twobit.split("C", 1), self.seq1_fasta.split("C", 1))

    def test_rsplit(self):
        self.assertEqual(self.seq1_twobit.rsplit(), self.seq1_fasta.rsplit())
        self.assertEqual(self.seq1_twobit.rsplit("C"), self.seq1_fasta.rsplit("C"))
        self.assertEqual(
            self.seq1_twobit.rsplit("C", 1), self.seq1_fasta.rsplit("C", 1)
        )

    def test_strip(self):
        self.assertEqual(self.seq1_twobit.strip("G"), self.seq1_fasta.strip("G"))

    def test_lstrip(self, chars=None):
        self.assertEqual(self.seq1_twobit.lstrip("G"), self.seq1_fasta.lstrip("G"))

    def test_rstrip(self, chars=None):
        self.assertEqual(self.seq1_twobit.rstrip("G"), self.seq1_fasta.rstrip("G"))

    def test_upper(self):
        seq1_twobit_upper = self.seq1_twobit.upper()
        seq1_fasta_upper = self.seq1_fasta.upper()
        self.assertEqual(seq1_twobit_upper, seq1_fasta_upper)
        self.assertEqual(seq1_twobit_upper[140:210], seq1_fasta_upper[140:210])
        seq2_twobit_upper = self.seq2_twobit.upper()
        seq2_fasta_upper = self.seq2_fasta.upper()
        self.assertEqual(seq2_twobit_upper, seq2_fasta_upper)
        self.assertEqual(seq2_twobit_upper[140:210], seq2_fasta_upper[140:210])

    def test_lower(self):
        seq1_twobit_lower = self.seq1_twobit.lower()
        seq1_fasta_lower = self.seq1_fasta.lower()
        self.assertEqual(seq1_twobit_lower, seq1_fasta_lower)
        self.assertEqual(seq1_twobit_lower[140:210], seq1_fasta_lower[140:210])
        seq2_twobit_lower = self.seq2_twobit.lower()
        seq2_fasta_lower = self.seq2_fasta.lower()
        self.assertEqual(seq2_twobit_lower, seq2_fasta_lower)
        self.assertEqual(seq2_twobit_lower[140:210], seq2_fasta_lower[140:210])

    def test_isupper(self):
        self.assertEqual(self.seq1_twobit.isupper(), self.seq1_fasta.isupper())
        self.assertEqual(self.seq2_twobit.isupper(), self.seq2_fasta.isupper())
        self.assertEqual(self.record1_twobit.isupper(), self.record1_fasta.isupper())
        self.assertEqual(self.record2_twobit.isupper(), self.record2_fasta.isupper())

    def test_islower(self):
        self.assertEqual(self.seq1_twobit.islower(), self.seq1_fasta.islower())
        self.assertEqual(self.seq2_twobit.islower(), self.seq2_fasta.islower())
        self.assertEqual(self.record1_twobit.islower(), self.record1_fasta.islower())
        self.assertEqual(self.record2_twobit.islower(), self.record2_fasta.islower())

    def test_replace(self):
        # seq.transcribe uses seq._data.replace
        self.assertEqual(self.seq1_twobit.transcribe(), self.seq1_fasta.transcribe())

    def test_translate(self):
        # seq.reverse_complement uses seq._data.translate
        self.assertEqual(
            self.seq1_twobit.reverse_complement(), self.seq1_fasta.reverse_complement()
        )
        record1_twobit = self.record1_twobit.reverse_complement()
        record1_fasta = self.record1_fasta.reverse_complement()
        self.assertEqual(record1_twobit.seq, record1_fasta.seq)

    def test_defined(self):
        self.assertTrue(self.seq1_twobit.defined)
        self.assertTrue(self.seq2_twobit.defined)
        self.assertEqual(self.seq1_twobit.defined_ranges, ((0, len(self.seq1_twobit)),))
        self.assertEqual(self.seq2_twobit.defined_ranges, ((0, len(self.seq2_twobit)),))


# Run only while the GIL is off. A free-threaded build turns it back on when
# an undeclared extension is imported, unless PYTHON_GIL=0 is set.
@unittest.skipIf(
    getattr(sys, "_is_gil_enabled", lambda: True)(), "requires the GIL to be disabled"
)
class ConcurrentConvert(unittest.TestCase):
    """Concurrent calls of the C helper on shared inputs match serial calls.

    The threads share the packed bytes and the N and mask block arrays. They
    call _twoBitIO.convert directly, because a twoBit-backed Seq seeks and
    reads its file on access, and threads sharing one would race on the
    file position; that is pure Python, and not what is tested here.
    """

    def test_threads_match_serial(self):
        calls = []
        path = support.DATA / "TwoBit" / "sequence.littleendian.2bit"
        with open(path, "rb") as stream:
            for record in SeqIO.parse(stream, "twobit"):
                data = record.seq._data
                n = len(data)
                stream.seek(data.offset)
                packed = stream.read((n + 3) // 4)
                # Forward and reverse, with steps that take both C paths.
                for start, end, step in (
                    (0, n, 1),
                    (0, n, 3),
                    (n - 1, -1, -1),
                    (n - 1, -1, -2),
                ):
                    calls.append(
                        (packed, start, end, step, data.nBlocks, data.maskBlocks)
                    )
        expected = [_twoBitIO.convert(*args) for args in calls]
        rounds = 50

        def work(_):
            return [_twoBitIO.convert(*args) for args in calls * rounds]

        with ThreadPoolExecutor(max_workers=4) as executor:
            for results in executor.map(work, range(16)):
                self.assertEqual(results, expected * rounds)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
