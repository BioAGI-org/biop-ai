"""Tests for SeqIO GFA module."""

import os
import subprocess
import sys
import unittest
import warnings
from io import StringIO

import support

import Bio
from Bio import BiopythonParserWarning
from Bio import BiopythonWarning
from Bio import SeqIO


class TestRead(unittest.TestCase):
    def test_read_GFA1(self):
        """Test parsing valid GFA 1.x files."""
        records = list(SeqIO.parse(support.DATA / "GFA" / "seq.gfa", "gfa1"))
        self.assertEqual(len(records), 8)
        self.assertEqual(records[6].id, "MTh13014")
        self.assertEqual(
            records[6].seq,
            "TTAGGTCTCCACCCCTGACTCCCCTCAGCCATAGAAGGCCCCACCCCAGTCTCAGCCCTACTCCACTCAAGCACTATAGTTGTAGCAGGAATCTTCTTACTCATCCGCTTCCACCCCCTAGCAGAAAATAGCCCACTAATCCAAACTCTAACACTATGCTTAGGCGCTATCACCACTCTGTTCGCAGCAGTCTGCGCCCTTACACAAAATGACATCAAAAAAATCGTAGCCTTCTCCACTTCAAGTCAACTAGGACTCATAATAGTTACAATCGGCATCAACCAACCACACCTAGCATTCCTGCACATCTGTACCCACGCCTTCTTCAAAGCCATACTATTTATGTGCTCCGGGTCCATCATCCACAACCTTAACAATGAACAAGATATTCGAAAAATAGGAGGACTACTCAAAACCATACCTCTCACTTCAACCTCCCTCACCATTGGCAGCCTAGCATTAGCAGGAATACCTTTCCTCACAGGTTTCTACTCCAAAGACC",
        )
        self.assertEqual(records[0].annotations["SN"], ("Z", "MT_human"))
        self.assertEqual(records[0].annotations["SO"], ("i", "0"))

        records = list(SeqIO.parse(support.DATA / "GFA" / "seq_with_len.gfa", "gfa1"))
        self.assertEqual(len(records), 9)
        self.assertEqual(
            records[8].seq,
            "GAAAAATTGCCCTTGGTTTTCGCTTCGCTCAAACTCTATTGAACTTCGCTTTCGCTCAGTTCGTCGGGGCAATTTTTTGGTTAATACTT",
        )

        records = list(
            SeqIO.parse(support.DATA / "GFA" / "fake_with_checksum.gfa", "gfa1")
        )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].seq, "AAA")

        records = list(SeqIO.parse(support.DATA / "GFA" / "no_seq.gfa", "gfa1"))
        self.assertEqual(len(records), 9)
        self.assertEqual(len(records[0]), 528)

    def test_read_GFA2(self):
        """Test parsing valid GFA 2.0 files."""
        records = list(SeqIO.parse(support.DATA / "GFA" / "fake_gfa2.gfa", "gfa2"))
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].seq, "AAA")

    def test_read_GFA2_no_seq(self):
        """Test a GFA 2.0 segment without sequence keeps its length."""
        record = SeqIO.read(StringIO("S\ts1\t100\t*\n"), "gfa2")
        self.assertEqual(len(record), 100)
        self.assertFalse(record.seq.defined)

    def test_read_GFA2_no_seq_LN(self):
        """Test LN gives the length of a GFA 2.0 segment without sequence."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            record = SeqIO.read(StringIO("S\ts1\t100\t*\tLN:i:528\n"), "gfa2")
        self.assertEqual(len(record), 528)
        self.assertFalse(record.seq.defined)

    def test_read_GFA2_no_seq_checksum(self):
        """Test a GFA 2.0 segment without sequence keeps an unverifiable checksum."""
        checksum = "ABCD" * 16
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            record = SeqIO.read(StringIO(f"S\ts1\t100\t*\tSH:H:{checksum}\n"), "gfa2")
        self.assertEqual(len(record), 100)
        self.assertEqual(record.annotations["SH"], ("H", checksum))


class TestCorrupt(unittest.TestCase):
    def test_corrupt_gfa2(self):
        """Check a GFA 1.x file does not parse in GFA 2."""
        with self.assertRaises(ValueError):
            list(SeqIO.parse(support.DATA / "GFA" / "seq.gfa", "gfa2"))

    def test_corrupt_segment_fields(self):
        """Check a GFA file with invalid fields on a segment line."""
        with self.assertRaises(ValueError):
            list(
                SeqIO.parse(support.DATA / "GFA" / "corrupt_segment_fields.gfa", "gfa1")
            )

    def test_corrupt_len(self):
        """Check a GFA file with an incorrect length."""
        with self.assertWarns(BiopythonParserWarning):
            list(SeqIO.parse(support.DATA / "GFA" / "corrupt_len.gfa", "gfa1"))

    def test_corrupt_repeated_len(self):
        """Check only the first LN tag sets the length of a segment without sequence."""
        for fmt, line in (
            ("gfa1", "S\ts1\t*\tLN:i:5\tLN:i:6\n"),
            ("gfa2", "S\ts1\t100\t*\tLN:i:5\tLN:i:6\n"),
        ):
            with self.subTest(fmt=fmt):
                with self.assertWarnsRegex(BiopythonParserWarning, "incorrect length"):
                    record = SeqIO.read(StringIO(line), fmt)
                self.assertEqual(len(record), 5)

    def test_corrupt_negative_len_gfa2(self):
        """Check a GFA 2.0 segment without sequence cannot have a negative length."""
        with self.assertRaisesRegex(ValueError, "non-negative length: S\ts1\t-5\t\\*"):
            SeqIO.read(StringIO("S\ts1\t-5\t*\n"), "gfa2")

    def test_corrupt_checksum(self):
        """Check a GFA file with an incorrect checksum."""
        with self.assertWarns(BiopythonParserWarning):
            list(SeqIO.parse(support.DATA / "GFA" / "corrupt_checksum.gfa", "gfa1"))

    def test_corrupt_tag_name(self):
        """Check a GFA file with an invalid tag name."""
        with self.assertWarns(BiopythonParserWarning):
            list(SeqIO.parse(support.DATA / "GFA" / "corrupt_tag_name.gfa", "gfa1"))

    def test_corrupt_tag_type(self):
        """Check a GFA file with an incorrect tag type."""
        with self.assertWarns(BiopythonParserWarning):
            list(SeqIO.parse(support.DATA / "GFA" / "corrupt_tag_type.gfa", "gfa1"))

    def test_parser_warnings(self):
        """Check each malformed-input warning is a BiopythonParserWarning."""
        for fmt, data, message in (
            ("gfa1", "\nS\ts1\tAAA\n", "blank line"),
            ("gfa2", "\nS\ts1\t3\tAAA\n", "blank line"),
            ("gfa1", "S\ts1\tAAA\tLN:i:5\n", "incorrect length"),
            ("gfa1", f"S\ts1\tAAA\tSH:H:{'0' * 64}\n", "incorrect checksum"),
            ("gfa1", "S\ts1\tAAA\t~~:i:0\n", "invalid name"),
            ("gfa1", "S\ts1\tAAA\tAB:X:0\n", "invalid type"),
            ("gfa1", "S\ts1\tAAA\tAB::0\n", "invalid type"),
            ("gfa1", "S\ts1\tAAA\tAB:if:0\n", "invalid type"),
            ("gfa1", "S\ts1\tAAA\tAB:A:ab\n", "printable character"),
            ("gfa1", "S\ts1\tAAA\tAB:i:C\n", "signed integer"),
            ("gfa1", "S\ts1\tAAA\tAB:f:C\n", "Expected float"),
            ("gfa1", "S\ts1\tAAA\tAB:Z:\u00e9\n", "printable string"),
            ("gfa1", "S\ts1\tAAA\tAB:J:\u00e9\n", "JSON"),
            ("gfa1", "S\ts1\tAAA\tAB:H:ab\n", "hex"),
            ("gfa1", "S\ts1\tAAA\tAB:B:C\n", "integers or floats"),
        ):
            with self.subTest(fmt=fmt, data=data):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    SeqIO.read(StringIO(data), fmt)
                self.assertEqual(len(caught), 1)
                self.assertIs(caught[0].category, BiopythonParserWarning)
                self.assertIn(message, str(caught[0].message))

    def test_float_tag_backtracking(self):
        """Check long f and B tags are type checked in linear time.

        The float pattern given by the GFA standard can match a run of digits
        in many ways, and Python's re tried them all on a value that did not
        match. That took time quadratic in the length of an f tag, and doubling
        with each element of a B tag, so a B tag of 40 elements took days. It
        still gave the right warnings in the end, so a long bad tag of each
        kind is checked in a child process, and the timeout fails the test on
        the old pattern.
        """
        # The old pattern took 8.2 seconds on an f tag of 20,000 digits and
        # "a", four times as long for each doubling, so 200,000 digits take
        # it about 14 minutes. It took 6.7 seconds on a B tag of 24 elements
        # "11" and "X", twice as long for each element, so 30 elements take
        # it about 7 minutes. The new pattern takes milliseconds on each.
        code = """\
import sys
import warnings
from Bio import SeqIO
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    for record in SeqIO.parse(sys.stdin, "gfa1"):
        print(record.id, [str(w.message).split(",")[0] for w in caught])
        caught.clear()
"""
        lines = [
            "S\tf\tAAA\tAB:f:" + "1" * 200_000 + "a\n",
            "S\tB\tAAA\tAB:B:f," + ",".join(["11"] * 30) + "X\n",
        ]
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join(
            filter(
                None,
                [os.path.dirname(os.path.dirname(Bio.__file__)), env.get("PYTHONPATH")],
            )
        )
        result = subprocess.run(
            [sys.executable, "-c", code],
            input="".join(lines),
            timeout=60,
            capture_output=True,
            text=True,
            env=env,
        )
        float_warning = "Tag has incorrect type. Expected float"
        array_warning = "Tag has incorrect type. Expected array of integers or floats"
        self.assertEqual(
            result.stdout,
            f"f {[float_warning]}\nB {[array_warning]}\n",
            result.stderr,
        )

        # Shorter bad tags give the same warnings, and long good tags none.
        array = ",".join(["11"] * 10)
        for value, message in (
            ("f:" + "1" * 1000 + "a", "Expected float"),
            ("B:f," + array + "X", "integers or floats"),
            ("f:" + "1" * 60000, None),
            ("B:f," + ",".join(["11"] * 28), None),
            ("B:i,-1.5,+.5,2e-3,7E+10," + array, None),
        ):
            with self.subTest(value=value[:12], message=message):
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter("always")
                    record = SeqIO.read(StringIO(f"S\ts1\tAAA\tAB:{value}\n"), "gfa1")
                self.assertEqual(record.annotations["AB"], tuple(value.split(":")))
                if message is None:
                    self.assertEqual(caught, [])
                else:
                    self.assertEqual(len(caught), 1)
                    self.assertIs(caught[0].category, BiopythonParserWarning)
                    self.assertIn(message, str(caught[0].message))

    def test_parser_warnings_filter(self):
        """Check a filter on BiopythonWarning still applies to GFA warnings."""
        with warnings.catch_warnings():
            warnings.simplefilter("error", BiopythonWarning)
            with self.assertRaises(BiopythonParserWarning):
                SeqIO.read(StringIO("S\ts1\tAAA\tLN:i:5\n"), "gfa1")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
