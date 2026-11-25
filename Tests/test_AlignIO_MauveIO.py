# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# Revisions copyright 2011 Brandon Invergo. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.AlignIO.MauveIO module."""

import os
import unittest
from io import StringIO
from types import MappingProxyType

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO.MauveIO import MauveIterator
from Bio.AlignIO.MauveIO import MauveWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class TestMauveIO(unittest.TestCase):
    MAUVE_TEST_DATA_DIR = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "Mauve"
    )
    SIMPLE_XMFA = os.path.join(MAUVE_TEST_DATA_DIR, "simple.xmfa")
    SIMPLE_FA = os.path.join(MAUVE_TEST_DATA_DIR, "simple.fa")

    def test_one(self):
        ids = []
        with open(self.SIMPLE_XMFA) as handle:
            for alignment in MauveIterator(handle):
                for record in alignment:
                    ids.append(record.id)

        self.assertEqual(
            ids,
            [
                "1/0-5670",
                "2/0-5670",
                "1/5670-9940",
                "2/7140-11410",
                "1/9940-14910",
                "2/5670-7140",
                "2/11410-12880",
            ],
        )

        expected = """ATTCGCACAT AAGAATGTAC CTTGCTGTAA TTTATACTCA
            GCAGGTGGTG CAGACATCAT AACAAAAGAA GACTCTTGTT GTACTAGATA TTGTGTAGCA
            TCACGACCAC ACACACATGG AATGGAAACA CCTGTCTTAA GATTATCATA AGATAGAGTA
            CCCATATACA TCACAGCTTC TACACCCGTT AAGGTAGTAG TTTTCTGACC ACAATGTTTA
            CACACCACAT TAAGAACTCG CTTTGCAGAT TCCAAATTAG CATGCTGTAG AAGATGGGTC
            ATAGTTTCTC TGACATCACC AAGCTCGCCA ACAGTTTTAT TACTGTAAGC GAGTATGAGT
            GCACAAAAGT TAGCAGCATC ACCAGCACGG GCTCTATAAT AAGCCTCTTG AAGTGCTGGT
            GCATTGAATT TGACTTCAAG CTGTTGAAGT GCTAATAAAA CACTAGACAA ATAACAATTG
            TTATCAGCCC ATTTAATTGA AGTTAAACCA CCAACTTGAG GAAATTTCCA TTTCTTTGTG
            TGGTTTAAAG CAGACATGTA CCTACCAAGA AAACTCTCAT CAAGAGTATG GTAGTACTCG
            AAAGCTTCAC TACGTAGTGT GTCATCACTA GGTAGTACAA AGAAAGTCTT ACCCTCATGA
            TTTACATGAG GTTTAATTTT TGTAACATCA GCACCATCCA AGTATGTTGG ACCAAACTGC
            TGTCCATATG TCATAGACAT ATCCACAAGC TGTGTGTGGA GATTAGTGTT GTCCACAGTT
            GTGAACACTT TTATAGTCTT AACCTCCCGC AGGGATAAGA GACTCTTTAG TTTGTCAAGT
            GAAAGAACCT CACCGTCAAG ATGAAACTCG ACGGGGCTCT CCAGAGTGTG GTACACAATT
            TTGTCACCAC GCTTAAGAAA TTCAACACCT AACTCTGTAC GCTGTCCTGA ATAGGACCAA
            TCTCTGTAAG AGCCAGCCAA AGAAACTGTT TCTACAAAGT GCTCCTCAGA TGTCTTTGAT
            GACGAAGTGA GGTATCCATT ATATGTAGTA ACAGCATCTG GTGATGATAC TGACACTACG
            GCAGGAGCTT TAAGAGAACG CATACAGCGC GCAGCCTCTT CAAGATTAAA ACCATGTGTC
            ACATAACCAA TTGGCATTGT GACAAGCGGC TCATTTAGAG AGTTCAGCTT CGTAATAATA
            GAAGCTACAG GCTCTTTACT AGTATAAAAG AAGAATCGGA CACCATAGTC AACGATGCCC
            TCTTGAATTT TAATTCCTTT ATACTTACGT TGGATGGTTG CCATTATGGC TCTAACATCC
            ATGCATATAG GCATTAATTT TCTTGTCTCT TCAGCATGAG CAAGCATTTC TCTCAAATTC
            CAGGATACAG TTCCTAGAAT CTCTTCCTTA GCATTAGGTG CTTCTGAAGG TAGTACATAA
            AATGCAGATT TGCATTTCTT AAGAGCAGTC TTAGCTTCCT CAAGTGTATA """
        self.assertEqual(
            str(record.seq).replace("-", ""),
            expected.replace(" ", "").replace("\n", ""),
        )

    def test_sequence_positions(self):
        with open(self.SIMPLE_FA) as handle:
            seqs = list(SeqIO.parse(handle, "fasta"))

        with open(self.SIMPLE_XMFA) as handle:
            aln_list = list(MauveIterator(handle))

        for aln in aln_list:
            for record in aln:
                if not record.seq.startswith("-"):
                    expected = record.seq[0:10]
                    # seqs 0, 1 are ids 1, 2
                    actual = seqs[int(record.name) - 1].seq
                    # Slice out portion mentioned in file
                    actual = actual[record.annotations["start"] : record.annotations["end"]]  # fmt: skip

                    if record.annotations["strand"] < 0:
                        actual = actual.reverse_complement()
                    # Slice first 10 chars for comparison, don't want to
                    # get any '-'s by accident
                    actual = actual[0:10]
                    if len(actual) == 0:
                        # We can't test sequences which don't provide
                        # proper annotation start/end/strand information
                        continue
                    self.assertEqual(expected, actual)

    def test_write_read(self):
        with open(self.SIMPLE_XMFA) as handle:
            aln_list = list(MauveIterator(handle))

        handle = StringIO()
        MauveWriter(handle).write_file(aln_list)
        handle.seek(0)
        aln_list_out = list(MauveIterator(handle))

        for a1, a2 in zip(aln_list, aln_list_out):
            self.assertEqual(len(a1), len(a2))
            for r1, r2 in zip(a1, a2):
                self.assertEqual(r1.id, r2.id)
                self.assertEqual(r1.seq, r2.seq)

    def test_write_first_block_lists_every_sequence(self):
        """The first block written lists every sequence, so rows keep order."""
        with open(self.SIMPLE_XMFA) as handle:
            aln_list = list(MauveIterator(handle))
        handle = StringIO()
        # The first block written holds only sequence 2
        MauveWriter(handle).write_file(aln_list[::-1])
        self.assertIn("#Sequence2Entry\t2\n", handle.getvalue())
        handle.seek(0)
        aln_list_out = list(MauveIterator(handle))
        self.assertEqual([r.id for r in aln_list_out[0]], ["1", "2/11410-12880"])
        self.assertEqual(set(aln_list_out[0][0].seq), {"-"})
        self.assertEqual(
            [r.id for r in aln_list_out[3]], ["1/5670-9940", "2/7140-11410"]
        )
        # Two placeholders in the first block are both kept
        data = (
            "> 1:0-0 + a.fa\n> 2:0-0 + b.fa\n> 3:1-4 + c.fa\nACGT\n=\n"
            "> 1:1-4 + a.fa\nAAAA\n> 2:1-4 + b.fa\nCCCC\n=\n"
        )
        handle = StringIO()
        MauveWriter(handle).write_file(MauveIterator(StringIO(data)))
        handle.seek(0)
        alignment = next(MauveIterator(handle))
        self.assertEqual([r.id for r in alignment], ["a.fa", "b.fa", "c.fa/0-4"])

    def test_write_reused_alignment(self):
        """Each alignment from an iterator that reuses one object is kept."""
        alignment = MultipleSeqAlignment([SeqRecord(Seq("AAAA"), id="1", name="1")])

        def blocks():
            yield alignment
            alignment[0].seq = Seq("CCCC")
            yield alignment

        handle = StringIO()
        MauveWriter(handle).write_file(blocks())
        handle.seek(0)
        self.assertEqual(
            [str(a[0].seq) for a in MauveIterator(handle)], ["AAAA", "CCCC"]
        )

    def test_write_uncopyable_annotation(self):
        """An unused annotation that cannot be copied does not stop writing."""
        record = SeqRecord(Seq("AAAA"), id="1", name="1")
        record.annotations["metadata"] = MappingProxyType({"key": "value"})
        handle = StringIO()
        MauveWriter(handle).write_file([MultipleSeqAlignment([record])])
        handle.seek(0)
        self.assertEqual([str(a[0].seq) for a in MauveIterator(handle)], ["AAAA"])

    def test_write_not_an_alignment(self):
        """Something that is not an alignment is rejected by write_alignment."""
        # Not an UndefinedSequenceError from numbering its "rows" first.
        record = SeqRecord(Seq(None, length=4), id="a")
        with self.assertRaisesRegex(AttributeError, "get_alignment_length"):
            MauveWriter(StringIO()).write_file([record])

    def test_malformed_header(self):
        """A malformed header line is named in the error message."""
        handle = StringIO("#FormatVersion Mauve1\n> nonsense\nACGT\n=\n")
        with self.assertRaisesRegex(ValueError, "^Malformed header line: > nonsense$"):
            list(MauveIterator(handle))

    def test_ids_not_shared_between_iterators(self):
        """Each iterator orders records by its own file, not earlier ones."""
        a = "> 1:1-4 + a.fa\nACGT\n> 2:1-4 + b.fa\nACGA\n=\n"
        b = "> 2:1-4 + b.fa\nTTTT\n> 1:1-4 + a.fa\nGGGG\n=\n"
        next(MauveIterator(StringIO(a)))
        iterator = MauveIterator(StringIO(b))
        alignment = next(iterator)
        self.assertEqual([r.id for r in alignment], ["b.fa/0-4", "a.fa/0-4"])
        self.assertEqual(iterator.ids, ["2", "1"])

    def test_seq_count(self):
        """A block whose record count differs from seq_count is an error."""
        data = "> 1:1-4 + a.fa\nACGT\n> 2:1-4 + b.fa\nACGA\n=\n"
        self.assertEqual(len(next(MauveIterator(StringIO(data), seq_count=2))), 2)
        with self.assertRaisesRegex(
            ValueError, "^Found 2 records in this alignment, told to expect 3$"
        ):
            next(MauveIterator(StringIO(data), seq_count=3))


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
