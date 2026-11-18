# Copyright 2006-2014 by Peter Cock.  All rights reserved.
# Revisions copyright 2011 Brandon Invergo. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.AlignIO.PhylipIO module."""

import unittest
from io import StringIO

import support

from Bio.AlignIO.PhylipIO import PhylipIterator
from Bio.AlignIO.PhylipIO import PhylipWriter
from Bio.AlignIO.PhylipIO import RelaxedPhylipIterator
from Bio.AlignIO.PhylipIO import SequentialPhylipIterator


class TestPhylipIO(unittest.TestCase):
    def test_one(self):
        input_file = support.DATA / "Phylip" / "one.dat"
        with open(input_file) as handle:
            ids = []
            for alignment in PhylipIterator(handle):
                for record in alignment:
                    ids.append(record.id)
        self.assertEqual(
            ids,
            [
                "V_Harveyi_",
                "B_subtilis",
                "B_subtilis",
                "YA80_HAEIN",
                "FLIY_ECOLI",
                "E_coli_Gln",
                "Deinococcu",
                "HISJ_E_COL",
            ],
        )

        expected = (
            """mkklvlslsl vlafssataa faaipqniri gtdptyapfe sknsqgelvg
        fdidlakelc krintqctfv enpldalips lkakkidaim sslsitekrq qeiaftdkly
        aadsrlvvak nsdiqptves lkgkrvgvlq gttqetfgne hwapkgieiv syqgqdniys
        dltagridaafqdevaaseg flkqpvgkdy kfggpsvkde klfgvgtgmg lrkednelre
        alnkafaemradgtyeklak kyfdfdvygg""".replace(
                " ", ""
            )
            .replace("\n", "")
            .upper()
        )
        self.assertEqual(str(record.seq).replace("-", ""), expected)

    def test_two_and_three(self):
        path = support.DATA / "Phylip" / "two.dat"
        # derived from http://atgc.lirmm.fr/phyml/usersguide.html
        with open(path) as handle:
            list2 = list(PhylipIterator(handle))
        self.assertEqual(len(list2), 1)
        self.assertEqual(len(list2[0]), 5)

        path = support.DATA / "Phylip" / "three.dat"
        with open(path) as handle:
            list3 = list(PhylipIterator(handle))
        self.assertEqual(len(list3), 1)
        self.assertEqual(len(list3[0]), 5)

        for i in range(5):
            self.assertEqual(list2[0][i].id, list3[0][i].id)
            self.assertEqual(list2[0][i].seq, list3[0][i].seq)

    def test_four(self):
        path = support.DATA / "Phylip" / "four.dat"
        # File derived from here:
        # http://evolution.genetics.washington.edu/phylip/doc/sequence.html
        # Note the lack of any white space between names 2 and 3 and their seqs.
        with open(path) as handle:
            list4 = list(PhylipIterator(handle))
        self.assertEqual(len(list4), 1)
        self.assertEqual(len(list4[0]), 5)

    def test_five(self):
        # File derived rom here:
        # http://evolution.genetics.washington.edu/phylip/doc/sequence.html
        path = support.DATA / "Phylip" / "five.dat"
        with open(path) as handle:
            self.assertRaises(ValueError, list, PhylipIterator(handle))

    def test_six(self):
        # File derived rom here:
        # http://evolution.genetics.washington.edu/phylip/doc/sequence.html
        path = support.DATA / "Phylip" / "six.dat"
        with open(path) as handle:
            list5 = list(PhylipIterator(handle))
        self.assertEqual(len(list5), 1)

    def test_concatenation(self):
        path = support.DATA / "Phylip" / "one.dat"
        with open(path) as handle:
            phylip_text = handle.read()
        path = support.DATA / "Phylip" / "three.dat"
        with open(path) as handle:
            phylip_text3 = handle.read()
        path = support.DATA / "Phylip" / "four.dat"
        with open(path) as handle:
            phylip_text4 = handle.read()
        handle = StringIO(phylip_text4 + "\n" + phylip_text4)
        self.assertEqual(len(list(PhylipIterator(handle))), 2)

        handle = StringIO(phylip_text3 + "\n" + phylip_text4 + "\n\n\n" + phylip_text)
        self.assertEqual(len(list(PhylipIterator(handle))), 3)

    def test_write_read(self):
        path = support.DATA / "Phylip" / "six.dat"
        with open(path) as handle:
            list5 = list(PhylipIterator(handle))

        handle = StringIO()
        PhylipWriter(handle).write_file(list5)
        handle.seek(0)
        list6 = list(PhylipIterator(handle))

        self.assertEqual(len(list5), len(list6))
        for a1, a2 in zip(list5, list6):
            self.assertEqual(len(a1), len(a2))
            for r1, r2 in zip(a1, a2):
                self.assertEqual(r1.id, r2.id)
                self.assertEqual(r1.seq, r2.seq)


class EndOfFileGuard(StringIO):
    """A text handle that fails the test if read on and on at its end."""

    def __init__(self, text):
        """Hold the text, and count the reads made at its end."""
        super().__init__(text)
        self.reads_at_end = 0

    def readline(self, size=-1):
        """Return the next line, failing after 100 reads at the end."""
        line = super().readline(size)
        if not line:
            self.reads_at_end += 1
            if self.reads_at_end > 100:
                raise AssertionError("readline called 100 times at end of file")
        return line


class TestPhylipEndOfFile(unittest.TestCase):
    """Files ending before the sequences their header announces."""

    iterators = (PhylipIterator, RelaxedPhylipIterator, SequentialPhylipIterator)

    def check(self, iterator, text, message):
        with self.assertRaises(ValueError) as cm:
            list(iterator(EndOfFileGuard(text)))
        self.assertEqual(str(cm.exception), message)

    def read_text(self, filename):
        with open(support.DATA / "Phylip" / filename) as handle:
            return handle.read()

    def test_header_only(self):
        # Reading stops at the end of the file, not a billion lines later.
        for iterator in self.iterators:
            with self.subTest(iterator=iterator.__name__):
                self.check(
                    iterator,
                    "1000000000 10\n",
                    "End of file before sequence 1 of 1000000000",
                )

    def test_header_only_of_second_alignment(self):
        for iterator, filename in (
            (PhylipIterator, "four.dat"),
            (RelaxedPhylipIterator, "interlaced.phy"),
            (SequentialPhylipIterator, "sequential.phy"),
        ):
            with self.subTest(iterator=iterator.__name__):
                self.check(
                    iterator,
                    self.read_text(filename) + "\n1000000000 10\n",
                    "End of file before sequence 1 of 1000000000",
                )

    def test_truncated(self):
        # four.dat has 5 sequences on 5 lines in its first block,
        # interlaced.phy 3 sequences on 3 lines, and sequential.phy
        # 3 sequences on 7 lines each.
        for iterator, filename, lines, message in (
            (PhylipIterator, "four.dat", 3, "sequence 4 of 5"),
            (RelaxedPhylipIterator, "interlaced.phy", 2, "sequence 3 of 3"),
            (SequentialPhylipIterator, "sequential.phy", 7, "sequence 2 of 3"),
            (SequentialPhylipIterator, "sequential.phy", 10, "sequence 3 of 3"),
        ):
            with self.subTest(iterator=iterator.__name__, lines=lines):
                text = self.read_text(filename).splitlines(keepends=True)
                self.check(
                    iterator,
                    "".join(text[: lines + 1]),
                    "End of file before " + message,
                )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
