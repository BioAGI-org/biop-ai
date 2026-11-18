# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO fasta-m10 indexing."""

import unittest
from io import BytesIO
from unittest import mock

import support
from search_tests_common import CheckIndex

from Bio import SearchIO


class EndOfFileGuard(BytesIO):
    """Binary handle that fails if the indexer keeps reading at end of file.

    An indexer stuck in a loop at end of file would otherwise hang the suite.
    """

    def __init__(self, data):
        """Initialize the handle with the bytes to read."""
        super().__init__(data)
        self.reads_at_end = 0

    def readline(self, *args):
        """Read a line, failing after many reads at end of file."""
        line = super().readline(*args)
        if not line:
            self.reads_at_end += 1
            if self.reads_at_end > 100:
                raise AssertionError("indexer kept reading at end of file")
        return line


def lines(filename, start, stop):
    """Return lines start + 1 to stop of a test file, as bytes."""
    with open(support.DATA / "Fasta" / filename, "rb") as handle:
        return b"".join(handle.readlines()[start:stop])


class FastaM10IndexCases(CheckIndex):
    fmt = "fasta-m10"

    def test_output_002(self):
        """Test fasta-m10 indexing, fasta34, multiple queries."""
        filename = support.DATA / "Fasta" / "output002.m10"
        self.check_index(filename, self.fmt)

    def test_output_001(self):
        """Test fasta-m10 indexing, fasta35, multiple queries."""
        filename = support.DATA / "Fasta" / "output001.m10"
        self.check_index(filename, self.fmt)

    def test_output_005(self):
        """Test fasta-m10 indexing, ssearch35, multiple queries."""
        filename = support.DATA / "Fasta" / "output005.m10"
        self.check_index(filename, self.fmt)

    def test_output_008(self):
        """Test fasta-m10 indexing, tfastx36, multiple queries."""
        filename = support.DATA / "Fasta" / "output008.m10"
        self.check_index(filename, self.fmt)

    def test_output_009(self):
        """Test fasta-m10 indexing, fasta36, multiple queries."""
        filename = support.DATA / "Fasta" / "output009.m10"
        self.check_index(filename, self.fmt)

    def test_output_010(self):
        """Test fasta-m10 indexing, fasta36, single query, no hits."""
        filename = support.DATA / "Fasta" / "output010.m10"
        self.check_index(filename, self.fmt)

    def test_output_011(self):
        """Test fasta-m10 indexing, fasta36, single query, hits with single hsp."""
        filename = support.DATA / "Fasta" / "output011.m10"
        self.check_index(filename, self.fmt)

    def test_output_012(self):
        """Test fasta-m10 indexing, fasta36, single query with multiple hsps."""
        filename = support.DATA / "Fasta" / "output012.m10"
        self.check_index(filename, self.fmt)


class TruncatedCases(unittest.TestCase):
    """A file with no query, or none after its preamble, raises ValueError."""

    def index(self, data):
        """Index the bytes as fasta-m10, through a handle guarding end of file."""
        with mock.patch(
            "Bio.SearchIO._index._open_for_random_access",
            return_value=EndOfFileGuard(data),
        ):
            index = SearchIO.index("test.m10", "fasta-m10")
        self.addCleanup(index.close)
        return index

    def test_empty_file(self):
        """Test an empty file, which holds no query results."""
        self.assertEqual(len(self.index(b"")), 0)

    def test_no_query(self):
        """Test a file that ends before its first query line."""
        # Lines 1 to 7 of output001.m10 are its preamble, ending with
        # "Query: NC_002127.faa"; line 8 is the first query line, "  1>>>gi|...".
        # A FASTA sequence file, given under the wrong format name, has no
        # query line at all.
        for name, data in (
            ("output001.m10, 3 lines", lines("output001.m10", 0, 3)),
            ("output001.m10, 7 lines", lines("output001.m10", 0, 7)),
            ("FASTA sequence", b">seq1\nACGT\n"),
        ):
            with self.subTest(name):
                with self.assertRaises(ValueError) as cm:
                    self.index(data)
                self.assertEqual(
                    str(cm.exception),
                    "File ended before the first query line of the FASTA output",
                )

    def test_query_without_preamble(self):
        """Test reading back a lone query that has no preamble before it."""
        # Lines 8 to 105 of output001.m10 are its first query. get_raw reads
        # the preamble from the start of the file to the first query line
        # after line 1, and here there is none.
        index = self.index(lines("output001.m10", 7, 105))
        key = "gi|10955263|ref|NP_052604.1|"
        self.assertEqual(list(index), [key])
        for name, get in (("get_raw", index.get_raw), ("getitem", index.__getitem__)):
            with self.subTest(name):
                with self.assertRaises(ValueError) as cm:
                    get(key)
                self.assertEqual(
                    str(cm.exception),
                    "File ended before the first query line after the FASTA"
                    " preamble",
                )


class LongLineCases(unittest.TestCase):
    """Long query lines are indexed in linear time."""

    def test_query_line_without_length(self):
        """Test a query line without its '- <length> aa' ending."""
        # Line 8 of output001.m10 is the first query line, "  1>>>gi|...
        # [Escherichia coli O157:H7 s 107 aa - 107 aa". With its "- 107 aa"
        # replaced by spaces, the old regular expression took about 7 seconds
        # to reject it, then the indexer failed with AttributeError.
        line = lines("output001.m10", 7, 8).replace(b" - 107 aa\n", b" " * 600 + b"\n")
        data = lines("output001.m10", 0, 7) + line
        with mock.patch(
            "Bio.SearchIO._index._open_for_random_access",
            return_value=EndOfFileGuard(data),
        ):
            with self.assertRaises(ValueError) as cm:
                SearchIO.index("test.m10", "fasta-m10")
        self.assertEqual(
            str(cm.exception),
            "Expected a query line '>>>ID DESCRIPTION - LENGTH aa' (or nt),"
            f" found:\n{line!r}",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
