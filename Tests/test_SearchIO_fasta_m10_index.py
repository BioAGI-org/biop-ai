# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO fasta-m10 indexing."""

import os
import unittest
from unittest import mock

import support
from search_tests_common import CheckIndex
from search_tests_common import CheckRaw
from search_tests_common import EndOfFileGuard
from search_tests_common import run_python

from Bio import SearchIO


def lines(filename, start, stop):
    """Return lines start + 1 to stop of a test file, as bytes."""
    with open(support.DATA / "Fasta" / filename, "rb") as handle:
        return b"".join(handle.readlines()[start:stop])


class FastaM10RawCases(CheckRaw):
    """Check FASTA m10 get_raw method."""

    fmt = "fasta-m10"

    def test_multi_block_bgzf(self):
        """Test get_raw on a BGZF copy with many blocks, as on the plain file."""
        self.check_raw_multi_block_bgzf(support.DATA / "Fasta" / "output002.m10")

    def test_output_001_many_hits(self):
        """Test fasta-m10 raw string retrieval, fasta35, query with many hits (output001.m10)."""
        filename = support.DATA / "Fasta" / "output001.m10"
        # The query's first line in its table of best scores, 40000 times
        self.check_raw_repeated(
            filename,
            "gi|10955263|ref|NP_052604.1|",
            b"gi|152973457|ref|YP_001338508.1| ATPase with chape (",
            b"\n",
            40000,
        )


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
        # replaced by spaces, the old regular expression took time growing
        # with about the fourth power of the number of spaces to reject the
        # line, as in test_SearchIO_fasta_m10.py: with 2400 spaces, well over
        # 20 minutes. The indexer matches the line as bytes, so time that too
        # in a child process, with the same 60 second timeout.
        filename = os.fspath(support.DATA / "Fasta" / "output001.m10")
        code = f"""\
from io import BytesIO
from unittest import mock
from Bio import SearchIO
with open({filename!r}, "rb") as handle:
    lines = handle.readlines()
line = lines[7].replace(b" - 107 aa\\n", b" " * 2400 + b"\\n")
with mock.patch(
    "Bio.SearchIO._index._open_for_random_access",
    return_value=BytesIO(b"".join(lines[:7]) + line),
):
    try:
        SearchIO.index("test.m10", "fasta-m10")
    except ValueError:
        print("ValueError")
"""
        result = run_python(code)
        self.assertEqual(result.stdout, "ValueError\n", result.stderr)
        # With 600 spaces, in this process, the error quotes the line.
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
