# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO blat-psl indexing."""

import unittest

import support
from search_tests_common import CheckIndex
from search_tests_common import CheckRaw


class BlatPslRawCases(CheckRaw):
    """Check BLAT PSL get_raw method."""

    fmt = "blat-psl"

    def test_psl_34_001_many_hits(self):
        """Test blat-psl raw string retrieval, query with many hits (psl_34_001.psl)."""
        filename = support.DATA / "Blat" / "psl_34_001.psl"
        # The query's first hit line, 40000 times
        self.check_raw_repeated(
            filename,
            "hg19_dna",
            b"38\t3\t0\t0\t0\t0\t0\t0\t+\thg19_dna\t",
            b"\n",
            40000,
        )


class BlatPslIndexCases(CheckIndex):
    fmt = "blat-psl"

    def test_psl_34_001(self):
        """Test blat-psl indexing, multiple queries."""
        filename = support.DATA / "Blat" / "psl_34_001.psl"
        self.check_index(filename, self.fmt)

    def test_psl_34_002(self):
        """Test blat-psl indexing, single query, no hits."""
        filename = support.DATA / "Blat" / "psl_34_002.psl"
        self.check_index(filename, self.fmt)

    def test_psl_34_003(self):
        """Test blat-psl indexing, single query, single hit."""
        filename = support.DATA / "Blat" / "psl_34_003.psl"
        self.check_index(filename, self.fmt)

    def test_psl_34_004(self):
        """Test blat-psl indexing, single query, multiple hits with multiple hsps."""
        filename = support.DATA / "Blat" / "psl_34_004.psl"
        self.check_index(filename, self.fmt)

    def test_psl_34_005(self):
        """Test blat-psl indexing, multiple queries, no header."""
        filename = support.DATA / "Blat" / "psl_34_005.psl"
        self.check_index(filename, self.fmt)

    def test_psl_34_006(self):
        """Test blat-pslx indexing, multiple queries."""
        filename = support.DATA / "Blat" / "pslx_34_001.pslx"
        self.check_index(filename, self.fmt, pslx=True)

    def test_psl_34_007(self):
        """Test blat-pslx indexing, single query, no hits."""
        filename = support.DATA / "Blat" / "pslx_34_002.pslx"
        self.check_index(filename, self.fmt, pslx=True)

    def test_psl_34_008(self):
        """Test blat-pslx indexing, single query, single hit."""
        filename = support.DATA / "Blat" / "pslx_34_003.pslx"
        self.check_index(filename, self.fmt, pslx=True)

    def test_psl_34_009(self):
        """Test blat-pslx indexing, single query, multiple hits with multiple hsps."""
        filename = support.DATA / "Blat" / "pslx_34_004.pslx"
        self.check_index(filename, self.fmt, pslx=True)

    def test_psl_34_010(self):
        """Test blat-pslx indexing, multiple queries, no header."""
        filename = support.DATA / "Blat" / "pslx_34_005.pslx"
        self.check_index(filename, self.fmt, pslx=True)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
