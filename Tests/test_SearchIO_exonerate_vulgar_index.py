# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO exonerate-vulgar indexing."""

import unittest

import support
from search_tests_common import CheckIndex
from search_tests_common import CheckRaw


class ExonerateVulgarRawCases(CheckRaw):
    """Check exonerate-vulgar get_raw method."""

    fmt = "exonerate-vulgar"

    def test_exn_22_o_vulgar_many_hits(self):
        """Test exonerate-vulgar raw string retrieval, query with many hits (exn_22_o_vulgar.exn)."""
        filename = support.DATA / "Exonerate" / "exn_22_o_vulgar.exn"
        # The query's second vulgar line, 40000 times
        self.check_raw_repeated(
            filename, "sacCer3_dna", b"vulgar: sacCer3_dna 0 529 +", b"\n", 40000
        )


class ExonerateVulgarIndexCases(CheckIndex):
    fmt = "exonerate-vulgar"

    def test_exn_22_m_est2genome(self):
        """Test exonerate-vulgar indexing, single."""
        filename = support.DATA / "Exonerate" / "exn_22_o_vulgar.exn"
        self.check_index(filename, self.fmt)

    def test_exn_22_q_multiple(self):
        """Test exonerate-vulgar indexing, single."""
        filename = support.DATA / "Exonerate" / "exn_22_q_multiple_vulgar.exn"
        self.check_index(filename, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
