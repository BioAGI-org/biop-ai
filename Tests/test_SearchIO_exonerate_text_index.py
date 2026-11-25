# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for SearchIO exonerate-text indexing."""

import unittest

import support
from search_tests_common import CheckIndex
from search_tests_common import CheckRaw


class ExonerateTextRawCases(CheckRaw):
    """Check exonerate-text get_raw method."""

    fmt = "exonerate-text"

    def test_multi_block_bgzf(self):
        """Test get_raw on a BGZF copy with many blocks, as on the plain file."""
        self.check_raw_multi_block_bgzf(
            support.DATA / "Exonerate" / "exn_22_q_multiple.exn"
        )

    def test_exn_22_m_est2genome_long_alignment(self):
        """Test exonerate-text raw string retrieval, query with a long alignment (exn_22_m_est2genome.exn)."""
        filename = support.DATA / "Exonerate" / "exn_22_m_est2genome.exn"
        # The alignment's second block, 4 lines, 10000 times
        self.check_raw_repeated(
            filename,
            "gi|296143771|ref|NM_001180731.1|",
            b"      57 : GCAAATGTTG",
            b"\n\n",
            10000,
        )


class ExonerateTextIndexCases(CheckIndex):
    fmt = "exonerate-text"

    def test_exn_22_m_est2genome(self):
        """Test exonerate-text indexing, single."""
        filename = support.DATA / "Exonerate" / "exn_22_m_est2genome.exn"
        self.check_index(filename, self.fmt)

    def test_exn_22_q_multiple(self):
        """Test exonerate-text indexing, single."""
        filename = support.DATA / "Exonerate" / "exn_22_q_multiple.exn"
        self.check_index(filename, self.fmt)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
