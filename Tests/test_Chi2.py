# Copyright (C) 2011 by Brandon Invergo (b.invergo@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Chi2 module."""

import unittest
from unittest import mock

from Bio.Phylo.PAML import chi2


class ModTest(unittest.TestCase):
    def test_cdf_chi2(self):
        self.assertRaises(ValueError, chi2.cdf_chi2, df=0, stat=3.84)
        self.assertRaises(ValueError, chi2.cdf_chi2, df=1, stat=-3.84)
        self.assertRaises(TypeError, chi2.cdf_chi2, df="d", stat="stat")
        self.assertAlmostEqual(chi2.cdf_chi2(2, 3.84), 0.1466070, places=5)

    def test_ln_gamma(self):
        self.assertRaises(ValueError, chi2._ln_gamma_function, -1)
        self.assertAlmostEqual(chi2._ln_gamma_function(10), 12.80183, places=5)

    def test_incomplete_gamma(self):
        self.assertRaises(ValueError, chi2._incomplete_gamma, x=0.5, alpha=-1)
        self.assertAlmostEqual(chi2._incomplete_gamma(0.5, 0.5), 0.6826895, places=5)


class TestContinuedFraction(unittest.TestCase):
    """The continued fraction in _incomplete_gamma ends on every input."""

    def setUp(self):
        # Before the fix the continued fraction often looped forever when df
        # was not 2 or 4, and always for statistics so large that its terms
        # overflowed to nan. It calls abs twice a round; stop such a loop
        # early so that the test fails instead.
        def bounded_abs(value):
            self.abs_calls += 1
            if self.abs_calls > 10000:
                raise AssertionError("_incomplete_gamma called abs 10000 times")
            return abs(value)

        patcher = mock.patch.object(chi2, "abs", bounded_abs, create=True)
        patcher.start()
        self.addCleanup(patcher.stop)

    def cdf_chi2(self, df, stat):
        self.abs_calls = 0
        return chi2.cdf_chi2(df, stat)

    def test_critical_values(self):
        """The 5%, 1% and 0.1% critical values give p-values to 7 places."""
        critical_values = {
            1: (3.841458820694124, 6.634896601021214, 10.827566170662733),
            2: (5.991464547107979, 9.21034037197618, 13.815510557964274),
            3: (7.814727903251178, 11.344866730144373, 16.26623619623813),
            4: (9.487729036781154, 13.276704135987622, 18.46682695290317),
            5: (11.070497693516351, 15.08627246938899, 20.515005652432873),
        }
        for df, stats in critical_values.items():
            for stat, p_value in zip(stats, (0.05, 0.01, 0.001)):
                with self.subTest(df=df, stat=stat):
                    self.assertAlmostEqual(self.cdf_chi2(df, stat), p_value, places=7)

    def test_huge_statistic(self):
        """A statistic too large for the continued fraction gives p = 0."""
        for df, stat in ((1, 1e102), (1, 1e105), (2, 1e104), (1, 1e300)):
            with self.subTest(df=df, stat=stat):
                self.assertEqual(self.cdf_chi2(df, stat), 0.0)
        for df in (1, 2, 3):
            with self.subTest(df=df, stat="inf"):
                self.assertEqual(self.cdf_chi2(df, float("inf")), 0.0)

    def test_overflow(self):
        """An overflow in the continued fraction raises ValueError."""
        with self.assertRaises(ValueError) as cm:
            self.cdf_chi2(2e200, 2e200)
        self.assertEqual(
            str(cm.exception),
            "The continued fraction for the incomplete gamma ratio "
            "overflowed at x=1e+200, alpha=1e+200",
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
