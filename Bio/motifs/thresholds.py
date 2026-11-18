# Copyright 2008 by Norbert Dojer.  All rights reserved.
# Adapted by Bartek Wilczynski.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Approximate calculation of appropriate thresholds for motif finding."""

import numpy as np


class ScoreDistribution:
    """Class representing approximate score distribution for a given motif.

    Utilizes a dynamic programming approach to calculate the distribution of
    scores with a predefined precision. Provides a number of methods for calculating
    thresholds for motif occurrences.
    """

    def __init__(self, motif=None, precision=10**3, pssm=None, background=None):
        """Initialize the class."""
        if pssm is None:
            self.min_score = min(0.0, motif.min_score())
            self.interval = max(0.0, motif.max_score()) - self.min_score
            self.n_points = precision * motif.length
            self.ic = motif.ic()
        else:
            self.min_score = min(0.0, pssm.min)
            self.interval = max(0.0, pssm.max) - self.min_score
            self.n_points = precision * pssm.length
            self.ic = pssm.mean(background)
        self.step = self.interval / (self.n_points - 1)
        self.mo_density = [0.0] * self.n_points
        self.mo_density[-self._index_diff(self.min_score)] = 1.0
        self.bg_density = [0.0] * self.n_points
        self.bg_density[-self._index_diff(self.min_score)] = 1.0
        if pssm is None:
            for lo, mo in zip(motif.log_odds(), motif.pwm()):
                self.modify(lo, mo, motif.background)
        else:
            mo_density = np.array(self.mo_density)
            bg_density = np.array(self.bg_density)
            for position in range(pssm.length):
                mo_new = np.zeros(self.n_points)
                bg_new = np.zeros(self.n_points)
                lo = pssm[:, position]
                for letter, score in lo.items():
                    bg = background[letter]
                    mo = pow(2, pssm[letter, position]) * bg
                    d = self._index_diff(score)
                    self._shift(mo_new, mo_density, d, mo)
                    self._shift(bg_new, bg_density, d, bg)
                mo_density = mo_new
                bg_density = bg_new
            self.mo_density = mo_density.tolist()
            self.bg_density = bg_density.tolist()

    def _index_diff(self, x, y=0.0):
        return int((x - y + 0.5 * self.step) // self.step)

    def _shift(self, new, density, d, weight):
        """Add density times weight to new, moved up by d bins.

        This is ``new[max(0, min(n - 1, i + d))] += density[i] * weight`` for
        each bin i in turn, done by numpy. Bins moved past either end pile up
        in the end bin, added in the same order as by that loop, so that the
        sums are the same to the last bit.
        """
        n = self.n_points
        # A shift of n - 1 bins already moves every bin to the end bin.
        d = max(1 - n, min(n - 1, d))
        if d >= 0:
            new[d:] += density[: n - d] * weight
            np.add.at(new, np.full(d, n - 1), density[n - d :] * weight)
        else:
            np.add.at(new, np.zeros(-d, np.intp), density[:-d] * weight)
            new[: n + d] += density[-d:] * weight

    def modify(self, scores, mo_probs, bg_probs):
        """Modify motifs and background density."""
        mo_density = np.array(self.mo_density)
        bg_density = np.array(self.bg_density)
        mo_new = np.zeros(self.n_points)
        bg_new = np.zeros(self.n_points)
        for k, v in scores.items():
            d = self._index_diff(v)
            self._shift(mo_new, mo_density, d, mo_probs[k])
            self._shift(bg_new, bg_density, d, bg_probs[k])
        self.mo_density = mo_new.tolist()
        self.bg_density = bg_new.tolist()

    def threshold_fpr(self, fpr):
        """Approximate the log-odds threshold which makes the type I error (false positive rate)."""
        i = self.n_points
        prob = 0.0
        while prob < fpr:
            i -= 1
            prob += self.bg_density[i]
        return self.min_score + i * self.step

    def threshold_fnr(self, fnr):
        """Approximate the log-odds threshold which makes the type II error (false negative rate)."""
        i = -1
        prob = 0.0
        while prob < fnr:
            i += 1
            prob += self.mo_density[i]
        return self.min_score + i * self.step

    def threshold_balanced(self, rate_proportion=1.0, return_rate=False):
        """Approximate log-odds threshold making FNR equal to FPR times rate_proportion."""
        i = self.n_points
        fpr = 0.0
        fnr = 1.0
        while fpr * rate_proportion < fnr:
            i -= 1
            fpr += self.bg_density[i]
            fnr -= self.mo_density[i]
        if return_rate:
            return self.min_score + i * self.step, fpr
        else:
            return self.min_score + i * self.step

    def threshold_patser(self):
        """Threshold selection mimicking the behaviour of patser (Hertz, Stormo 1999) software.

        It selects such a threshold that the log(fpr)=-ic(M)
        note: the actual patser software uses natural logarithms instead of log_2, so the numbers
        are not directly comparable.
        """
        return self.threshold_fpr(fpr=2**-self.ic)
