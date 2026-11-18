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
            mo_type = bg_type = float
            for position in range(pssm.length):
                shifts = []
                mo_weights = []
                bg_weights = []
                lo = pssm[:, position]
                for letter, score in lo.items():
                    bg = background[letter]
                    mo = pow(2, pssm[letter, position]) * bg
                    shifts.append(self._index_diff(score))
                    mo_weights.append(mo)
                    bg_weights.append(bg)
                mo_density, mo_type = self._add_column(
                    mo_density, mo_type, shifts, mo_weights
                )
                bg_density, bg_type = self._add_column(
                    bg_density, bg_type, shifts, bg_weights
                )
            self.mo_density = self._as_list(mo_density, mo_type, shifts)
            self.bg_density = self._as_list(bg_density, bg_type, shifts)

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

    def _add_column(self, density, value_type, shifts, weights):
        """Return the density after one more column, and its value type.

        The new density sums, over the column's letters, the density times
        the letter's weight, moved up by the letter's shift (see _shift).
        value_type is the type of the density's values, Python float at the
        start. The per-bin loop summed ``density[i] * weight`` in the type
        that Python and numpy give that expression, so float32 weights gave
        float32 sums; the new density is summed in that type too, and so
        gives the same thresholds.
        """
        total = 0.0
        for weight in weights:
            total += value_type(1) * weight
        value_type = type(total)
        dtype = np.dtype(value_type)
        density = density.astype(dtype, copy=False)
        new = np.zeros(self.n_points, dtype)
        for d, weight in zip(shifts, weights):
            self._shift(new, density, d, dtype.type(weight))
        return new, value_type

    def _as_list(self, density, value_type, shifts):
        """Return the density as a list of values of value_type.

        The threshold methods add these up one by one, so numpy scalars,
        such as float32 ones, keep those sums in their own type. As in the
        per-bin loop, the bins that none of the last column's shifts reach
        hold the Python 0.0 they started as; for a rate too small for
        value_type, a sum still at Python 0.0 compares differently.
        """
        if value_type is float:
            return density.tolist()
        values = list(density)
        # Each shift d reaches bins 0 to n - 1 + d if d <= 0, or bins d to
        # n - 1 if d >= 0, so the bins reached by none of them lie between.
        n = self.n_points
        shifts = [max(1 - n, min(n - 1, d)) for d in shifts]
        start = max((n + d for d in shifts if d <= 0), default=0)
        stop = min((d for d in shifts if d >= 0), default=n)
        values[start:stop] = [0.0] * (stop - start)
        return values

    @staticmethod
    def _value_type(density):
        """Return the type of the density's values.

        At most one end of the density holds bins left as Python 0.0 (see
        _as_list), so this is the type of its first or last value.
        """
        first, last = type(density[0]), type(density[-1])
        return last if first is float else first

    def modify(self, scores, mo_probs, bg_probs):
        """Modify motifs and background density."""
        shifts = [self._index_diff(v) for v in scores.values()]
        mo_density, mo_type = self._add_column(
            np.array(self.mo_density),
            self._value_type(self.mo_density),
            shifts,
            [mo_probs[k] for k in scores],
        )
        bg_density, bg_type = self._add_column(
            np.array(self.bg_density),
            self._value_type(self.bg_density),
            shifts,
            [bg_probs[k] for k in scores],
        )
        self.mo_density = self._as_list(mo_density, mo_type, shifts)
        self.bg_density = self._as_list(bg_density, bg_type, shifts)

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
