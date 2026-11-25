# Copyright 2008 by Norbert Dojer.  All rights reserved.
# Adapted by Bartek Wilczynski.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Approximate calculation of appropriate thresholds for motif finding."""

import contextlib
import math
import warnings

import numpy as np

from Bio import BiopythonWarning


def _cast(values, dtype):
    """Return values as dtype, rounding a value too small for it to zero.

    numpy converts a Python float operand to a narrower type, such as
    float32, without checking for underflow, whatever numpy.seterr says.
    """
    with np.errstate(under="ignore"):
        return values.astype(dtype, copy=False)


def _where(kinds, present):
    """Yield each kind in kinds, with where it is: everywhere, or a mask.

    present holds every kind in kinds, and maybe others.
    """
    if len(present) == 1:
        yield next(iter(present)), ...
    else:
        for kind in np.flatnonzero(np.bincount(kinds)):
            yield kind, kinds == kind


class _Types:
    """The types of a density's values, numbered as kinds.

    The per-bin loop computed ``new[j] += density[i] * weight`` with Python
    floats and numpy scalars, so each product and sum has the type that
    Python and numpy give it, and is rounded in that type. With weights of
    different types, such as float32 and float64, those types differ from
    letter to letter and from bin to bin, and bins that no letter reaches
    keep the Python 0.0 they started as.
    """

    def __init__(self, types):
        self.types = []
        self.index = {}
        for value_type in types:
            self.kind(value_type)
        self._sums = None

    def kind(self, value_type):
        """Return the kind of value_type, numbering it if it is new."""
        if value_type not in self.index:
            self.index[value_type] = len(self.types)
            self.types.append(value_type)
        return self.index[value_type]

    def dtype(self, kind):
        """Return the dtype of the given kind."""
        return np.dtype(self.types[kind])

    def wide(self):
        """Return a dtype that holds a value of each kind exactly."""
        return np.result_type(np.float64, *map(np.dtype, self.types))

    def products(self, weight):
        """Return the kind of value * weight for a value of each kind."""
        types = list(self.types)
        return np.array([self.kind(type(t(1) * weight)) for t in types], np.int8)

    def sums(self):
        """Return the kind of a + b for an a and a b of each kind."""
        while self._sums is None or len(self._sums) < len(self.types):
            types = list(self.types)
            self._sums = np.array(
                [[self.kind(type(a(0) + b(0))) for b in types] for a in types],
                np.int8,
            )
        return self._sums

    def errstate(self, kind):
        """Return the numpy error state for a product of the given kind.

        Python floats round a product too small for a float to zero,
        whatever numpy.seterr says; numpy scalars follow numpy.seterr.
        """
        if self.types[kind] is float:
            return np.errstate(under="ignore")
        return contextlib.nullcontext()


class _Bins:
    """The values of a density's bins, each kept in its own type.

    values holds them in a dtype wide enough for all of them, kinds the
    kind of each one's type in types, and present each kind in kinds (it
    may hold others).
    """

    def __init__(self, values, kinds, types, present):
        self.values = values
        self.kinds = kinds
        self.types = types
        self.present = present

    @classmethod
    def from_list(cls, density):
        """Return the bins of a list of Python floats or numpy scalars."""
        types = _Types(dict.fromkeys(map(type, density)))
        kinds = np.array([types.index[type(value)] for value in density], np.int8)
        values = np.array(density, types.wide())
        return cls(values, kinds, types, set(range(len(types.types))))

    def tolist(self):
        """Return the values as a list, each value of its own type."""
        density = [0.0] * len(self.values)
        for kind, where in _where(self.kinds, self.present):
            value_type = self.types.types[kind]
            values = self.values[where].astype(value_type)
            values = values.tolist() if value_type is float else list(values)
            if where is ...:
                return values
            for i, value in zip(np.flatnonzero(where).tolist(), values):
                density[i] = value
        return density


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
            for letter in pssm.alphabet:
                if not all(map(math.isfinite, pssm[letter])):
                    raise ValueError(
                        "Cannot compute the score distribution of a matrix with "
                        "an infinite or NaN value; a log-odds matrix has -inf "
                        "wherever a count is zero, which pseudocounts avoid, and "
                        "+inf or NaN wherever a background frequency is zero"
                    )
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
            mo_density = _Bins.from_list(self.mo_density)
            bg_density = _Bins.from_list(self.bg_density)
            mismatched = False
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
                mo_density = self._add_column(mo_density, shifts, mo_weights)
                bg_density = self._add_column(bg_density, shifts, bg_weights)
                if not math.isclose(math.fsum(mo_weights), 1, rel_tol=1e-3):
                    mismatched = True
            self.mo_density = mo_density.tolist()
            self.bg_density = bg_density.tolist()
            if mismatched:
                # 2**score * background gives the motif's probabilities only
                # for the background the PSSM was built with.
                warnings.warn(
                    "The background is not the one the PSSM was built with, so "
                    "the motif score density does not add up to 1, and "
                    "threshold_fnr and threshold_balanced will be wrong.",
                    BiopythonWarning,
                )

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

    def _add_column_in(self, new, density, shifts, weights, kind):
        """Add the column to new, with every product and sum in one kind."""
        types = density.types
        dtype = types.dtype(kind)
        values = np.zeros(self.n_points, dtype)
        density = _cast(density.values, dtype)
        with types.errstate(kind):
            for d, weight in zip(shifts, weights):
                self._shift(values, density, d, dtype.type(weight))
        new.values = values.astype(types.wide(), copy=False)
        if types.types[kind] is not float:
            # Each shift d reaches bins 0 to n - 1 + d if d <= 0, or bins d
            # to n - 1 if d >= 0, so the bins reached by none of them lie
            # between.
            n = self.n_points
            shifts = [max(1 - n, min(n - 1, d)) for d in shifts]
            start = max((n + d for d in shifts if d <= 0), default=0)
            stop = min((d for d in shifts if d >= 0), default=n)
            new.kinds[:start] = kind
            new.kinds[stop:] = kind
            new.present.add(kind)

    def _add_column(self, density, shifts, weights):
        """Return the density after one more column.

        This is the per-bin loop, done by numpy: for each letter in turn,
        with its shift d and weight, and for each bin i in turn,
        ``new[max(0, min(n - 1, i + d))] += density[i] * weight``. Each
        product and sum is done in the type that Python and numpy give it,
        which can differ from letter to letter and from bin to bin, so the
        values and their types are the same to the last bit. Bins moved
        past either end pile up in the end bin, added in the loop's order.
        """
        types = density.types
        start = types.kind(float)
        tables = [types.products(weight) for weight in weights]
        sums = types.sums()
        wide = types.wide()
        density.values = density.values.astype(wide, copy=False)
        n = self.n_points
        new = _Bins(None, np.full(n, start, np.int8), types, {start})
        column_kinds = {table[kind] for table in tables for kind in density.present}
        if len(column_kinds) == 1 and column_kinds == {
            sums[start, kind] for kind in column_kinds
        }:
            # One type, as with a background of one type, for every product
            # and sum; but bins that no letter reaches stay Python 0.0.
            self._add_column_in(new, density, shifts, weights, column_kinds.pop())
            return new
        new.values = np.zeros(n, wide)
        for d, weight, table in zip(shifts, weights, tables):
            # Each bin times the weight, in the type of that product
            present = {table[kind] for kind in density.present}
            if len(present) == 1:
                kinds = next(iter(present))
            else:
                kinds = table[density.kinds]
            values = np.empty(n, wide)
            for kind, where in _where(kinds, present):
                dtype = types.dtype(kind)
                with types.errstate(kind):
                    product = _cast(density.values[where], dtype) * dtype.type(weight)
                if where is ...:
                    values = product.astype(wide, copy=False)
                else:
                    values[where] = product
            # Bin i goes to bin i + d, or to the end bin past either end.
            d = max(1 - n, min(n - 1, d))
            if d >= 0:
                target, source = slice(d, n - 1), slice(0, n - 1 - d)
                end, pile = n - 1, slice(n - 1 - d, n)
            else:
                target, source = slice(1, n + d), slice(1 - d, n)
                end, pile = 0, slice(0, 1 - d)
            sum_present = {sums[a, b] for a in new.present for b in present}
            if len(sum_present) == 1:
                new_kinds = next(iter(sum_present))
            else:
                new_kinds = sums[new.kinds[target], np.broadcast_to(kinds, n)[source]]
                sum_present = set(np.flatnonzero(np.bincount(new_kinds)))
            new_values = new.values[target]
            for kind, where in _where(new_kinds, sum_present):
                # A sum cannot underflow; this only lets numpy round a
                # Python float to a narrower type unchecked, as it does.
                dtype = types.dtype(kind)
                with np.errstate(under="ignore"):
                    if where is ...:
                        np.add(new_values, values[source], new_values, dtype=dtype)
                    else:
                        new_values[where] = np.add(
                            new_values[where], values[source][where], dtype=dtype
                        )
            if len(sum_present) > 1 or new.present != sum_present:
                new.kinds[target] = new_kinds
                new.present |= sum_present
            if len(present) > 1:
                kinds = kinds[pile]
            new.values[end], new.kinds[end] = self._add_in_order(
                types, new.values[end], new.kinds[end], values[pile], kinds
            )
            new.present.add(new.kinds[end])
        return new

    @staticmethod
    def _add_in_order(types, total, kind, values, kinds):
        """Return total plus each of values in turn, and the sum's kind.

        kinds holds the kind of each value, or is the kind of all of them.
        The sum's type changes only where a value's type widens it, so
        each run of values added in one type is added by np.add.at.
        """
        sums = types.sums()
        start = 0
        while start < len(values):
            if np.ndim(kinds):
                kind = sums[kind, kinds[start]]
                changes = np.flatnonzero(sums[kind, kinds[start + 1 :]] != kind)
                stop = start + 1 + changes[0] if len(changes) else len(values)
            else:
                kind = sums[kind, kinds]
                stop = len(values)
            dtype = types.dtype(kind)
            total = _cast(np.array([total]), dtype)
            indices = np.zeros(stop - start, np.intp)
            np.add.at(total, indices, _cast(values[start:stop], dtype))
            total = total[0]
            start = stop
        return total, kind

    def modify(self, scores, mo_probs, bg_probs):
        """Modify motifs and background density."""
        shifts = [self._index_diff(v) for v in scores.values()]
        mo_density = self._add_column(
            _Bins.from_list(self.mo_density),
            shifts,
            [mo_probs[k] for k in scores],
        )
        bg_density = self._add_column(
            _Bins.from_list(self.bg_density),
            shifts,
            [bg_probs[k] for k in scores],
        )
        self.mo_density = mo_density.tolist()
        self.bg_density = bg_density.tolist()

    def threshold_fpr(self, fpr):
        """Approximate the log-odds threshold which makes the type I error (false positive rate)."""
        if not 0 <= fpr <= 1:
            raise ValueError(f"fpr must be between 0 and 1, not {fpr}")
        i = self.n_points
        prob = 0.0
        # Rounding can make the density add up to just under 1.
        while prob < fpr and i > 0:
            i -= 1
            prob += self.bg_density[i]
        return self.min_score + i * self.step

    def threshold_fnr(self, fnr):
        """Approximate the log-odds threshold which makes the type II error (false negative rate)."""
        if not 0 <= fnr <= 1:
            raise ValueError(f"fnr must be between 0 and 1, not {fnr}")
        i = -1
        prob = 0.0
        # Rounding can make the density add up to just under 1.
        while prob < fnr and i < self.n_points - 1:
            i += 1
            prob += self.mo_density[i]
        return self.min_score + i * self.step

    def threshold_balanced(self, rate_proportion=1.0, return_rate=False):
        """Approximate log-odds threshold making FNR equal to FPR times rate_proportion."""
        if not 0 < rate_proportion < float("inf"):
            raise ValueError(
                "rate_proportion must be a finite number greater than 0, "
                f"not {rate_proportion}"
            )
        i = self.n_points
        fpr = 0.0
        fnr = 1.0
        # Rounding can leave fnr just above 0 at the lowest bin.
        while fpr * rate_proportion < fnr and i > 0:
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
