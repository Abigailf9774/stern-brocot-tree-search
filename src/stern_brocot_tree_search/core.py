"""Core Stern-Brocot tree search for best rational approximation."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class Approximation:
    """A rational approximation result.

    Fields:
        numerator:   integer p
        denominator: integer q (always >= 1; reduced with numerator)
        value:       Fraction(p, q)
        error:       |target - p/q| as a Fraction (>= 0)
    """

    numerator: int
    denominator: int
    value: Fraction
    error: Fraction


def best_rational_approximation(target: float, max_denominator: int) -> Approximation:
    """Find the best rational approximation to *target*.

    Returns the reduced fraction p/q with 1 <= q <= max_denominator that
    minimises |target - p/q|.

    Algorithm: iterative Stern-Brocot descent between adjacent Farey
    neighbours. We maintain two bounding fractions left and right (initialised
    to 0/1 and 1/0). At each step we form their mediant m = (a+c)/(b+d) and
    compare it with the target:
        - if m == target exactly, we are done;
        - if target < m, move the right bound to m;
        - if target > m, move the left bound to m.
    After each move we consider the newly formed mediant as a candidate. We
    stop when the mediant denominator exceeds max_denominator; the best among
    all candidates seen is returned.

    Interpretation chosen (stated plainly): *max_denominator* is the inclusive
    upper bound on q. Raising it past the point where the denominator stops
    changing means there is no better fraction within the bound, and we return
    the same result.

    Edge cases handled:
        - target is negative, zero, or an exact integer.
        - target is infinity or NaN (ValueError).
        - max_denominator < 1 (ValueError).
    """
    if max_denominator < 1:
        raise ValueError("max_denominator must be >= 1")

    # Reject non-finite targets early: the algorithm assumes an ordering.
    if target != target:  # NaN check (more portable than math.isnan for float)
        raise ValueError("target must be a finite real number, not NaN")
    if target in (float("inf"), float("-inf")):
        raise ValueError("target must be a finite real number, not infinity")

    t = Fraction(target).limit_denominator(10**12)

    # Exact small denominator check for integers.
    if t.denominator == 1:
        return Approximation(t.numerator, 1, t, Fraction(0, 1))

    a_lo, a_hi = 0, 1  # numerators
    b_lo, b_hi = 1, 0  # denominators  (0/1 and 1/0 enclose everything positive)

    # If t < 0 we work in |t| and restore the sign at the end.
    negate = t < 0
    tt = -t if negate else t

    best = None  # (error, num, den, val)

    # Walk down the tree: narrow [lo, hi] around tt until the next mediant
    # would exceed the denominator bound.
    while True:
        med_num = a_lo + a_hi
        med_den = b_lo + b_hi
        if med_den > max_denominator:
            break
        med = Fraction(med_num, med_den)
        err = abs(tt - med)
        if best is None or err < best[0]:
            best = (err, med_num, med_den, med)
        if med == tt:
            break
        if med < tt:
            a_lo, b_lo = med_num, med_den
        else:
            a_hi, b_hi = med_num, med_den

    # Also consider the surviving bounds as candidates.
    lo = Fraction(a_lo, b_lo)
    hi = Fraction(a_hi, b_hi)
    err_lo = abs(tt - lo)
    err_hi = abs(tt - hi)
    if best is None or err_lo < best[0]:
        best = (err_lo, a_lo, b_lo, lo)
    if best is None or err_hi < best[0]:
        best = (err_hi, a_hi, b_hi, hi)

    _, chosen_num, chosen_den, chosen_val = best
    return _build(chosen_num, chosen_den, chosen_val, tt, negate)


def _build(num: int, den: int, val: Fraction, tt: Fraction, negate: bool) -> Approximation:
    """Assemble an Approximation, restoring the sign if the target was negative."""
    if negate:
        num = -num
        val = Fraction(num, den)
    err = abs(tt - val) if not negate else abs(-tt - val)
    return Approximation(num, den, val, err)
