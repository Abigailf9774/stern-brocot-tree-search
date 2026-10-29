"""stern-brocot-tree-search: best rational approximation within a denominator bound.

The Stern-Brocot tree is a binary search tree over the positive rationals.
Every reduced fraction p/q appears exactly once; an in-order traversal yields
them in ascending order. This module navigates that tree to find the fraction
p/q with 1 <= q <= max_denominator that minimises |target - p/q|.

Public API:
    best_rational_approximation(target, max_denominator)
        -> Approximation (with .numerator, .denominator, .value, .error)
    Approximation
"""

from .core import Approximation, best_rational_approximation

__all__ = ["Approximation", "best_rational_approximation"]
