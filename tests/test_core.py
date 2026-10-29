"""Tests for stern_brocot_tree_search.core."""

import math
import unittest
from fractions import Fraction

from stern_brocot_tree_search import Approximation, best_rational_approximation


class TestBestRationalApproximation(unittest.TestCase):

    def test_pi_max_10(self):
        # 22/7 is the classic best approximation of pi for q <= 10.
        r = best_rational_approximation(3.14159265358979, 10)
        self.assertEqual(r.numerator, 22)
        self.assertEqual(r.denominator, 7)
        self.assertEqual(r.value, Fraction(22, 7))

    def test_pi_max_100(self):
        # For q <= 100 the best is 311/99.
        r = best_rational_approximation(3.14159265358979, 100)
        self.assertEqual(r.numerator, 311)
        self.assertEqual(r.denominator, 99)

    def test_zero(self):
        r = best_rational_approximation(0.0, 10)
        self.assertEqual(r.numerator, 0)
        self.assertEqual(r.denominator, 1)
        self.assertEqual(r.error, Fraction(0, 1))

    def test_integer(self):
        r = best_rational_approximation(5.0, 10)
        self.assertEqual(r.numerator, 5)
        self.assertEqual(r.denominator, 1)
        self.assertEqual(r.error, Fraction(0, 1))

    def test_one_half(self):
        r = best_rational_approximation(0.5, 10)
        self.assertEqual(r.numerator, 1)
        self.assertEqual(r.denominator, 2)
        self.assertEqual(r.error, Fraction(0, 1))

    def test_negative(self):
        r = best_rational_approximation(-3.14159265358979, 10)
        self.assertEqual(r.numerator, -22)
        self.assertEqual(r.denominator, 7)
        self.assertEqual(r.value, Fraction(-22, 7))

    def test_negative_small(self):
        r = best_rational_approximation(-0.3333, 10)
        self.assertEqual(r.numerator, -1)
        self.assertEqual(r.denominator, 3)

    def test_denominator_one(self):
        # target close to an integer with max_denominator=1
        r = best_rational_approximation(2.7, 1)
        self.assertEqual(r.denominator, 1)
        self.assertIn(r.numerator, (2, 3))

    def test_reduced_output(self):
        r = best_rational_approximation(0.5, 10)
        self.assertEqual(math.gcd(r.numerator, r.denominator), 1)

    def test_max_denominator_one_integer(self):
        r = best_rational_approximation(7.0, 1)
        self.assertEqual(r.numerator, 7)
        self.assertEqual(r.denominator, 1)

    def test_golden_ratio(self):
        # phi ~ 1.6180339887; 13/8 is a convergent for small bounds.
        r = best_rational_approximation(1.61803398874989, 10)
        self.assertEqual(r.numerator, 13)
        self.assertEqual(r.denominator, 8)

    def test_sqrt2(self):
        # sqrt(2) ~ 1.41421356; 24/17 is best for q <= 20.
        r = best_rational_approximation(1.4142135623730951, 20)
        self.assertEqual(r.numerator, 24)
        self.assertEqual(r.denominator, 17)

    def test_error_field(self):
        r = best_rational_approximation(0.1, 10)
        # 0.1's float value equals Fraction(1, 10) after limit_denominator,
        # so the best fraction for q<=10 is exactly 1/10 with zero error.
        self.assertEqual(r.numerator, 1)
        self.assertEqual(r.denominator, 10)
        self.assertEqual(r.error, Fraction(0))

    def test_nan_raises(self):
        with self.assertRaises(ValueError):
            best_rational_approximation(float("nan"), 10)

    def test_inf_raises(self):
        with self.assertRaises(ValueError):
            best_rational_approximation(float("inf"), 10)

    def test_neg_inf_raises(self):
        with self.assertRaises(ValueError):
            best_rational_approximation(float("-inf"), 10)

    def test_bad_max_denominator(self):
        with self.assertRaises(ValueError):
            best_rational_approximation(1.5, 0)

    def test_large_target(self):
        # target ~ 100.5 should give 201/2.
        r = best_rational_approximation(100.5, 10)
        self.assertEqual(r.numerator, 201)
        self.assertEqual(r.denominator, 2)

    def test_just_above_one(self):
        r = best_rational_approximation(1.0001, 1000)
        self.assertEqual(r.denominator, 1)
        self.assertEqual(r.numerator, 1)

    def test_monotone_improvement(self):
        # Larger denominator bound should not increase the error.
        r1 = best_rational_approximation(3.14159265358979, 5)
        r2 = best_rational_approximation(3.14159265358979, 50)
        self.assertLessEqual(r2.error, r1.error)


if __name__ == "__main__":
    unittest.main()
