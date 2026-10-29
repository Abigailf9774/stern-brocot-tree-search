# stern-brocot-tree-search

Find the best rational approximation to a real number within a given denominator bound using Stern-Brocot tree navigation.

## Usage

```python
from stern_brocot_tree_search import best_rational_approximation, Approximation

r = best_rational_approximation(3.14159265358979, 100)
print(r.numerator, r.denominator)  # 311 99
print(r.value)                      # 311/99
print(r.error)                      # absolute error as a Fraction
```

`best_rational_approximation(target: float, max_denominator: int) -> Approximation` returns the reduced fraction `p/q` with `1 <= q <= max_denominator` that minimises `|target - p/q|`. The result is a frozen dataclass with four fields: `numerator` (int), `denominator` (int), `value` (Fraction), and `error` (Fraction, the absolute difference between target and value).

## Why this exists

When you need a fraction with a bounded denominator that is as close as possible to a given real number, the Stern-Brocot tree gives a clean O(d) descent where d is the denominator of the result. The alternatives — enumerating every Farey fraction up to the bound, or rolling your own continued-fraction bookkeeping — are either slower or more error-prone.

The trade-off: the algorithm is exact in exact arithmetic but the input is a Python `float`, so the first step converts the float to a `Fraction` via `limit_denominator(10**12)`. This means the library approximates the real number that the float actually represents, which for values like `0.1` is the binary float nearest to 0.1, not 1/10 exactly. The returned `error` reflects this.

## Edge cases

- `max_denominator` must be >= 1; smaller values raise `ValueError`.
- `NaN`, `+inf`, and `-inf` targets raise `ValueError` — the tree assumes a finite, ordered input.
- Negative targets are supported: the search runs on `|target|` and the sign is restored.
- Integers return `n/1` with zero error.
- The output fraction is always in lowest terms.
