# Numerical Verification of Extreme Zero Bounds for Meixner Polynomials

This repository contains a Python implementation for computing the zeros of monic Meixner polynomials and testing several analytical bounds for their smallest and largest zeros.

The code generates numerical tables comparing the exact extreme zeros with bounds obtained from different mixed recurrence relations.

## Features

- Construction of monic Meixner polynomials through the three-term recurrence relation.
- Computation of polynomial zeros using NumPy.
- Implementation of seven families of bounds:
  
  - F₀
  - F₁
  - F₂
  - F₃
  - F₄
  - F₅
  - F₆

- Numerical comparison between:
  
  - smallest zero \(x_{1,n}\),
  - largest zero \(x_{n,n}\),
  - lower bounds \(F_n^{-}(k)\),
  - upper bounds \(F_n^{+}(k)\).

- Reproduction of numerical tables used in the study of Meixner polynomial zeros.

---

## Mathematical Background

The monic Meixner polynomials satisfy the three-term recurrence relation

\[
M_{n+1}(x)
=
(x-C_n)M_n(x)
-\lambda_n M_{n-1}(x),
\]

where

\[
C_n=\frac{c(\beta+n-1)+(n-1)}{1-c},
\]

and

\[
\lambda_n=
\frac{c(n-1)(\beta+n-2)}
     {(1-c)^2}.
\]

The implemented bounds are derived from mixed recurrence relations and interlacing properties of orthogonal polynomials.

---

## Requirements

Install the required packages:

```bash
pip install numpy sympy
```

---

## Usage

Run

```bash
python meixner_bounds.py
```

The program outputs a table of the form

```text
------------------------------------------------------------------------------------
    beta      c   k         x_{1,n}       F_n^-(k)       F_n^+(k)         x_{n,n}
------------------------------------------------------------------------------------
...
```

where

- `x_{1,n}` is the smallest zero,
- `x_{n,n}` is the largest zero,
- `F_n^-(k)` is the lower bound,
- `F_n^+(k)` is the upper bound.

---

## Example Parameter Sets

The current implementation evaluates

- \(n=8\),
- \(\beta=0.09\), \(c=0.02\),
- \(\beta=0.09\), \(c=0.50\),
- \(\beta=0.09\), \(c=0.99\),
- \(\beta=20\), \(c=0.50\),
- \(\beta=20\), \(c=0.99\).

Additional parameter values can easily be added by modifying the `param_sets` list.

---

## Repository Structure

```text
.
├── meixner_bounds.py
├── README.md
└── LICENSE
```

---

## Applications

This code can be used for

- numerical verification of theoretical bounds,
- research on orthogonal polynomials,
- studies of zero distributions,
- validation of recurrence-based estimates,
- reproducibility of numerical experiments appearing in research articles.

---

## Citation

If you use this code in academic work, please cite the corresponding article or preprint where the bounds are derived.
