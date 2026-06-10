import numpy as np
from itertools import product
from numpy.polynomial import polynomial as P
import sympy as sp

# ===Pochhammer symbol===
def poch(a, n):
    result = 1.0
    for k in range(n):
        result *= (a + k)
    return result

# ===Meixner polynomial coefficients via 3-term recurrence===
# M_{n+1}(x) = (x - C_n)*M_n(x) - Lambda_n*M_{n-1}(x)
def cn(n, beta, c):
    return (c*(beta + n - 1) + (n - 1)) / (1 - c)

def lambda_n(n, beta, c):
    return c*(n - 1)*(beta + n - 2) / (1 - c)**2

def meixner_coeffs(n, beta, c):
    """Build monic Meixner polynomial via recurrence, return numpy coeffs."""
    if n == 0:
        return np.array([1.0])
    p_prev = np.array([1.0])           # M_0 = 1
    p_curr = np.array([-cn(1, beta, c), 1.0])  # M_1 = x - C_1
    for k in range(1, n):
        C  = cn(k + 1, beta, c)
        Lm = lambda_n(k + 1, beta, c)
        # M_{k+1} = (x - C)*M_k - Lm*M_{k-1}
        p_next = np.polymul([-C, 1.0], p_curr)
        p_next[: len(p_prev)] -= Lm * p_prev
        p_prev = p_curr
        p_curr = p_next
    return p_curr

def meixner_zeros(n, beta, c):
    coeffs = meixner_coeffs(n, beta, c)
    roots = np.roots(coeffs[::-1])
    roots = np.sort(roots.real)
    return roots

# ===Bounds F_k (each returns [upper, lower])===
def Fn0(n, b, c):
    num  = (2*n - 3)*(1 + c) + 2*c*b
    disc = (1 + c)**2 + 4*c*(n - 1)*(b + n - 2)
    return [(num + np.sqrt(disc)) / (2*(1 - c)),
            (num - np.sqrt(disc)) / (2*(1 - c))]

def A1(n,b,c): return (c-1)**3 / (c**2*(n-2)*(n-1)*(n+b-2))
def B1(n,b,c): return (c-1)**2*(2*n-3+c*(n+2*b-1)) / (c**2*(n-2)*(n-1)*(n+b-2))
def D1(n,b,c): return (c-1)*(n**2-3*n+2+c*(-2+c*(b-1))*b+n*(c*(1+c)*b)) / (c**2*(n-2)*(n-1)*(n+b-2))
def Fn1(n,b,c):
    A,B,D = A1(n,b,c), B1(n,b,c), D1(n,b,c)
    disc = B**2 - 4*A*D
    return [(-B + np.sqrt(disc))/(2*A), (-B - np.sqrt(disc))/(2*A)]

def A2(n,b,c): return -(c-1)**2 / (c**2*(n-2)*(n-1))
def B2(n,b,c): return -(c-1)*(2*n+c+2*c*b-3) / (c**2*(n-2)*(n-1))
def D2(n,b,c): return ((c-1)*(n-2)*(n-1) - c*(n+c-2)*b - c**2*b**2) / (c**2*(n-2)*(n-1))
def Fn2(n,b,c):
    A,B,D = A2(n,b,c), B2(n,b,c), D2(n,b,c)
    disc = B**2 - 4*A*D
    return [(-B + np.sqrt(disc))/(2*A), (-B - np.sqrt(disc))/(2*A)]

def A3(n,b,c): return (c-1)*(n+b-1) / (c**2*(n-2)*(n-1))
def B3(n,b,c): return -(n+b-1)*(-2*n+3+c*(n-2*b-3)) / (c**2*(n-2)*(n-1))
def D3(n,b,c): return (n+b-1)*((c-1)**2*(n-2)*(n-1)+c*(n-c*(n-3)-2)*b+c**2*b**2) / ((c-1)*c**2*(n-2)*(n-1))
def Fn3(n,b,c):
    A,B,D = A3(n,b,c), B3(n,b,c), D3(n,b,c)
    disc = B**2 - 4*A*D
    return [(-B + np.sqrt(disc))/(2*A), (-B - np.sqrt(disc))/(2*A)]

def N4(n,b,c): return ((c-1)**3*n**3 - (c-1)**2*n**2*(3*c+b-3)
    + (c-1)*n*(2*c**2+c*(b**2+4*b-4)-3*b+2)
    - b*(c*(b+3)*(c*b+2*c-2)+2))
def A4(n,b,c): return -(n+b-1)*(n+b) / (c**2*(n-2)*(n-1))
def B4(n,b,c): return (n+b-1)*(n+b)*(-2*n+3+c*(2*n-2*b-5)) / ((c-1)*c**2*(n-2)*(n-1))
def D4(n,b,c): return (n+b-1)*N4(n,b,c) / ((c-1)**2*c**2*(n-2)*(n-1))
def Fn4(n,b,c):
    A,B,D = A4(n,b,c), B4(n,b,c), D4(n,b,c)
    disc = B**2 - 4*A*D
    return [(-B + np.sqrt(disc))/(2*A), (-B - np.sqrt(disc))/(2*A)]

def N5B(n,b,c): return ((2-3*c+c**3)*n**3
    + n**2*(-1-3*c**3-4*c*(-1+b)+4*b)
    + b*(1+b)*(-3+c*(7+2*b))
    + n*(-3+2*c**3+2*(-2+b)*b+c*(7+b*(13+b))))
def N5D(n,b,c): return ((-1+c)**4*n**4 + (-1+c)**3*n**3*(2+c*(-2+b)-2*b)
    - (-1+c)**2*n**2*((-1+c)**2+(5+3*(-3+c)*c)*b-b**2)
    + (-1+c)*n*(2*(-1+c)**3+(-1+c*(5+2*(-4+c)*c))*b+3*(1-2*c)*b**2-c*b**3)
    + b*(1+b)*(2+c*(4+b)*(-2+c*(3+b))))
def A5(n,b,c): return (n+b-1)*(n+b)*(n+b+1) / ((c-1)*c**2*(n-2)*(n-1))
def B5(n,b,c): return (n+b-1)*N5B(n,b,c) / ((c-1)**2*c**2*(n-2)*(n-1))
def D5(n,b,c): return (n+b-1)*N5D(n,b,c) / ((c-1)**3*c**2*(n-2)*(n-1))
def Fn5(n,b,c):
    A,B,D = A5(n,b,c), B5(n,b,c), D5(n,b,c)
    disc = B**2 - 4*A*D
    return [(-B + np.sqrt(disc))/(2*A), (-B - np.sqrt(disc))/(2*A)]

def N6A(n,b,c): return ((c**3-1)*n**3 - 3*(c**3+b+1)*n**2
    + (2*c**3-3*b*(b+2)-2)*n - b*(b+1)*(b+2))
def N6B(n,b,c): return (2*(c-1)**3*(c+1)*n**4
    + (c-1)**2*(-3*(c**2+c+1)+2*(c*(c-1)-3)*b)*n**3
    - (-5+19*c+c**3*(5*c-19)+3*b+9*c*b+6*c**3*(c-3)*b-6*(c-1)*b**2)*n**2
    + (6-b*(b+2)*(2*b-7)-6*c**3*(2*b+3)+c**4*(4*b+6)-c*(18+b*(2*b**2+27*b+50)))*n
    - b*(b+1)*(b+2)*(c*(2*b+9)-3))
def N6D(n,b,c): return ((-1+c)**5*n**5 + (-1+c)**4*(-3+2*c)*n**4*b
    + (-1+c)**3*n**3*(-5*(-1+c)**2+(-3+(7-3*c)*c)*b+(3+(-3+c)*c)*b**2)
    + (-1+c)*n*(4*(-1+c)**4+2*(3+(-1+c)*c*(14+c*(-11+3*c)))*b
        +(-3+2*c*(1+c*(10+(-5+c)*c)))*b**2+3*(-1+3*c)*b**3+c*b**4)
    - (-1+c)**2*n**2*b*(c*(26+c*(-19+5*c))+3*c*(6+(-4+c)*c)*b+b**2-2*(5+3*b))
    - b*(1+b)*(2+b)*(2+c*(5+b)*(-2+c*(4+b))))
def A6(n,b,c): return (n+b-1)*N6A(n,b,c) / ((c-1)**2*c**2*(n-2)*(n-1))
def B6(n,b,c): return (n+b-1)*N6B(n,b,c) / ((c-1)**3*c**2*(n-2)*(n-1))
def D6(n,b,c): return (n+b-1)*N6D(n,b,c) / ((c-1)**4*c**2*(n-2)*(n-1))
def Fn6(n,b,c):
    A,B,D = A6(n,b,c), B6(n,b,c), D6(n,b,c)
    disc = B**2 - 4*A*D
    return [(-B + np.sqrt(disc))/(2*A), (-B - np.sqrt(disc))/(2*A)]

def get_bounds(k, n, b, c):
    funcs = [Fn0, Fn1, Fn2, Fn3, Fn4, Fn5, Fn6]
    return funcs[k](n, b, c)

# ===Print table matching the LaTeX output===
param_sets = [
    (8, 0.09,  0.02),
    (8, 0.09,  0.50),
    (8, 0.09,  0.99),
    (8, 20.0,  0.50),
    (8, 20.0,  0.99),
]

sep = "-" * 84
print(sep)
print(f"{'beta':>8} {'c':>6} {'k':>3}  {'x_{1,n}':>14} {'F_n^-(k)':>14} {'F_n^+(k)':>14} {'x_{n,n}':>14}")
print(sep)

for (n, beta, c) in param_sets:
    zeros = meixner_zeros(n, beta, c)
    x1n   = zeros[0]
    xnn   = zeros[-1]
    first = True
    for k in range(7):
        vals  = get_bounds(k, n, beta, c)
        Fplus = max(vals)
        Fminus= min(vals)
        x1_str = f"{x1n:14.6f}" if first else f"{'':14}"
        xn_str = f"{xnn:14.6f}" if first else f"{'':14}"
        b_str  = f"{beta:>8.4f}" if first else f"{'':>8}"
        c_str  = f"{c:>6.4f}"   if first else f"{'':>6}"
        print(f"{b_str} {c_str} {k:>3}  {x1_str} {Fminus:14.6f} {Fplus:14.6f} {xn_str}")
        first = False
    print(sep)


