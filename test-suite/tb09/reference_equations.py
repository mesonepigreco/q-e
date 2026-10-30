#!/usr/bin/env python3
"""Independent BR/TB09 reference for XClib/test_tb09.f90 (Hartree/bohr).

Tran and Blaha, PRL 102, 226401 (2009), DOI:10.1103/PhysRevLett.102.226401;
Becke and Roussel, PRA 39, 3761 (1989), DOI:10.1103/PhysRevA.39.3761.
Uses bisection of the BR equation; imports no QE or Libxc implementation.
Only regular positive-density fixtures are evaluated here.
"""
import math


def spin_potential(n, grad_squared, laplacian, tau, c):
    """Single-spin density and tau=(1/2) sum f_i |grad psi_i|^2."""
    d = 2*tau - grad_squared/(4*n)
    q = (laplacian - 2*0.8*d)/6
    if abs(q) < 1.e-18:
        x = 2.0
    else:
        rhs = (2/3)*math.pi**(2/3)*n**(5/3)/q
        # Solve x exp(-2x/3) = rhs (x-2), using the physical root.
        lo, hi = (0., 2.) if q < 0 else (2., max(4., 2.+2./rhs))
        def f(x):
            return x*math.exp(-2*x/3) - rhs*(x-2)
        flo = f(lo)
        for _ in range(200):
            mid = (lo+hi)/2
            if f(mid)*flo > 0:
                lo, flo = mid, f(mid)
            else:
                hi = mid
        x = (lo+hi)/2
    b = (x**3*math.exp(-x)/(8*math.pi*n))**(1/3)
    vbr = -(1-math.exp(-x)*(1+x/2))/b
    return c*vbr + (3*c-2)/math.pi*math.sqrt(5/12)*math.sqrt(2*tau/n)


if __name__ == "__main__":
    # The unpolarized Libxc inputs split equally between the two spins.
    for sigma, lapl in [(0.0025, 0.03), (0.0025, -0.03), (0., 0.), (0., 0.03), (0.0025, 1.)]:
        print(f"{spin_potential(.05, sigma/4, lapl/2, .04, 1.2):.15f}")
    print("Unequal spin channels:")
    for n, sigma, lapl, tau in [(.06, .001, .015, .05), (.04, .0005, -.01, .03)]:
        print(f"{spin_potential(n, sigma, lapl, tau, 1.2):.15f}")
