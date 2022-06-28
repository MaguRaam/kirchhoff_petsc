from sympy import *

r, t, c0, f0, t0 = symbols('r t c0 f0 t0 pi_')

init_printing(use_unicode=True)


# monopole source:
source = -2.0*(t - t0)*f0**2*exp( -1.0*f0**2*(t - t0)**2)

# derivative of source
print (diff (source, t))


