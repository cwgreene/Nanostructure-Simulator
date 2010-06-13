"""
drift diffusion

Primary equations are:
div grad(V) - (h(x) - n(x) ) = 0
h(x) - exp(V) = 0
n(x) - exp(V) = 0

The exp(V) is a nonlinear term, so we need to 
focus on the first (linear) term, and use a dummy function
to handle the exp(V).
"""

from dolfin import *
V=FunctionSpace(mesh,"Langrange",2)

v = TrialFunction(V)
a = grad(
