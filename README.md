# 1-D-Schrodinger
The script is an example of solving eigenvalue problem via the Secant and Numerov method.
The system under study is one-dimensional Schroedinger equation.

a brief description of code : 
(1) Choose the region of the numerical solution. This region should be large enough
compared with the effective region of the potential to have a negligible effect on the
solution.
(2) Provide a reasonable guess for the lowest eigenvalue ε0. This can be found approximately
from the analytical result of the case with an infinite well and the same range
of well width.
4.9 The one-dimensional Schro¨ dinger equation 107
(3) Integrate the equation for φl(x) from the left to the point xr + h and the one for φr(x)
from the right to xr − h.We can choose zero to be the value of the first points of φl(x)
and φr(x), and a small quantity to be the value of the second points of φl(x) and φr(x),
to start the integration, for example, with the Numerov algorithm. Before matching
the solutions, rescale one of them to ensure that φl(xr) = φr(xr). For example, we can
multiply φl(x) by φr(xr)/φl(xr) up to x = xr + h. This rescaling also ensures that the
solutions have the correct nodal structure, that is, changing the sign of φl(x) if it is
incorrect.
(4) Evaluate f (ε0) = [φr(xr − h) − φr(xr + h) − φl(xr − h) + φl(xr + h)]/2hφr(xr).
(5) Apply a root search method to obtain ε0 from f (ε0) = 0 within a given tolerance.
(6) Carry out the above steps for the next eigenvalue. We can start the search with a
slightly higher value than the last eigenvalue. We need to make sure that no eigenstate
is missed. This can easily be done by counting the nodes in the solution; the
nth state has a total number of n nonboundary nodes, with n = 0 for the ground
state. A node is where φ(x) = 0. This also provides a way of pinpointing a specific
eigenstate.
