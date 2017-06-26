Blaze-Iterative
===============

This is a set of iterative linear system solvers intended for inclusion 
into the Blaze library, a high-performance C++ linear algebra library
(https://bitbucket.org/blaze-lib/).
The API is currently based on a tag-dispatch system to choose a particular algorithm.

Currently implemented algorithms:
- Conjugate Gradient (CG)
- BiCGSTAB


Planned algorithms:
- Preconditioned CG
- Preconditioned BiCGSTAB
- (Preconditioned) BiCGSTAB(l)
- GMRES
- LSQR
- LSMR

Please open an issue to discuss/request features.
