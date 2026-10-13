# HC-ts-diagonal

Code and data supporting the note *Linear subspaces and the Hodge classes of very general hypersurfaces
F(x)+F(y)=0* (Kevin Saller, 2026), DOI [10.5281/zenodo.23329495](https://doi.org/10.5281/zenodo.23329495).

**Subject.** For a form F of degree d in v variables, the note studies the hypersurface
X_F = {F(x) + F(y) = 0} in P^{2v-1}, of dimension 2k = 2v - 2. It contains the d linear spaces
L_η = {y = η x}, η^d = -1. For very general F the note bounds the rational Hodge classes of degree 2k by a
first-order (infinitesimal Noether–Lefschetz) computation in the Jacobian ring, which by Thom–Sebastiani splits into
chains of tensor products of graded pieces of the Jacobian ring of F, and compares the bound with the classes of
h^k and the L_η. It also treats Noether–Lefschetz subfamilies (F vanishing on a fixed subvariety), two-block
hypersurfaces F(x) + G(y) = 0 and the integral lattice spanned by h^k and the L_η.

**Role of the computations.** The proofs in the note are by hand. The computations in this repository are checks
of statements of the note (kernel dimensions modulo primes at random members, exact Gram determinants, a lattice
index); none of them is used in a proof. [`README.txt`](README.txt) gives the requirements, the commands, what
each script checks, and the result of the rerun made when this repository was assembled.

**Status:** not yet peer reviewed.

**Disclosure:** the arguments and computations were produced with the assistance of an AI system (Claude, by
Anthropic). The computations were checked with separately written code (a Python and a Singular implementation of
the first-order kernels) and the arguments in separate internal reviews, which were also carried out with AI
assistance. What was not checked is listed in README.txt.

## Licences

- **Code** (`code/`): MIT License, see [`LICENSE`](LICENSE).
- **Everything else** (README files and computed outputs in `outputs/`): Creative Commons Attribution 4.0
  International (CC BY 4.0), see [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).
- No third-party files are included.
