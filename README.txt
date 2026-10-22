Supporting code and data for the note
"Linear subspaces and the Hodge classes of very general hypersurfaces F(x)+F(y)=0"
==================================================================================

Status: accompanies the note (K. Saller, 2026), DOI 10.5281/zenodo.23329495
(https://doi.org/10.5281/zenodo.23329495); not peer reviewed.

Setting. F is a form of degree d in v variables x = (x_1..x_v), and X_F = {F(x) + F(y) = 0} in P^{2v-1}, of
dimension 2k = 2v - 2. R = R_F is the Jacobian ring of F, of socle degree sigma = v(d - 2). By Thom-Sebastiani
the Jacobian ring of F(x) + F(y) is R (x) R, and the primitive (k,k)-part of H^{2k}(X_F) is its piece of degree
N = sigma, the sum of R_a (x) R_{N-a}. A tangent direction A(x) + A(y) (A of degree d) maps it to degree N + d.
The first-order kernel K(F) is the set of elements killed by all tangent directions of the family considered;
dim K(F) bounds the dimension of the primitive Hodge classes of the very general member from above. K(F) splits
into "chains" (residue classes of a modulo d). X_F contains the d linear spaces L_eta = {y = eta x}, eta^d = -1.
eps_j = (1 - (1 - d)^j)/d; L_eta . L_eta = eps_v.

Role of the computations. The proofs in the note are by hand. The computations below are CHECKS of statements of
the note; none of them is used in a proof. They compute kernel dimensions exactly over F_p (python-flint or
Singular) at random members (or at the Fermat point), Gram determinants exactly over Q, and a lattice index
exactly over Z.

Not checked by any script: the arguments of the note are proofs, not computations. In particular the step from
a kernel dimension computed modulo p at one member to the bound for the very general complex member (lower
semicontinuity of the rank, reduction modulo p), the identification of the first-order kernel with an upper
bound for the Hodge classes, the intersection numbers used as input (L_eta . L_eta = eps_v and the other entries
of the Gram matrices; the intersection formula for the d^v linear spaces at the Fermat point), and the
algebraicity of the classes are arguments of the note and its references.

Labels. The script headers use working labels of the author's notes (for example S3, E1, E1*, E4, U1, U2, J3,
P1, "Theorem A*", "Lemma K", "NOTES"). The notes are not part of this repository, and the labels are not
numbers of the note.

Licence: the code in code/ is under the MIT licence; the outputs in outputs/ and this text are under CC BY 4.0
(see the LICENSE files of the repository). No third-party files are included.


1. Requirements
---------------
Python 3 with python-flint (all scripts except gram_check.py and sat_check.py) and numpy (sat_check.py only);
gram_check.py uses the standard library only. Tested with Python 3.12.9, python-flint 0.9.0 and numpy 2.1.3.
The three .sing scripts need Singular; they were run with Singular 4.4.1 as shipped in the docker image
sagemath/sagemath (SageMath 10.10). Run times on a desktop PC: under 2 minutes per command, except
ivhs_gen.py for (v,d) = (3,9) (about 4 minutes), (3,10) (about 15 minutes), (4,6) and (8,3) (more than 15
minutes each; see section 4).


2. How to run
-------------
The reference results are in outputs/. Every Python script takes the name of its output file (or folder) as
its first argument, except sat_check.py, whose screen output is the result; the .log files in outputs/ are
screen output saved with ">". From the top folder of this repository, after creating the folders rerun/,
rerun/koszul, rerun/formula, rerun/gen, rerun/j3_integral, rerun/u1_data, rerun/u2_twoblock and rerun/gram:

    python -I code/gram_check.py rerun/gram/gram_check.txt

    python -I code/koszul_h1.py rerun/koszul/h1_v3_d3_p32003.json 32003 3 3
        and likewise h1_v<v>_d<d>_p32003.json for (v,d) = (3,4), (3,5), (3,6), (4,3), (4,4), (5,3), (7,3)

    python -I code/formula_check.py rerun/formula/formula_p32003.json 32003 "2,3 2,4 2,5 2,6 2,7 2,8 2,9 2,10 3,3 3,4 4,3"

    python -I code/ivhs_ts.py rerun/ivhs_ts_p32003.json 32003 4
    python -I code/ivhs_ts.py rerun/ivhs_ts_p1000003.json 1000003 3

    python -I code/ivhs_gen.py rerun/gen/A_v<v>_d<d>_p32003.json 32003 <v> <d> <seeds> <mode>
        with (v, d, seeds, mode) =
        (3,3,2,diag) (3,4,2,diag) (3,5,2,diag) (3,6,2,diag) (3,7,2,diag) (3,8,2,diag)
        (3,9,1,diagiter) (3,10,1,diagiter) (4,3,2,diag) (4,4,2,diag) (4,5,1,diagiter) (4,6,1,diagiter)
        (5,3,1,diagiter) (5,4,1,diagiter) (6,3,1,diag) (7,3,1,diag) (8,3,1,diag)

    python -I code/j3_integral/sat_check.py > rerun/j3_integral/sat_check.log

    python -I code/u1_data/u1_families.py rerun/u1_data 32003 ln:4 2l:4 tc:4 pc:4 ml:4 lc:4 sk:4 pp:3 sc:3 pq:4 ln:3 tc:5 pc:5 ml:5 sk:5 ln:5 > rerun/u1_data/u1_run_p32003.log
    python -I code/u1_data/u1_families.py rerun/u1_data 31991 ln:4 2l:4 tc:4 pc:4 ml:4 lc:4 sk:4 pp:3 sc:3 pq:4 ln:3 tc:5 pc:5 ml:5 sk:5 ln:5 > rerun/u1_data/u1_run_p31991.log
    python -I code/u1_data/u1_crosscheck.py rerun/u1_data 32003 ln:4 2l:4 tc:4 pc:4 ml:4 lc:4 sk:4 pp:3 sc:3 pq:4 ln:3 > rerun/u1_data/u1_crosscheck_p32003.log
    python -I code/u1_data/u1_crosscheck.py rerun/u1_data 31991 sk:5 > rerun/u1_data/u1_crosscheck_sk5_p31991.log
        (u1_crosscheck.py reads the JSON files written by u1_families.py, so it runs after it)

    for P in 32003 31991, the seven runs (D, FAMF, FAMG) = (4, line4, line4), (4, 2lines4, line4),
    (4, 2lines4, 2lines4), (4, line4, gen), (4, conic4, 2lines4), (3, plane6, plane6), (3, plane6, gen), in this
    order, with the screen output of all fourteen runs appended to one file:
    python -I code/u2_twoblock/u2.py rerun/u2_twoblock/u2_FAMF_FAMG_dD_pP.json P D FAMF FAMG >> rerun/u2_twoblock/u2_run.log

    Singular (the image's entry point runs its arguments in a shell, so the Singular command is one argument;
    <repository> is the absolute path of this folder):
    docker run --rm --network none -v <repository>/code:/w:ro -w /w sagemath/sagemath:10.10 "/home/sage/sage/local/bin/Singular -q ivhs_singular.sing" > rerun/ivhs_singular_p32003.txt
    docker run --rm --network none -v <repository>/code:/w:ro -w /w sagemath/sagemath:10.10 "/home/sage/sage/local/bin/Singular -q planes_singular.sing" > rerun/planes_singular_p31981.txt
    docker run --rm --network none -v <repository>/code:/w:ro -w /w sagemath/sagemath:10.10 "/home/sage/sage/local/bin/Singular -q fermat_singular.sing" > rerun/gen/B_fermat_p32003.txt
    (with a local Singular: Singular -q code/ivhs_singular.sing > rerun/ivhs_singular_p32003.txt, and so on)

Then compare the files in rerun/ with those of the same name in outputs/ (for instance with "diff" or "fc").
They should be identical, except for the run times recorded by u1_families.py (the fields time_data_s and
time_kernel_s of its JSON files and the "t=...s" at the end of each line of u1_run_p*.log). The .log files in
outputs/ were saved on Windows and have CRLF line endings; all other files have LF line endings
(.gitattributes keeps every file byte for byte).


3. What each script checks
--------------------------
Each script has a header with the details. All Python scripts except gram_check.py and sat_check.py are exact
over F_p (python-flint nmod_mat); random members are drawn from fixed seeds, so reruns are deterministic.

ivhs_gen.py: library (polynomials as dictionaries, the Jacobian ring R by row reduction degree by degree,
  kernel_dims and kernel_dims_iter: the first-order kernel chain by chain) and script. As a script: the
  first-order kernel of the family of all X_F (tangents A(x) + A(y) for all monomials A of degree d) at random
  F modulo 32003, for the 17 pairs (v,d) of section 2. Expected (diagonal family kernels): dim K(F) = d,
  one dimension per chain (spanned by the classes of the L_eta), for all (v,d) run except (3,3), (3,4), (4,3),
  where the computed values are 19 = 1 + 9 + 9, 52 = 1 + 1 + 1 + 49 and 38 = 1 + 1 + 36 (chains of length one
  contribute (dim R_r)^2). Mode diag stacks all tangent directions; mode diagiter adds them one at a time and
  stops a chain once its kernel has dimension 1 (the lower bound given by the L_eta), which gives the same
  answer faster. The JSON files record, per seed, dim K, dim H^{k,k}_prim, and per chain the number of unknowns
  and the kernel dimension ("L.L" is eps_v). (The script also has a two-block mode F(x) + G(y); its outputs
  are not part of this repository.)
ivhs_ts.py: an earlier, separate implementation for sextics in three variables (fourfolds in P^5, v = 3,
  d = 6), modulo 32003 (4 seeds) and 1000003 (3 seeds). Diagonal family X_F: K = 6, one per chain. Controls in
  the same run: two-block F(x) + G(y) with F, G independent (K = 0); F with a tritangent line (K = 7, one extra
  class from the Noether-Lefschetz locus); the seven-dimensional family of van Geemen and Schuett (K = 333, far
  above its lower bound 8); no tangent directions (K = all of R_12, 1751).
ivhs_singular.sing: second, independent implementation (Singular, Groebner bases in the six variables, no
  Thom-Sebastiani splitting, Singular's own random generator, p = 32003, 3 seeds): dim K = dim of degree 12 of
  S/(J + (tangents)) by Gorenstein duality. Two-block 0, diagonal 6, tritangent 7: agrees with ivhs_ts.py.
planes_singular.sing: the class polynomials of the six planes L_eta on X_F (v = 3, d = 6, p = 31981): the
  determinant formula for the class of a linear cycle, in R_12. All six lie in K, they have rank 6 modulo J,
  and dim K = 6, so K is spanned by the plane classes. Negative controls: a random degree-12 form is not in K,
  and the "plane" y = x (not on X_F) is not in K; the line "decomposition check failed" printed before each
  seed is this intended control (y = x does not give a decomposition of F(x) + F(y)).
fermat_singular.sing: the first-order kernel at the special smooth points F = Fermat and F = Fermat +
  x_1^(d-1) x_2, modulo 32003, 11 runs: Fermat at (v,d) = (3,5), (3,6), (3,7), (4,4), (5,3), (3,3), (3,4),
  (4,3), and the perturbed point at (3,5), (4,4), (5,3): the same values as for random F (d in general; 19, 52, 38 at
  (3,3), (3,4), (4,3)). Output in outputs/gen/B_fermat_p32003.txt.
formula_check.py: the exact formula for dim K(F): #{r <= min(d-1, sigma) : r + d <= sigma} plus the sum of
  (dim R_r)^2 over sigma - d < r <= min(d-1, sigma), against the computed kernel (ivhs_gen.kernel_dims, all
  tangents), for (v,d) = (2,3..10), (3,3), (3,4), (4,3): the cases in which chains of length one occur. The
  JSON file records the total and the formula; the kernel of each chain is printed to the screen (in the
  rerun it was also compared with the formula chain by chain: 1 if r + d <= sigma, else (dim R_r)^2).
koszul_h1.py: the Koszul step used in the proof of the uniform chain lemma: for the Jacobian ring R of a random
  F, the first Koszul cohomology at level j (v-tuples (P_a) in R_j with x_a P_b = x_b P_a, modulo (x_a q))
  vanishes for 1 <= j <= sigma except at j = sigma - d + 2, where it has dimension v. Cases (v,d) = (3,3),
  (3,4), (3,5), (3,6), (4,3), (4,4), (5,3), (7,3), modulo 32003 (seed argument left at its default 0).
gram_check.py: Gram determinants, exactly over Q (fractions), for v = 2..8, d = 3..10: (a) h^k and the d
  spaces L_eta (h^k.h^k = d, h^k.L_eta = 1, L_eta.L_eta = eps_v, distinct L_eta orthogonal): determinant
  d eps_v^(d-1) (eps_v - 1), nonzero for v >= 2; (b) for v = 2c even, h^k, the L_eta and a class Lambda with
  Lambda.Lambda = eps_v, L_eta.Lambda = eps_c, h^k.Lambda = 1: determinant -u^2 eps_v^d, u = (1-d)^c;
  (c) eps_v is the coefficient of H^(v-1) in (1+H)^v/(1+dH). The determinants are computed by Gaussian
  elimination and compared with the closed forms. All 56 cases pass.
sat_check.py (in code/j3_integral/): lattice saturation at the Fermat point, exactly over Z: the d^v linear
  spaces P_a = {y_i = w_i x_i} with the intersection formula P_a.P_b = eps_m, m = #{i : a_i = b_i} (taken as
  input); the Gram rank is (d-1)^v + 1, the claimed basis relations hold in the Gram matrix, and the
  lattice A spanned by h^k and the d spaces L_eta (the P_a with all a_i equal) is saturated in the lattice
  spanned by all d^v spaces (index 1, via Hermite normal form). Controls: replacing one L_eta by twice it gives
  index 2 at (4,4) and (6,3). Cases (v,d) = (2,4), (3,3), (3,5), (4,4), (5,3), (6,3), (4,5), (5,4). It also
  prints det Gram(A), which equals d eps_v^(d-1) (eps_v - 1) (compare gram_check.py). Saturation inside
  H^{2k}(X, Z) is an argument of the note, not checked here.
u1_families.py (in code/u1_data/): Noether-Lefschetz families: F in I_D(d) for a subscheme D of P^{v-1} of
  dimension v/2 - 1 (v = 4: a line ln, two skew lines 2l, a twisted cubic tc, a plane cubic pc, two meeting
  lines ml, a line and a disjoint conic lc, three skew lines sk, a plane quartic pq; v = 6: two disjoint planes
  pp, a cubic scroll sc), tangents A(x) + A(y) for A in I_D(d). For each case: the Hilbert function of
  N = S/(J_F + I_D) and n_D = dim N_{sigma/2}; the first-order kernel per chain; Gorenstein duality
  (0 :_R I_D)_t = N_{sigma - t}; the Gram matrix of h and the components of D on {F = 0} and its rank rho.
  Found at both primes (32003, 31991): dim K = d + n_D^2 (1 + n_D^2 on the chain of residue -v mod d, 1 on the
  others) in the 15 cases with sigma >= 2d; the case ln:3 (v = 4, d = 3, sigma = 4) is outside this range and
  gives 38 = 1 + 1 + 36, the same value as ivhs_gen.py gives for random F at (v,d) = (4,3). n_D = rho - 1 in
  15 cases; in the case pq:4 (a plane quartic, a hyperplane section of {F = 0}) the ideal of D has a generator
  of degree d, the generator condition
  recorded in the JSON (gen_holds: all generators of degree <= d - 1) fails, and n_D = 3 while rho - 1 = 0.
  Each JSON file records these tests as booleans.
u1_crosscheck.py (in code/u1_data/): rebuilds the same F and recomputes the chain kernels with the dense routine
  ivhs_gen.kernel_dims: agrees with u1_families.py in all cases run.
u2.py (in code/u2_twoblock/): two-block families X_{F,G} = {F(x) + G(y) = 0} with F in I_D(d), G in I_D'(d)
  independent, tangents A(x) and B(y) separately: dim K = n_D n_D' (all on the middle chain), and 0 when one
  block is general. Cases: line/line 1, two skew lines/line 2, two skew lines/two skew lines 4, line/general 0,
  conic/two skew lines 2 (d = 4, v = 4); plane/plane 1, plane/general 0 (d = 3, v = 6); primes 32003 and 31991.
ivhs_nl_gen.py: library only here (the families line4, 2lines4, conic4, plane6 used by u2.py); its own script
  mode is not used in this repository.


4. Rerun made when this repository was assembled
------------------------------------------------
All commands of section 2 were rerun in a copy of this repository: the Python scripts with Python 3.12.9,
python-flint 0.9.0 and numpy 2.1.3 (up to eight commands in parallel, with a limit of 15 minutes per
command), and the three Singular scripts in the docker image sagemath/sagemath:10.10 (Singular 4.4.1), without
network access. Comparison of the fresh files with outputs/:
- byte for byte identical (64 files): the outputs of ivhs_ts.py (2), formula_check.py (1), koszul_h1.py (8),
  gram_check.py (1), sat_check.py (1), u2.py (14 JSON files and u2_run.log), u1_crosscheck.py (2 logs), 15 of
  the 17 ivhs_gen.py files, the three Singular outputs, and 16 of the 32 JSON files of u1_families.py;
- identical except for the run times (18 files): the other 16 JSON files of u1_families.py (fields
  time_data_s and time_kernel_s) and u1_run_p32003.log, u1_run_p31991.log (the "t=...s" values);
- not rerun: ivhs_gen.py for (v,d) = (4,6) (mode diagiter) and (8,3) (mode diag) did not finish within the
  15-minute limit and were stopped; outputs/gen/A_v4_d6_p32003.json and A_v8_d3_p32003.json are the original
  results and were not reproduced. The other 15 ivhs_gen.py runs, among them (3,9), (3,10), (4,5), (5,3),
  (5,4) in mode diagiter and (6,3), (7,3) in mode diag, were reproduced; (3,10) took about 15 minutes.
- formula_check.py, chain by chain: the kernel of each chain printed to the screen agrees with the formula
  (1 if r + d <= sigma, else (dim R_r)^2, r = 0..min(d-1, sigma)) in all 11 cases.
