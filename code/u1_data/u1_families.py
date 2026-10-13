"""U1 data: for a subscheme D of P^{v-1} (dim v/2 - 1, ideal I_D) and F random in I_D(d) with C_F smooth:
  - Hilbert function of N = S/(J_F + I_D) in degrees 0..sigma (sigma = v(d-2)), s_D, n_D = dim N_{sigma/2};
  - first-order kernel K(F) = {lambda in (R(x)R)_sigma : (A(x)+A(y)) lambda = 0 for all A in I_D(d)} per chain
    (exact: kernel intersected tangent by tangent over a basis of I_D(d));
  - minimal generator degrees of I_D, HF of S/I_D, Gorenstein duality dim (0 :_R I_D)_t vs dim N_{sigma-t};
  - Gram matrix (hand formulas, see GRAM) and its rank.
I_D(t) is computed degree-wise as the intersection of the component ideals (each component ideal saturated).
Imports ../ivhs_gen.py read-only. Exact over F_p (python-flint).
Usage: python -I u1_families.py OUTDIR prime case1 [case2 ...]   with case = family:d  (e.g. tc:4)
Families: tc twisted cubic; pc plane cubic; ml two meeting lines; lc line + disjoint conic; sk three skew lines;
          pp (v=6) two disjoint planes; sc (v=6) cubic scroll in x5=0; pq plane quartic; ln line; 2l two skew lines.
Output: OUTDIR/u1_<family>_d<d>_p<prime>.json
"""
import json
import os
import random
import sys
import time

import flint

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ivhs_gen import JacRing, monomials, pmul  # noqa: E402

MONS = {}


def mons(nv, n):
    if (nv, n) not in MONS:
        ms = monomials(nv, n)
        MONS[(nv, n)] = (ms, {e: i for i, e in enumerate(ms)})
    return MONS[(nv, n)]


def unit(nv, i):
    return tuple(1 if t == i else 0 for t in range(nv))


def var(nv, i):
    return {unit(nv, i): 1}


def poly(nv, terms, p):
    """terms: list of (coef, exponent tuple)"""
    P = {}
    for c, e in terms:
        P[e] = (P.get(e, 0) + c) % p
    return {e: c for e, c in P.items() if c}


def mono(nv, *idx):
    e = [0] * nv
    for i in idx:
        e[i] += 1
    return tuple(e)


def embed(P3, nv, slots):
    """embed a polynomial in len(slots) variables into nv variables (variable k -> slot[k])."""
    out = {}
    for e, c in P3.items():
        f = [0] * nv
        for k, a in enumerate(e):
            f[slots[k]] += a
        out[tuple(f)] = c
    return out


def rref_rows(rows, ncols, p):
    rows = [r for r in rows if any(r)]
    if not rows:
        return []
    M, rk = flint.nmod_mat(rows, p).rref()
    ent = [int(x) for x in M.entries()]
    return [ent[i * ncols:(i + 1) * ncols] for i in range(rk)]


def vec(P, nv, t):
    ms, idx = mons(nv, t)
    v = [0] * len(ms)
    for e, c in P.items():
        v[idx[e]] = c
    return v


def span_basis(gens, nv, t, p):
    rows = []
    for g, dg in gens:
        for m in mons(nv, t - dg)[0] if t >= dg else []:
            rows.append(vec(pmul(g, {m: 1}, p), nv, t))
    return rref_rows(rows, len(mons(nv, t)[0]), p)


def intersect(U, W, n, p):
    if not U or not W:
        return []
    u = len(U)
    M = flint.nmod_mat(U + W, p).transpose()
    X, nul = M.nullspace()
    if nul == 0:
        return []
    xs = [int(x) for x in X.entries()]
    w = X.ncols()
    A = flint.nmod_mat(nul, u, [xs[i * w + c] for c in range(nul) for i in range(u)], p)
    V = A * flint.nmod_mat(U, p)
    ent = [int(x) for x in V.entries()]
    return rref_rows([ent[i * n:(i + 1) * n] for i in range(nul)], n, p)


def ideal_basis(comps, nv, t, p):
    B = None
    n = len(mons(nv, t)[0])
    for gens in comps:
        S = span_basis(gens, nv, t, p)
        B = S if B is None else intersect(B, S, n, p)
    return B


def rank_rows(rows, p):
    rows = [r for r in rows if any(r)]
    if not rows:
        return 0
    return flint.nmod_mat(rows, p).rank()


def rand_form(nv, deg, p, rng):
    return {e: rng.randrange(p) for e in mons(nv, deg)[0]}


def family(name, d, p, rng):
    """returns (nv, components (list of generator lists [(poly, deg)]), info dict)."""
    info = {}
    if name in ("pp", "sc"):
        nv = 6
    else:
        nv = 4
    X = [var(nv, i) for i in range(nv)]
    if name == "tc":
        q1 = poly(nv, [(1, mono(nv, 0, 2)), (p - 1, mono(nv, 1, 1))], p)
        q2 = poly(nv, [(1, mono(nv, 0, 3)), (p - 1, mono(nv, 1, 2))], p)
        q3 = poly(nv, [(1, mono(nv, 1, 3)), (p - 1, mono(nv, 2, 2))], p)
        comps = [[(q1, 2), (q2, 2), (q3, 2)]]
    elif name in ("pc", "pq"):
        dc = 3 if name == "pc" else 4
        for _ in range(10):
            c3 = rand_form(3, dc, p, rng)
            if JacRing(c3, p, 3, dc).smooth():
                break
        info["plane_curve_smooth"] = JacRing(c3, p, 3, dc).smooth()
        comps = [[(X[0], 1), (embed(c3, nv, [1, 2, 3]), dc)]]
    elif name == "ml":
        comps = [[(X[0], 1), (X[1], 1)], [(X[0], 1), (X[2], 1)]]
    elif name == "lc":
        for _ in range(10):
            q = rand_form(3, 2, p, rng)
            if JacRing(q, p, 3, 2).smooth() and q.get((0, 0, 2), 0):
                break
        info["conic_smooth"] = JacRing(q, p, 3, 2).smooth()
        info["conic_coef_x3^2 (line meets plane x0=0 at (0:0:0:1); disjoint iff nonzero)"] = q.get((0, 0, 2), 0)
        comps = [[(X[0], 1), (embed(q, nv, [1, 2, 3]), 2)], [(X[1], 1), (X[2], 1)]]
    elif name == "sk":
        # lines of one ruling of x0 x3 - x1 x2: {x0 = t x1, x2 = t x3}, t = 0, 1, infinity
        comps = [[(X[0], 1), (X[2], 1)],
                 [(poly(nv, [(1, mono(nv, 0)), (p - 1, mono(nv, 1))], p), 1),
                  (poly(nv, [(1, mono(nv, 2)), (p - 1, mono(nv, 3))], p), 1)],
                 [(X[1], 1), (X[3], 1)]]
    elif name == "pp":
        comps = [[(X[0], 1), (X[1], 1), (X[2], 1)], [(X[3], 1), (X[4], 1), (X[5], 1)]]
    elif name == "sc":
        # 2x2 minors of [[x0, x1, x3], [x1, x2, x4]] in the hyperplane x5 = 0
        m1 = poly(nv, [(1, mono(nv, 0, 2)), (p - 1, mono(nv, 1, 1))], p)
        m2 = poly(nv, [(1, mono(nv, 0, 4)), (p - 1, mono(nv, 1, 3))], p)
        m3 = poly(nv, [(1, mono(nv, 1, 4)), (p - 1, mono(nv, 2, 3))], p)
        comps = [[(X[5], 1), (m1, 2), (m2, 2), (m3, 2)]]
    elif name == "ln":
        comps = [[(X[0], 1), (X[1], 1)]]
    elif name == "2l":
        comps = [[(X[0], 1), (X[1], 1)], [(X[2], 1), (X[3], 1)]]
    else:
        raise ValueError(name)
    return nv, comps, info


def gram(name, d):
    """Gram matrix of <h^{v/2-1}, components of D> on C_F, by the hand formulas of the brief."""
    if name == "tc":
        return [[d, 3], [3, -2 - 3 * (d - 4)]]
    if name == "pc":
        return [[d, 3], [3, -3 * (d - 4)]]
    if name == "pq":
        return [[d, 4], [4, 4 - 4 * (d - 4)]]
    if name == "ml":
        return [[d, 1, 1], [1, 2 - d, 1], [1, 1, 2 - d]]
    if name == "lc":
        return [[d, 1, 2], [1, 2 - d, 0], [2, 0, -2 - 2 * (d - 4)]]
    if name == "sk":
        return [[d, 1, 1, 1], [1, 2 - d, 0, 0], [1, 0, 2 - d, 0], [1, 0, 0, 2 - d]]
    if name == "pp":
        return [[3, 1, 1], [1, 3, 0], [1, 0, 3]]
    if name == "sc":
        return [[3, 3], [3, 7]]
    if name == "ln":
        return [[d, 1], [1, 2 - d]]
    if name == "2l":
        return [[d, 1, 1], [1, 2 - d, 0], [1, 0, 2 - d]]
    return None


class NF:
    def __init__(self, R, p):
        self.R, self.p, self.nv, self.top = R, p, R.nv, R.top
        self.T = {}
        for t in range(self.top + 1):
            r, n = len(R.std[t]), len(R.mons[t])
            if r == 0:
                self.T[t] = None
                continue
            M = flint.nmod_mat(r, n, p)
            for k, j in enumerate(R.std[t]):
                M[k, j] = 1
            for row, pj in zip(R.rows[t], R.piv[t]):
                for k, js in enumerate(R.std[t]):
                    c = row[js]
                    if c:
                        M[k, pj] = (p - c) % p
            self.T[t] = M

    def dim(self, t):
        return self.R.dim(t)

    def image_rank(self, t, rows):
        """rank of the image in R_t of the span of the given S_t vectors."""
        if self.dim(t) == 0 or not rows:
            return 0
        B = flint.nmod_mat(rows, self.p).transpose()
        return (self.T[t] * B).rank()

    def mult(self, A, e, a):
        """matrix (dim R_{a+e} x dim R_a) of multiplication by A (degree e), as nmod_mat; None if zero."""
        if a < 0 or a + e > self.top or self.dim(a) == 0 or self.dim(a + e) == 0:
            return None
        ms, idx = mons(self.nv, a + e)
        P = flint.nmod_mat(len(ms), self.dim(a), self.p)
        for k, js in enumerate(self.R.std[a]):
            m = self.R.mons[a][js]
            for ex, c in A.items():
                i = idx[tuple(x + y for x, y in zip(ex, m))]
                P[i, k] = (int(P[i, k]) + c) % self.p
        return self.T[a + e] * P


def tolist(M):
    nr, nc = M.nrows(), M.ncols()
    ent = [int(x) for x in M.entries()]
    return [ent[i * nc:(i + 1) * nc] for i in range(nr)]


def chain_kernels(nf, tangents, d, p):
    N = nf.top
    out = {}
    mcache = {}

    def mm(ti, a):
        if (ti, a) not in mcache:
            M = nf.mult(tangents[ti], d, a)
            mcache[(ti, a)] = tolist(M) if M is not None else None
        return mcache[(ti, a)]

    for r0 in range(d):
        chain = [a for a in range(r0, N + 1, d) if nf.dim(a) and nf.dim(N - a)]
        if not chain:
            continue
        off, nvar = {}, 0
        for a in chain:
            off[a] = nvar
            nvar += nf.dim(a) * nf.dim(N - a)
        state = {"Kb": None, "k": nvar}
        pending = []

        def flush():
            if not pending or state["k"] == 0:
                pending.clear()
                return
            E = flint.nmod_mat(len(pending), nvar, p)
            for r, row in enumerate(pending):
                for c, v in row.items():
                    E[r, c] = v
            pending.clear()
            M = E if state["Kb"] is None else E * state["Kb"]
            X, nul = M.nullspace()
            k_old = X.nrows()
            xs = [int(x) for x in X.entries()]
            w = X.ncols()
            Xs = flint.nmod_mat(k_old, nul, [xs[i * w + c] for i in range(k_old) for c in range(nul)], p) if nul else None
            if nul == 0:
                state["Kb"], state["k"] = None, 0
            else:
                state["Kb"] = Xs if state["Kb"] is None else state["Kb"] * Xs
                state["k"] = nul

        used = 0
        for ti in range(len(tangents)):
            used += 1
            for a in range(r0, N + 1, d):
                b = a + d
                ra, rb = nf.dim(a), nf.dim(b)
                sa, sb = nf.dim(N - a), nf.dim(N - b)
                if rb == 0 or sa == 0:
                    continue
                MA = mm(ti, a) if a in off else None
                MB = mm(ti, N - b) if b in off else None
                if MA is None and MB is None:
                    continue
                for i in range(rb):
                    for j in range(sa):
                        row = {}
                        if MA is not None:
                            for k in range(ra):
                                c = MA[i][k]
                                if c:
                                    row[off[a] + k * sa + j] = c
                        if MB is not None:
                            for l in range(sb):
                                c = MB[j][l]
                                if c:
                                    row[off[b] + i * sb + l] = c
                        if row:
                            pending.append(row)
            if len(pending) >= state["k"]:
                flush()
            if state["k"] == 0:
                break
        flush()
        out[r0] = {"chain_degrees": chain, "unknowns": nvar, "kernel": state["k"], "tangents_used": used}
    return out


def run_case(outdir, p, name, d, seed=0):
    t0 = time.time()
    rng = random.Random(4242 + 1009 * d + 7 * seed + sum(map(ord, name)) * 31 + p)
    nv, comps, info = family(name, d, p, rng)
    sigma = nv * (d - 2)
    # ideal degree pieces
    Ib = {t: ideal_basis(comps, nv, t, p) for t in range(0, sigma + 1)}
    dimI = {t: len(Ib[t]) for t in Ib}
    hf_SI = [len(mons(nv, t)[0]) - dimI[t] for t in range(sigma + 1)]
    mingens = {}
    for t in range(1, sigma + 1):
        prev = []
        for b in Ib[t - 1]:
            P = {e: c for e, c in zip(mons(nv, t - 1)[0], b) if c}
            for i in range(nv):
                prev.append(vec(pmul(P, var(nv, i), p), nv, t))
        g = dimI[t] - rank_rows(prev, p)
        if g:
            mingens[t] = g
    # F random in I_D(d)
    basis_d = [{e: c for e, c in zip(mons(nv, d)[0], b) if c} for b in Ib[d]]
    for attempt in range(10):
        F = {}
        for P in basis_d:
            c = rng.randrange(p)
            for e, a in P.items():
                F[e] = (F.get(e, 0) + c * a) % p
        R = JacRing(F, p, nv, d)
        if R.smooth():
            break
    smooth = R.smooth()
    nf = NF(R, p)
    dimR = [nf.dim(t) for t in range(sigma + 1)]
    hfN = [dimR[t] - nf.image_rank(t, Ib[t]) for t in range(sigma + 1)]
    sD = max([t for t in range(sigma + 1) if hfN[t]], default=-1)
    half = sigma // 2
    nD = hfN[half]
    # Gorenstein duality: (0 :_R I_D)_t vs N_{sigma - t}; generators = bases of I_D(t) for t <= max min-gen degree
    gdeg = max(mingens) if mingens else 0
    gens = []
    for t in range(1, gdeg + 1):
        if t in mingens:
            gens += [({e: c for e, c in zip(mons(nv, t)[0], b) if c}, t) for b in Ib[t]]
    ann = []
    for t in range(sigma + 1):
        rows = []
        for g, e in gens:
            M = nf.mult(g, e, t)
            if M is not None:
                rows += tolist(M)
        ann.append(dimR[t] - rank_rows(rows, p))
    gor = all(ann[t] == hfN[sigma - t] for t in range(sigma + 1))
    t1 = time.time()
    # kernel
    kd = chain_kernels(nf, basis_d, d, p)
    K = sum(c["kernel"] for c in kd.values())
    inv = (-nv) % d
    t2 = time.time()
    G = gram(name, d)
    rho = flint.fmpz_mat(G).rank() if G else None
    gen_ok = (max(mingens) <= d - 1) if mingens else True
    W_ok = sD <= half
    sig_ok = sigma >= 2 * d
    res = {"prime": p, "family": name, "v": nv, "d": d, "seed": seed, "sigma": sigma, "C_F_smooth": smooth,
           "info": info, "HF_S/I_D": hf_SI, "min_generators_by_degree": mingens, "dim_I_D(d)": len(basis_d),
           "HF_R": dimR, "HF_N": hfN, "s_D": sD, "n_D": nD, "dim_R_half": dimR[half],
           "ann_0:I_D": ann, "gorenstein_duality_ann_t_eq_N_sigma_minus_t": gor,
           "K": K, "chains": {str(r): c for r, c in kd.items()}, "invariant_residue": inv,
           "invariant_chain_value": kd[inv]["kernel"] if inv in kd else None,
           "gram": G, "gram_rank": rho, "rho_minus_1": (rho - 1) if rho else None,
           "gen_holds": gen_ok, "W_holds": W_ok, "sigma_ge_2d": sig_ok,
           "NL_nD_eq_rho_minus_1": (nD == rho - 1) if rho else None,
           "K_eq_d_plus_nD2": K == d + nD * nD,
           "pred_chains_ok": all((c["kernel"] == (1 + nD * nD if r == inv else 1)) for r, c in kd.items()),
           "time_data_s": round(t1 - t0, 1), "time_kernel_s": round(t2 - t1, 1)}
    fn = os.path.join(outdir, "u1_%s_d%d_p%d.json" % (name, d, p))
    with open(fn, "w", newline="\n") as fh:
        fh.write(json.dumps(res, indent=1) + "\n")
    print("%s d=%d p=%d smooth=%s sigma=%d HF_N=%s s_D=%d n_D=%d R_half=%d K=%d chains=%s inv=r%d rho-1=%s "
          "gen=%s W=%s gor=%s K==d+n^2:%s chains_ok:%s mingens=%s t=%.1f+%.1fs"
          % (name, d, p, smooth, sigma, hfN, sD, nD, dimR[half], K, [kd[r]["kernel"] for r in sorted(kd)], inv,
             res["rho_minus_1"], gen_ok, W_ok, gor, res["K_eq_d_plus_nD2"], res["pred_chains_ok"], mingens,
             t1 - t0, t2 - t1), flush=True)
    return res


def main():
    outdir, p = sys.argv[1], int(sys.argv[2])
    assert flint.fmpz(p).is_prime()
    for case in sys.argv[3:]:
        name, d = case.split(":")
        run_case(outdir, p, name, int(d))


if __name__ == "__main__":
    main()
