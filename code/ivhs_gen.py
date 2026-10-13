"""E1 implementation A: IVHS bound for X_{F,F} = {F(x) + F(y) = 0} in P^{2v-1}, F a form of degree d in v variables
(n = 2v - 2 = 2k, k = v - 1). Thom-Sebastiani: R_X = R_F (x) R_F, R_F of socle degree s = v(d - 2).
H^{k,k}_prim(X) = (R_X)_N with N = (k + 1)(d - 2); a tangent direction A(x) + A(y) (A of degree d) maps it to
(R_X)_{N+d} = H^{k-1,k+1}_prim. K = {lambda : (A(x) + A(y)) lambda = 0 for all A}. With lambda = sum_a L_a,
L_a in (R_F)_a (x) (R_F)_{N-a}, the component in bidegree (a+d, N-a) is M_A^{(a)} L_a + L_{a+d} (M_A^{(N-a-d)})^T;
chains are residue classes of a mod d. Lower bound: the d disjoint linear spaces {y = c x}, c^d = -1 (Gram rank d + 1
for d >= 3). Also: two-block (F, G independent) mode for P1 (expect 0).
Exact over F_p (python-flint). Usage: python -I ivhs_gen.py OUT.json prime v d [seeds] [mode: diag|two|diagrandT (T random tangents)]
A subset of tangent directions gives a larger kernel, so diagrandT is still a valid upper bound.
"""
import itertools
import json
import random
import sys

import flint


def monomials(nv, n):
    if n < 0:
        return []
    return [e for e in itertools.product(range(n + 1), repeat=nv) if sum(e) == n]


def pmul(P, Q, p):
    R = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = tuple(a + b for a, b in zip(e1, e2))
            R[e] = (R.get(e, 0) + c1 * c2) % p
    return {e: v for e, v in R.items() if v}


def pderiv(P, i, p):
    R = {}
    for e, c in P.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            f = tuple(f)
            R[f] = (R.get(f, 0) + c * e[i]) % p
    return {e: v for e, v in R.items() if v}


class JacRing:
    def __init__(self, F, p, nv, d):
        self.p, self.nv, self.d = p, nv, d
        self.top = nv * (d - 2)
        parts = [pderiv(F, i, p) for i in range(nv)]
        self.mons, self.idx, self.rows, self.piv, self.std = {}, {}, {}, {}, {}
        for n in range(self.top + 2):
            mons = monomials(nv, n)
            idx = {e: k for k, e in enumerate(mons)}
            gens = []
            for g in parts:
                for m in monomials(nv, n - (d - 1)):
                    row = [0] * len(mons)
                    for e, c in g.items():
                        row[idx[tuple(a + b for a, b in zip(e, m))]] = c
                    gens.append(row)
            if gens:
                Rm, rk = flint.nmod_mat(gens, p).rref()
                ent = [int(x) for x in Rm.entries()]
                nc = len(mons)
                rows = [ent[i * nc:(i + 1) * nc] for i in range(rk)]
            else:
                rows = []
            piv = [next(j for j, v in enumerate(r) if v) for r in rows]
            ps = set(piv)
            self.mons[n], self.idx[n], self.rows[n], self.piv[n] = mons, idx, rows, piv
            self.std[n] = [j for j in range(len(mons)) if j not in ps]

    def dim(self, n):
        return len(self.std[n]) if 0 <= n <= self.top else 0

    def smooth(self):
        return len(self.std[self.top + 1]) == 0 and len(self.std[self.top]) == 1

    def reduce(self, n, P):
        p = self.p
        v = [0] * len(self.mons[n])
        for e, c in P.items():
            v[self.idx[n][e]] = (v[self.idx[n][e]] + c) % p
        for r, j in zip(self.rows[n], self.piv[n]):
            c = v[j]
            if c:
                v = [(a - c * b) % p for a, b in zip(v, r)]
        return [v[j] for j in self.std[n]]

    def mult_matrix(self, A, a):
        """(r_{a+d} x r_a) matrix of multiplication by A: (R_F)_a -> (R_F)_{a+d}; None if the target is 0."""
        if a + self.d > self.top or self.dim(a) == 0:
            return None
        cols = [self.reduce(a + self.d, pmul(A, {self.mons[a][j]: 1}, self.p)) for j in self.std[a]]
        return [[cols[k][i] for k in range(len(cols))] for i in range(self.dim(a + self.d))]


def kernel_dims(RF, RG, tangents, N, p):
    """Kernel dimension per chain. The equation matrix is filled sparsely into a zero nmod_mat."""
    d = RF.d
    out = {}
    mcache = {}

    def mm(R, A, a, key):
        kk = (key, a)
        if kk not in mcache:
            mcache[kk] = R.mult_matrix(A, a)
        return mcache[kk]

    for r0 in range(d):
        chain = [a for a in range(r0, N + 1, d) if RF.dim(a) and RG.dim(N - a)]
        if not chain:
            continue
        off, nv = {}, 0
        for a in chain:
            off[a] = nv
            nv += RF.dim(a) * RG.dim(N - a)
        blocks = []  # (rows, entries)
        nrows = 0
        entries = []
        for ti, (A, B) in enumerate(tangents):
            for a in range(r0, N + 1, d):
                b = a + d
                ra, rb = RF.dim(a), RF.dim(b)
                sa, sb = RG.dim(N - a), RG.dim(N - b)
                if rb == 0 or sa == 0:
                    continue
                MA = mm(RF, A, a, ("F", ti)) if (A and a in off) else None
                MB = mm(RG, B, N - b, ("G", ti)) if (B and b in off) else None
                if MA is None and MB is None:
                    continue
                for i in range(rb):
                    for j in range(sa):
                        r = nrows
                        nrows += 1
                        if MA is not None:
                            for k in range(ra):
                                c = MA[i][k]
                                if c:
                                    entries.append((r, off[a] + k * sa + j, c))
                        if MB is not None:
                            for l in range(sb):
                                c = MB[j][l]
                                if c:
                                    entries.append((r, off[b] + i * sb + l, c))
        if nrows:
            M = flint.nmod_mat(nrows, nv, p)
            for r, c, v in entries:
                M[r, c] = (int(M[r, c]) + v) % p
            rk = M.rank()
        else:
            rk = 0
        out[str(chain)] = {"unknowns": nv, "kernel": nv - rk}
    return out


def kernel_dims_iter(RF, RG, tangents, N, p, floor=1):
    """Same kernel as kernel_dims, computed incrementally: tangents are added one at a time to a row-echelon basis
    of the equation space; a chain stops once its kernel reaches `floor` (the proven lower bound per chain), since
    more tangents can only shrink the kernel. Returns kernel per chain and the number of tangents used."""
    d = RF.d
    out = {}
    for r0 in range(d):
        chain = [a for a in range(r0, N + 1, d) if RF.dim(a) and RG.dim(N - a)]
        if not chain:
            continue
        off, nv = {}, 0
        for a in chain:
            off[a] = nv
            nv += RF.dim(a) * RG.dim(N - a)
        basis = []  # list of dense rows (echelon basis)
        used = 0
        for (A, B) in tangents:
            used += 1
            new_rows = []
            for a in range(r0, N + 1, d):
                b = a + d
                ra, rb = RF.dim(a), RF.dim(b)
                sa, sb = RG.dim(N - a), RG.dim(N - b)
                if rb == 0 or sa == 0:
                    continue
                MA = RF.mult_matrix(A, a) if (A and a in off) else None
                MB = RG.mult_matrix(B, N - b) if (B and b in off) else None
                if MA is None and MB is None:
                    continue
                for i in range(rb):
                    for j in range(sa):
                        row = [0] * nv
                        if MA is not None:
                            for k in range(ra):
                                c = MA[i][k]
                                if c:
                                    row[off[a] + k * sa + j] = c
                        if MB is not None:
                            for l in range(sb):
                                c = MB[j][l]
                                if c:
                                    row[off[b] + i * sb + l] = (row[off[b] + i * sb + l] + c) % p
                        if any(row):
                            new_rows.append(row)
            if not new_rows:
                continue
            Rm, rk = flint.nmod_mat(basis + new_rows, p).rref()
            ent = [int(x) for x in Rm.entries()]
            basis = [ent[i * nv:(i + 1) * nv] for i in range(rk)]
            if nv - rk <= floor:
                break
        out[str(chain)] = {"unknowns": nv, "kernel": nv - len(basis), "tangents_used": used}
    return out


def gram_s(k, d):
    return (1 - (1 - d) ** (k + 1)) // d


def main():
    out = sys.argv[1]
    p, nv, d = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    seeds = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    mode = sys.argv[6] if len(sys.argv) > 6 else "diag"
    k = nv - 1
    N = (k + 1) * (d - 2)
    res = {"prime": p, "v": nv, "d": d, "k": k, "N": N, "mode": mode, "L.L": gram_s(k, d), "runs": []}
    alld = monomials(nv, d)
    for s in range(seeds):
        rng = random.Random(31337 + 1000 * nv + 10 * d + s)
        F = {e: rng.randrange(p) for e in alld}
        RF = JacRing(F, p, nv, d)
        assert RF.smooth(), "not smooth"
        if mode in ("diag", "diagiter"):
            RG, tans = RF, [({m: 1}, {m: 1}) for m in alld]
            expect = d
        elif mode.startswith("diagrand"):
            nt = int(mode[len("diagrand"):])
            RG = RF
            tans = []
            for _ in range(nt):
                A = {e: rng.randrange(p) for e in alld}
                tans.append((A, A))
            expect = d
        else:
            G = {e: rng.randrange(p) for e in alld}
            RG = JacRing(G, p, nv, d)
            assert RG.smooth()
            tans = [({m: 1}, {}) for m in alld] + [({}, {m: 1}) for m in alld]
            expect = 0
        if mode == "diagiter":
            kd = kernel_dims_iter(RF, RF, [({m: 1}, {m: 1}) for m in alld], N, p, floor=1)
        else:
            kd = kernel_dims(RF, RG, tans, N, p)
        K = sum(c["kernel"] for c in kd.values())
        tot = sum(c["unknowns"] for c in kd.values())
        print("v=%d d=%d k=%d mode=%s seed=%d dim H^{k,k}_prim=%d K=%d expected(lower bound)=%d chains=%s"
              % (nv, d, k, mode, s, tot, K, expect, [c["kernel"] for c in kd.values()]))
        res["runs"].append({"seed": s, "K": K, "h_kk_prim": tot, "chains": kd})
    with open(out, "w", newline="\n") as fh:
        fh.write(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
