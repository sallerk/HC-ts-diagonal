"""Check of the Koszul step in the uniform chain lemma (NOTES E1, Lemma K). R = R_F = S/J_F, F random smooth of degree d
in v variables. Level j: Z_j = {(P_1..P_v) in R_j^v : x_a P_b - x_b P_a = 0 in R_{j+1}}, B_j = {(x_a q) : q in R_{j-1}}.
Prediction (self-duality of the Koszul complex + Tor^S(R, C) = Koszul complex on the partials):
dim Z_j - dim B_j = v if j = sigma - d + 2 (sigma = v(d-2)), else 0 (for 1 <= j <= sigma).
Exact over F_p (python-flint). Usage: python -I koszul_h1.py OUT.json prime v d [seed]
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint  # noqa: E402
from ivhs_gen import JacRing, monomials, pmul  # noqa: E402


def lin_mult(R, a, j):
    """matrix of multiplication by x_a: R_j -> R_{j+1} (rows = std basis of R_{j+1})."""
    e = tuple(1 if i == a else 0 for i in range(R.nv))
    cols = [R.reduce(j + 1, pmul({e: 1}, {R.mons[j][c]: 1}, R.p)) for c in R.std[j]]
    return [[cols[k][i] for k in range(len(cols))] for i in range(R.dim(j + 1))]


def rank(rows, p, ncols):
    if not rows:
        return 0
    return flint.nmod_mat(rows, p).rank() if ncols else 0


def h1(R, j):
    v, p = R.nv, R.p
    rj, rj1 = R.dim(j), R.dim(j + 1)
    M = [lin_mult(R, a, j) for a in range(v)] if rj1 else None
    ncols = v * rj
    rows = []
    if rj1:
        for a in range(v):
            for b in range(a + 1, v):
                for i in range(rj1):
                    row = [0] * ncols
                    for c in range(rj):
                        row[b * rj + c] = M[a][i][c]
                        row[a * rj + c] = (-M[b][i][c]) % p
                    rows.append(row)
    dimZ = ncols - rank(rows, p, ncols)
    if j >= 1 and R.dim(j - 1):
        Mq = [lin_mult(R, a, j - 1) for a in range(v)]
        brows = [[Mq[a][i][c] for c in range(R.dim(j - 1))] for a in range(v) for i in range(rj)]
        dimB = rank(brows, p, R.dim(j - 1))
    else:
        dimB = 0
    return dimZ, dimB


def main():
    out, p, v, d = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    rng = random.Random(4242 + 100 * v + d + 7 * seed)
    F = {e: rng.randrange(p) for e in monomials(v, d)}
    R = JacRing(F, p, v, d)
    assert R.smooth()
    sigma = R.top
    res = {"prime": p, "v": v, "d": d, "sigma": sigma, "predicted_j": sigma - d + 2, "levels": []}
    ok = True
    for j in range(1, sigma + 1):
        dz, db = h1(R, j)
        pred = v if j == sigma - d + 2 else 0
        ok &= (dz - db == pred)
        res["levels"].append({"j": j, "dimZ": dz, "dimB": db, "H1": dz - db, "predicted": pred})
        print("v=%d d=%d j=%d dimZ=%d dimB=%d H1=%d predicted=%d" % (v, d, j, dz, db, dz - db, pred))
    res["all_match"] = ok
    print("all match:", ok)
    with open(out, "w", newline="\n") as fh:
        fh.write(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
