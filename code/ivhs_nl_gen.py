"""E4 implementation A: IVHS bound for X_{F,F} = {F(x) + F(y) = 0}, F very general in a Noether-Lefschetz family
"F in I_d", I the ideal of a fixed subvariety D of P^{v-1} (so C_F = {F = 0} contains D). Tangent directions of the
family modulo J_F: A(x) + A(y), A in I_d. Prediction (NOTES E4): non-invariant chains 1 each; the invariant chain
r = -v mod d (= H^{v-2}_prim(C_F) (x) H^{v-2}_prim(C_F)) carries the extra classes from D.
Families (name: v, ideal generators):
  line4   : v = 4, D = line {x0 = x1 = 0}                     -> I = (x0, x1)
  2lines4 : v = 4, D = two skew lines {x0=x1=0} u {x2=x3=0}    -> I = (x0 x2, x0 x3, x1 x2, x1 x3)
  conic4  : v = 4, D = conic {x0 = 0, q = 0}, q random quadric  -> I = (x0, q)
  plane6  : v = 6, D = plane {x0 = x1 = x2 = 0}               -> I = (x0, x1, x2)
Exact over F_p (python-flint). Usage: python -I ivhs_nl_gen.py OUT.json prime family d [seed]
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ivhs_gen import JacRing, kernel_dims, monomials, pmul  # noqa: E402


def unit(nv, i):
    return tuple(1 if t == i else 0 for t in range(nv))


def family(name, d, p, rng):
    """returns (nv, spanning list of I_d as polynomial dicts)."""
    if name in ("line4", "2lines4", "conic4"):
        nv = 4
    elif name == "plane6":
        nv = 6
    else:
        raise ValueError(name)
    gens = []
    if name == "line4":
        gens = [({unit(nv, 0): 1}, 1), ({unit(nv, 1): 1}, 1)]
    elif name == "2lines4":
        for i in (0, 1):
            for j in (2, 3):
                e = [0] * nv
                e[i] += 1
                e[j] += 1
                gens.append(({tuple(e): 1}, 2))
    elif name == "conic4":
        q = {e: rng.randrange(p) for e in monomials(nv, 2) if e[0] == 0}
        gens = [({unit(nv, 0): 1}, 1), (q, 2)]
    elif name == "plane6":
        gens = [({unit(nv, i): 1}, 1) for i in range(3)]
    span = []
    for g, dg in gens:
        for m in monomials(nv, d - dg):
            span.append(pmul(g, {m: 1}, p))
    return nv, span


def main():
    out, p, name, d = sys.argv[1], int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    rng = random.Random(9000 + 17 * d + seed + sum(map(ord, name)))
    nv, span = family(name, d, p, rng)
    for attempt in range(5):
        F = {}
        for P in span:
            c = rng.randrange(p)
            for e, a in P.items():
                F[e] = (F.get(e, 0) + c * a) % p
        R = JacRing(F, p, nv, d)
        if R.smooth():
            break
    else:
        raise SystemExit("no smooth member found")
    N = R.top
    kd = kernel_dims(R, R, [(A, A) for A in span], N, p)
    chains = [c["kernel"] for c in kd.values()]
    inv = (-nv) % d
    K = sum(chains)
    print("family=%s v=%d d=%d sigma=%d K=%d chains=%s invariant_chain=r%d (value %d)"
          % (name, nv, d, N, K, chains, inv, chains[inv]))
    res = {"prime": p, "family": name, "v": nv, "d": d, "seed": seed, "K": K, "chains": kd,
           "invariant_chain": inv, "num_tangent_spanners": len(span)}
    with open(out, "w", newline="\n") as fh:
        fh.write(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
