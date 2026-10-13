"""U2 check: two independent Noether-Lefschetz blocks X_{F,G} = {F(x) + G(y) = 0}, F in I_D(d), G in I_D'(d) chosen
independently (different random seeds). Tangent directions: A(x) for A in I_D(d) and B(y) for B in I_D'(d)
(separately). Hand prediction (NOTES U2): K = n_D * n_D' (all in the middle bidegree (sigma/2, sigma/2)); with
D' = 'gen' (G general, tangents all of S_d) K = 0.
Families from ivhs_nl_gen.family: line4, 2lines4, conic4, plane6; 'gen' = all degree-d forms.
Exact over F_p (python-flint). Usage: python -I u2.py OUT.json prime d famF famG [seed]
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ivhs_gen import JacRing, kernel_dims, monomials  # noqa: E402
from ivhs_nl_gen import family  # noqa: E402

NV = {"line4": 4, "2lines4": 4, "conic4": 4, "plane6": 6}


def block(name, d, p, rng, nv_hint):
    if name == "gen":
        nv = nv_hint
        span = [{m: 1} for m in monomials(nv, d)]
    else:
        nv, span = family(name, d, p, rng)
    for _ in range(5):
        F = {}
        for P in span:
            c = rng.randrange(p)
            for e, a in P.items():
                F[e] = (F.get(e, 0) + c * a) % p
        R = JacRing(F, p, nv, d)
        if R.smooth():
            return nv, span, R
    raise SystemExit("no smooth member for " + name)


def main():
    out, p, d, fa, fb = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    nv = NV.get(fa, NV.get(fb))
    nvF, spanF, RF = block(fa, d, p, random.Random(7001 + seed), nv)
    nvG, spanG, RG = block(fb, d, p, random.Random(9103 + 31 * seed), nv)
    assert nvF == nvG
    N = (RF.top + RG.top) // 2
    tangents = [(A, {}) for A in spanF] + [({}, B) for B in spanG]
    kd = kernel_dims(RF, RG, tangents, N, p)
    chains = [c["kernel"] for c in kd.values()]
    K = sum(chains)
    print("p=%d d=%d F=%s G=%s seed=%d N=%d K=%d chains=%s" % (p, d, fa, fb, seed, N, K, chains))
    res = {"prime": p, "d": d, "famF": fa, "famG": fb, "seed": seed, "N": N, "K": K, "chains": kd}
    with open(out, "w", newline="\n") as fh:
        fh.write(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
