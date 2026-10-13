"""Cross-check of u1_families.py: rebuild the same F (same rng sequence) and recompute the chain kernels with the
original dense routine ivhs_gen.kernel_dims (all tangents stacked at once). Also prints HF(S/I_D), info and Gram
data from the u1 JSON files. Usage: python -I u1_crosscheck.py OUTDIR prime case1 [case2 ...] (case = family:d)
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ivhs_gen import JacRing, kernel_dims  # noqa: E402
from u1_families import family, ideal_basis, mons  # noqa: E402


def main():
    outdir, p = sys.argv[1], int(sys.argv[2])
    for case in sys.argv[3:]:
        name, d = case.split(":")
        d = int(d)
        seed = 0
        rng = random.Random(4242 + 1009 * d + 7 * seed + sum(map(ord, name)) * 31 + p)
        nv, comps, info = family(name, d, p, rng)
        Bd = ideal_basis(comps, nv, d, p)
        basis_d = [{e: c for e, c in zip(mons(nv, d)[0], b) if c} for b in Bd]
        for attempt in range(10):
            F = {}
            for P in basis_d:
                c = rng.randrange(p)
                for e, a in P.items():
                    F[e] = (F.get(e, 0) + c * a) % p
            R = JacRing(F, p, nv, d)
            if R.smooth():
                break
        kd = kernel_dims(R, R, [(A, A) for A in basis_d], R.top, p)
        dense = [c["kernel"] for c in kd.values()]
        with open(os.path.join(outdir, "u1_%s_d%d_p%d.json" % (name, d, p))) as fh:
            js = json.load(fh)
        u1 = [js["chains"][r]["kernel"] for r in sorted(js["chains"], key=int)]
        print("%s d=%d p=%d dense_chains=%s u1_chains=%s agree=%s | HF_S/I_D=%s | info=%s | gram=%s rank=%s"
              % (name, d, p, dense, u1, dense == u1, js["HF_S/I_D"][:7], js["info"], js["gram"], js["gram_rank"]),
              flush=True)


if __name__ == "__main__":
    main()
