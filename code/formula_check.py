"""Theorem A* formula check (NOTES E1*): dim K(F) = #{r <= min(d-1,sigma) : r + d <= sigma}
+ sum_{sigma-d < r <= min(d-1,sigma)} (dim R_r)^2, against the computed kernel (ivhs_gen.kernel_dims, diag mode),
for random smooth F mod p. Usage: python -I formula_check.py OUT.json prime "v,d v,d ..."
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ivhs_gen import JacRing, kernel_dims, monomials  # noqa: E402


def predicted(R, d):
    s = R.top
    return sum(1 for r in range(min(d - 1, s) + 1) if r + d <= s) + \
        sum(R.dim(r) ** 2 for r in range(min(d - 1, s) + 1) if r + d > s)


def main():
    out, p, cases = sys.argv[1], int(sys.argv[2]), sys.argv[3].split()
    res = []
    for c in cases:
        v, d = map(int, c.split(","))
        rng = random.Random(777 + 31 * v + d)
        alld = monomials(v, d)
        F = {e: rng.randrange(p) for e in alld}
        R = JacRing(F, p, v, d)
        assert R.smooth()
        kd = kernel_dims(R, R, [({m: 1}, {m: 1}) for m in alld], R.top, p)
        K = sum(x["kernel"] for x in kd.values())
        pr = predicted(R, d)
        print("v=%d d=%d sigma=%d K=%d predicted=%d chains=%s %s" % (v, d, R.top, K, pr,
              [x["kernel"] for x in kd.values()], "OK" if K == pr else "MISMATCH"))
        res.append({"v": v, "d": d, "sigma": R.top, "K": K, "predicted": pr})
    with open(out, "w", newline="\n") as fh:
        fh.write(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
