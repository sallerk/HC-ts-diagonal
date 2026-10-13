"""Gram determinants of the linear spaces on X_F = {F(x) + F(y) = 0} in P^{2v-1} (dim 2k = 2v - 2, F of degree d).
Exact rational arithmetic (Python fractions; standard library only). With eps_j = (1 - (1 - d)^j)/d:
 (a) the Gram matrix of h^k and the d linear spaces L_eta = {y = eta x} (eta^d = -1), with h^k.h^k = d,
     h^k.L_eta = 1, L_eta.L_eta = eps_v and L_eta.L_eta' = 0 (eta != eta'), has determinant
     d * eps_v^(d-1) * (eps_v - 1), and this is nonzero for v >= 2;
 (b) for v = 2c even, adding a class Lambda with Lambda.Lambda = eps_v, L_eta.Lambda = eps_c and h^k.Lambda = 1,
     the Gram matrix of h^k, the L_eta and Lambda has determinant -u^2 * eps_v^d with u = (1 - d)^c;
 (c) eps_v is the coefficient of H^(v-1) in (1 + H)^v / (1 + d H).
Ranges: v = 2..8, d = 3..10. The determinants are computed by Gaussian elimination over Q, independently of the
closed forms. Usage: python -I gram_check.py OUT.txt
"""
import sys
from fractions import Fraction
from math import comb


def eps(j, d):
    num = 1 - (1 - d) ** j
    assert num % d == 0
    return Fraction(num // d)


def det(M):
    A = [[Fraction(x) for x in row] for row in M]
    n, sign, prod = len(A), 1, Fraction(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if A[r][c] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            sign = -sign
        prod *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            if f:
                A[r] = [A[r][i] - f * A[c][i] for i in range(n)]
    return sign * prod


def gram_a(v, d):
    e = eps(v, d)
    M = [[d] + [1] * d]
    for i in range(d):
        M.append([1] + [e if j == i else 0 for j in range(d)])
    return M


def gram_b(v, d):
    c = v // 2
    e, ec = eps(v, d), eps(c, d)
    M = [[d] + [1] * d + [1]]
    for i in range(d):
        M.append([1] + [e if j == i else 0 for j in range(d)] + [ec])
    M.append([1] + [ec] * d + [e])
    return M


def coeff_c(v, d):
    """coefficient of H^(v-1) in (1 + H)^v * sum_i (-d H)^i."""
    return sum(comb(v, v - 1 - i) * (-d) ** i for i in range(v))


def main():
    out = sys.argv[1]
    lines, ok = [], True
    for v in range(2, 9):
        for d in range(3, 11):
            e = eps(v, d)
            Da = det(gram_a(v, d))
            pa = d * e ** (d - 1) * (e - 1)
            oka = (Da == pa) and Da != 0
            okc = (coeff_c(v, d) == e)
            s = "v=%d d=%d eps_v=%s (a) det=%s formula=%s nonzero=%s %s | (c) coeff=%s %s" % (
                v, d, e, Da, pa, Da != 0, "OK" if oka else "FAIL", coeff_c(v, d), "OK" if okc else "FAIL")
            ok &= oka and okc
            if v % 2 == 0:
                u = (1 - d) ** (v // 2)
                Db = det(gram_b(v, d))
                pb = -u * u * e ** d
                okb = (Db == pb)
                s += " | (b) eps_c=%s det=%s formula=%s %s" % (eps(v // 2, d), Db, pb, "OK" if okb else "FAIL")
                ok &= okb
            lines.append(s)
    lines.append("all checks pass: %s" % ok)
    for s in lines:
        print(s)
    with open(out, "w", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
