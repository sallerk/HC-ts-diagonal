# J3 sanity check: is A = Z<h^k, L_c> saturated in L_{J0} (lattice of the d^v spaces P_a = {y_i = w_i x_i})?
# Fermat point X = {sum x_i^d + sum y_i^d = 0} in P^{2v-1}, dim 2k = 2v-2.
# P_a, a in (Z/d)^v, w_i = c0*zeta^{a_i}.  P_a . P_b = sigma(#{i: a_i=b_i} - 1), sigma(r) = (1-(1-d)^{r+1})/d (AMV (6)).
# Claimed Z-basis of L_{J0}: {P_a : a in [0,d-2]^v} u {h}, with P_b = sum c_b[a] P_a + c_h h (S-model expansion).
import itertools, sys
import numpy as np

def sigma(r, d):
    return (1 - (1 - d) ** (r + 1)) // d

def hnf_index(rows):
    # rows: list of integer lists, full row rank r.  Returns gcd of r x r minors via column HNF (exact, python ints).
    M = [list(r) for r in rows]
    r, n = len(M), len(M[0])
    piv = 0
    diag = []
    for i in range(r):
        # column ops on columns piv..n-1 to make M[i][piv] = gcd, others zero in row i
        for j in range(piv + 1, n):
            a, b = M[i][piv], M[i][j]
            if b == 0:
                continue
            # extended gcd
            x0, y0, x1, y1, aa, bb = 1, 0, 0, 1, a, b
            while bb != 0:
                q = aa // bb
                aa, bb = bb, aa - q * bb
                x0, x1 = x1, x0 - q * x1
                y0, y1 = y1, y0 - q * y1
            g = aa  # = x0*a + y0*b
            if g == 0:
                continue
            u, w = -b // g, a // g
            for t in range(r):
                cp, cj = M[t][piv], M[t][j]
                M[t][piv], M[t][j] = x0 * cp + y0 * cj, u * cp + w * cj
        if M[i][piv] == 0:
            return None  # rank deficient
        diag.append(abs(M[i][piv]))
        piv += 1
    out = 1
    for x in diag:
        out *= x
    return out

def rank_modp(G, p=1000003):
    A = (G % p).astype(np.int64)
    n, m = A.shape
    rk, row = 0, 0
    for col in range(m):
        nz = np.nonzero(A[row:, col])[0]
        if len(nz) == 0:
            continue
        pr = row + nz[0]
        A[[row, pr]] = A[[pr, row]]
        inv = pow(int(A[row, col]), p - 2, p)
        A[row] = (A[row] * inv) % p
        others = np.nonzero(A[:, col])[0]
        others = others[others != row]
        if len(others):
            A[others] = (A[others] - np.outer(A[others, col], A[row])) % p
        row += 1
        rk += 1
        if row == n:
            break
    return rk

def run(v, d, control=False):
    T = list(itertools.product(range(d), repeat=v))
    idx = {a: i for i, a in enumerate(T)}
    N = len(T)
    Arr = np.array(T)
    agree = (Arr[:, None, :] == Arr[None, :, :]).sum(axis=2)
    sig = np.array([sigma(r - 1, d) for r in range(v + 1)], dtype=np.int64)
    G = sig[agree]
    rk = rank_modp(G)
    basis = list(itertools.product(range(d - 1), repeat=v))
    bpos = {a: i for i, a in enumerate(basis)}
    nb = len(basis) + 1  # + h
    hpre = np.zeros(N, dtype=np.int64)
    for t in range(d):
        hpre[idx[(0,) * (v - 1) + (t,)]] += 1
    # S-model expansion of each P_b
    def expand(b):
        polys = []
        for bi in b:
            polys.append({bi: 1} if bi <= d - 2 else {e: -1 for e in range(d - 1)})
        coeff = {}
        for combo in itertools.product(*[list(p.items()) for p in polys]):
            a = tuple(e for e, _ in combo)
            c = 1
            for _, cc in combo:
                c *= cc
            coeff[a] = coeff.get(a, 0) + c
        tot = sum(coeff.values())
        assert (1 - tot) % d == 0
        ch = (1 - tot) // d
        vec = [0] * nb
        for a, c in coeff.items():
            vec[bpos[a]] += c
        vec[-1] = ch
        return vec
    bad = 0
    for b in T:
        vec = expand(b)
        w = np.zeros(N, dtype=np.int64)
        w[idx[b]] += 1
        for a, i in bpos.items():
            w[idx[a]] -= vec[i]
        w -= vec[-1] * hpre
        if np.any(G @ w != 0):
            bad += 1
    # A generators in basis coordinates
    Arows = [expand((j,) * v) for j in range(d)]
    Arows.append([0] * (nb - 1) + [1])
    s = sigma(v - 1, d)
    ind = hnf_index(Arows)
    # Gram of A from G (sanity vs project data)
    gens = [np.eye(N, dtype=np.int64)[idx[(j,) * v]] for j in range(d)] + [hpre]
    GA = np.array([[int(x @ G @ y) for y in gens] for x in gens])
    detA = int(round(np.linalg.det(GA.astype(float))))
    out = f"(v,d)=({v},{d}) N={N} rankG_modp={rk} expected={(d-1)**v+1} badRelations={bad} s={s} detGramA={detA} index[Sat:A]={ind}"
    if control:
        Ac = [list(r) for r in Arows]
        Ac[0] = [2 * x for x in Ac[0]]
        out += f" | control(2*L_c0) index={hnf_index(Ac)}"
        # control 2: replace L_{c_{d-1}} by L_{c0}+L_{c1}+... : check Z<h, L_0..L_{d-2}, sum> same lattice
    print(out, flush=True)

if __name__ == "__main__":
    cases = [(2, 4), (3, 3), (3, 5), (4, 4), (5, 3), (6, 3), (4, 5), (5, 4)]
    for (v, d) in cases:
        run(v, d, control=(v, d) in [(4, 4), (6, 3)])
