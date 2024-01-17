from sage.all import legendre_symbol, next_prime, GF, crt
from rings    import QQ, ZZx

def algebraic_legendre_symbols(
    elem: tuple[int, int],
    qbases: list[tuple[int, int]]
) -> list[int, int]:
    a, b = elem
    qchr = []
    for (s, q) in qbases:
        if (a + b*s) % q == 0:
            raise ValueError(f"Hmm... Implementation Error! Quadratic Character of this element {elem} modulo prime element {(s, q)} should not be 0.")  
        qchr.append(legendre_symbol(a + b*s, q))
    return qchr

def zfill(
    x: list[int],
    l: int
):
    return x + [0] * (l-len(x))

def ZO_sqrt(
    g2: ZZx, 
    norm_g: int,
    f: ZZx,
    ub_q_safety_limit: int = 50000
) -> ZZx:
    # Iterate over primes q,
    # compute sqrt of h(O)^2
    # over Fq[O] instead of
    # Z[O] since it's hard
    # there.
    d = f.degree()

    q  = 1
    Q  = None
    gQ = None
    gQ_comp = None
    while True:
        q = int(next_prime(q))
        if not f.change_ring(GF(q)).is_irreducible():
            continue
        Fqd = GF(q**d, modulus=f, names=('y'))

        # Thank God Sage has this!
        g2q = Fqd(g2)
        if not g2q.is_square():
            continue
        gq = Fqd(g2).sqrt()
        
        # Use supplied norm to correct
        # the values g(O) mod q,
        # since there are two
        # square roots that exists
        # in F(p^d).
        norm_gq = int(norm_map_Fpd(gq, q, d))
        if not (
            norm_g % q == norm_gq or
            norm_g % q == q - norm_gq
        ):
            print("[warn] Quadratic Character Bases was not big enough to decide whether the values in Z[O] are squares or not! Skipping this psuedo-square element in tears... Or maybe the norm was calculated incorrectly.")
            return None
    
        if norm_g % q != norm_gq:
            gq = -gq
            
        # Convert to integer
        # list to perform CRT.
        gq = ZZx([*gq]).list()
        gq = zfill(gq, d)
        if gQ == None:
            gQ = gq
            gQ_comp = [cgQ-q for cgQ in gQ]
            gQ_comp = zfill(gQ_comp, d)
            Q = q
            continue

        # If after we CRT 2 polymonials in 
        # different rings, the result polynomial
        # is the same as the one in the large
        # ring Z/QZ, it means that all coefficients
        # is < Q in Z and hence they stay the
        # same after mod Q.
        # 
        # In this case, we have successfully found
        # the square root of the element in Z[O]!
        g = ZZx([crt([cgQ, cgq], [Q, q]) 
                    for cgQ, cgq in zip(gQ, gq)]).list()
        g = zfill(g, d)

        # However, the values we obtain
        # after CRT is just in the range
        # [0, Q*q), while coeffients in
        # Z[O] can be negative. So we should
        # keep track the complement values
        # as well.
        Q *= q
        g_comp = [cg-Q for cg in g]
        g_comp = zfill(g_comp, d)

        gZ = []
        for i in range(d):
            if g[i] == gQ[i]:
                gZ.append(g[i])
            elif g_comp[i] == gQ_comp[i]:
                gZ.append(g_comp[i])
            else:
                break
        else:
            return ZZx(gZ)
        
        if q > ub_q_safety_limit:
            raise ValueError("Urg, it seems like sqrt_ZO is going to burn your RAM down. If you don't think this is the case, raise the parameter ub_q_safety_limit.")
        
        gQ = g
        gQ_comp = g_comp


def norm_map_ZO(
    a: int, b: int,
    f
):
    raise NotImplementedError("Not implemented out of fear of not fast enough.")

def norm_map_Fpd(
    a, 
    p: int, d: int
):
    return a ** ((p**d - 1) // (p-1))