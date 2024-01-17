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

def ZO_sqrt(
    g2: ZZx, 
    f: ZZx
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
    while True:
        q = int(next_prime(q))
        if not f.change_ring(GF(q)).is_irreducible():
            continue
        Fqd = GF(q**d, modulus=f, names=('y'))

        # Thank God Sage has this!
        gq = Fqd(g2).sqrt()
        gq = ZZx([*gq])
        if gQ == None:
            gQ = gq
            Q  = q
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
        g = ZZx([crt(cgQ, cgq, Q, q) for cgQ, cgq in zip(gq, gQ)])
        print(f'[debug] {g = }')
        print(f'[debug] {Q = }')
        print(f'[debug] {q = }')
        print(f'[debug] ---------------------------------')
        if g == gQ:
            return g
        
        if q > 10000:
            break
        

        gQ = g
        Q *= q


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