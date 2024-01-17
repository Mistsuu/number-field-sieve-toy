from sage.all import vector, zero_vector, GF, prod, sqrt, prime_range, next_prime
from rings    import ZZx, ZZ, RR
from calc     import norm_map_Fpd

def recover_rational_square_then_sqrt_it_then_mod_N(
    rchooses: list[int],  # Note: This variable is not used, but rather put here to fits with the template of recover_xx() functions.
    rbases: list[int],
    rbaseexps: list[int],
    N: int
) -> int:
    # Find the exponents of
    # the primes appeared in
    # the factors of the
    # rational square g^2.
    sum_baseexp = zero_vector(ZZ, len(rbases))
    for baseexp in rbaseexps:
        sum_baseexp += vector(ZZ, baseexp)

    # Compute g directly
    # from the exponents of
    # the primes.
    g = 1
    for p, e in zip(rbases, sum_baseexp):
        g *= pow(p, int(e)//2, N)
        g %= N
    return g


def recover_algebraic_square_then_sqrt_it_then_do_a_norm_map_then_mod_N(
    achooses: list[tuple[int, int]], 
    abases: list[tuple[int, int]],
    abaseexps: list[int],
    f, m: int, N: int
) -> int:
    # Find the exponents of
    # the primes appeared in
    # the factors of the
    # algebraic square h(O)^2.
    sum_baseexp = zero_vector(ZZ, len(abases))
    for baseexp in abaseexps:
        sum_baseexp += vector(ZZ, baseexp)

    d = f.degree()
    QN     = 1  # (Q := q0q1q2...) mod N
    list_q = [] # List of qi
    list_T = [] # List of Ti := (Q/qi)^-1 mod qi
    list_h = [] # List of hi := phi(h(O)) mod qi
    
    # Iterate over primes q,
    # compute sqrt of h(O)^2
    # over Fq[O] instead of
    # Z[O] since it's hard
    # there.
    q = 1
    hN = None
    while True:
        q = int(next_prime(q))
        if not f.change_ring(GF(q)).is_irreducible():
            continue

        # Only compute sqrt if
        # h(O)^2 is square in Fq[O],
        # which is odd since h(O)^2 is
        # a square in Z[O], so the
        # property should propagate
        # into mod q as well, no?
        Fqd = GF(q**d, modulus=f, names=('y'))
        hO2q = prod(Fqd([a, b]) for (a, b) in achooses)
        if not hO2q.is_square():
            continue
        
        # Thank God Sage has this!
        hOq = sqrt(hO2q)

        # We can categorize which square
        # roots in Fp[O] corresponding
        # to the same square root in Z[O]
        # by checking the norm.
        # 
        # For a in Z[O] where O is the root
        # of an odd-degree f, N(a) = -N(-a).
        norm_hOq_1 = int(norm_map_Fpd(hOq, q, d))
        norm_hOq_2 = 1
        for (r, p), e in zip(abases, sum_baseexp):
            norm_hOq_2 *= pow(p, int(e)//2, q)
            norm_hOq_2 %= q

        if not (
            norm_hOq_1 == norm_hOq_2 or 
            norm_hOq_1 == q - norm_hOq_2
        ): 
            print("[warn] Quadratic Character Bases was not big enough to decide whether the values in Z[O] are squares or not! Skipping this psuedo-square element in tears...")
            return None

        # Correcting h(O) mod q.
        if norm_hOq_1 % q != norm_hOq_2:
            hOq = -hOq

        # There's always exists a
        # map phi: Z[O] -> Z
        # that sends O  -> m.
        # Here, we do mod q
        # instead, then do CRT...
        hq = GF(q)['x']([*hOq])(m)
        hq = int(hq)

        # Update CRT components
        # in O(N) time.
        T = 1
        for i in range(len(list_q)):
            qi = list_q[i]
            Ti = list_T[i]
            
            Ti *= pow(q, -1, qi)
            Ti %= qi
            list_T[i] = Ti

            T *= pow(qi, -1, q)
            T %= q

        QN *= q
        QN %= N
        list_q.append(q)
        list_T.append(T)
        list_h.append(hq)

        # If after we CRT 2 values in 
        # different rings, the result value
        # is the same as the one in the large
        # ring Z/QZ, it means that the value
        # is < Q in Z and hence they stay the
        # same after mod Q.
        #
        # If we only use the sum formula
        # to compute CRT, the result
        # might not < Q, so we need to 
        # compute r := quotient of that
        # value with Q.
        hN_ = 0
        r = RR(0)
        for i in range(len(list_q)):
            qi = list_q[i]
            Ti = list_T[i]
            hi = list_h[i]

            hN_ += Ti*hi*QN*pow(qi,-1,N)
            hN_ %= N

            r += Ti*hi/RR(qi)

        r = int(r.floor())
        hN_ -= r*QN
        hN_ %= N

        if hN != None and hN == hN_:
            return hN_
        hN = hN_