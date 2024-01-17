from sage.all import vector, zero_vector, GF, prod, sqrt
from rings    import ZZx, ZZ
from calc     import norm_map_Fpd

def recover_rational_square_then_sqrt_it_then_mod_N(
    rchooses: list[int],  # Note: This variable is not used, but rather put here to fits with the template of recover_xx() functions.
    rbases: list[int],
    rbaseexps: list[int],
    N: int
) -> int:
    sum_baseexp = zero_vector(ZZ, len(rbases))
    for baseexp in rbaseexps:
        sum_baseexp += vector(ZZ, baseexp)
    print(f'[debug] rational base exponents: {sum_baseexp}')

    g = 1
    for p, e in zip(rbases, sum_baseexp):
        g *= pow(p, int(e)//2, N)
        g %= N
    return g

def recover_algebraic_square_then_sqrt_it_then_do_a_norm_map_then_mod_N(
    achooses: list[tuple[int, int]], 
    abases: list[tuple[int, int]],
    abaseexps: list[int],
    f, N: int    
) -> int:
    sum_baseexp = zero_vector(ZZ, len(abases))
    for baseexp in abaseexps:
        sum_baseexp += vector(ZZ, baseexp)
    
    d = f.degree()
    for q in [2, 3, 5, 7]:
        if f.change_ring(GF(q)).is_irreducible():
            # Instead of computing square root
            # in Z[O], we do it in Fp[O] instead.
            Fqd = GF(q**d, modulus=f, names=('y'))
            hO2 = prod(Fqd([a, b]) for (a, b) in achooses)
            if not hO2.is_square():
                continue
            
            # Thank God Sage has this!
            hO = sqrt(hO2)

            # We can categorize which square
            # roots in Fp[O] corresponding
            # to the same square root in Z[O]
            # by checking the norm.
            # 
            # For a in Z[O] where O is the root
            # of an odd-degree f, N(a) = -N(-a).
            #
            norm_hO_Fp_1 = norm_map_Fpd(hO)
            norm_hO_Fp_2 = 1
            for (r, p), e in zip(abases, sum_baseexp):
                norm_hO_Fp_2 *= pow(p, int(e)//2, q)
                norm_hO_Fp_2 %= q
            # if norm_hO_Fp_1 == norm_hO_Fp_2:
            print(f'[debug] {q = }')
            print(f'[debug] {hO = }')
            print(f'[debug] {norm_hO_Fp_1 = }')
            print(f'[debug] {norm_hO_Fp_2 = }')


    exit(-1)