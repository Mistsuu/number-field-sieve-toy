from sage.all import vector, zero_vector, ZZ

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
    N: int    
) -> int:
    sum_baseexp = zero_vector(ZZ, len(abases))
    for baseexp in abaseexps:
        sum_baseexp += vector(ZZ, baseexp)
    print(f'[debug] algebraic base exponents: {sum_baseexp}')
    exit(-1)t