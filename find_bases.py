from sage.all import prime_range, GF
from rings    import ZZx, ZZ

def find_rational_factor_bases(
    BZ: int    
) -> list[int]:
    return [-1] + list(map(int, prime_range(2, BZ+1)))

def find_algebraic_factor_bases(
    f: ZZx,
    BA: int,
) -> list[tuple[int, int]]:
    algebraic_factor_bases = []
    for p in prime_range(2, BA+1):
        for r, e in f.change_ring(GF(p)).roots():
            algebraic_factor_bases.append((int(r), int(p)))
    return algebraic_factor_bases

def find_quadratic_character_bases(
    f: ZZx,
    BA: int,
    BQ: int,
) -> list[tuple[int, int]]:
    quadratic_character_bases = []
    for q in prime_range(BA+1, BQ):
        for s, e in f.change_ring(GF(q)).roots():
            quadratic_character_bases.append((int(s), int(q)))
    return quadratic_character_bases
