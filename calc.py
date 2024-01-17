from sage.all import legendre_symbol
from rings    import QQ

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