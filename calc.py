from sage.all import legendre_symbol

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

