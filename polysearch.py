"""
    polysearch.py:
        Find f(x) in ZZ[x] and m in ZZ 
        such that:
            f(m) = 0 mod n with degree d.
"""

from random   import randint
from rings    import ZZx, ZZ

def polysearch(
    N: ZZ,
    d: ZZ,
    upper_k = 1000,    # Could be an option but that's not necessary...
):
    k  = randint(1, upper_k)
    kN = k*N
    m  = kN.nth_root(d, truncate_mode=True)[0]
    return ZZx(kN.digits(m)), m