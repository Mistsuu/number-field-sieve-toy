from sortedcontainers import SortedDict
from rings            import QQ

def find_algebraic_and_rational_smooths(
    rbases: list[int],
    abases: list[tuple[int, int]],
    sieve_a_bound: tuple[int, int],
    sieve_b: int,
    f, m: int
) -> tuple[
        list[tuple[int, int]],  # smooths
        list[list[int]],        # exponents of rational factors
        list[list[int]]         # exponents of algebraic factors
    ]:

    lb_a, ub_a = sieve_a_bound
    b          = sieve_b
    nelems     = ub_a - lb_a + 1

    f = f.change_ring(QQ)
    d = f.degree()

    # Building sieve arrays
    rsieve_arr = []
    asieve_arr = []
    for a in range(lb_a, ub_a + 1):
        rsieve_arr.append(a + b*m)
        asieve_arr.append(int((-b)**d * f(QQ(-a)/b)))

    # Sieving rational array.
    rexps = [[0] * len(rbases) for _ in range(nelems)]
    for ibase, p in enumerate(rbases):
        ielem_start = (-b*m - lb_a) % p
        for ielem in range(ielem_start, nelems, p):
            while rsieve_arr[ielem] % p == 0:
                rsieve_arr[ielem] //= p
                rexps[ielem][ibase] += 1        

    # Sieving algebraic array.
    aexps = [[0] * len(abases) for _ in range(nelems)]
    for ibase, (r, p) in enumerate(abases):
        ielem_start = (-b*r - lb_a) % p
        for ielem in range(ielem_start, nelems, p):
            while asieve_arr[ielem] % p == 0:
                asieve_arr[ielem] //= p
                aexps[ielem][ibase] += 1

    # Filter indices i where rsieve_arr[i] 
    # and asieve_arr[i] == 1 at the same
    # time.
    smooths = []
    rexps_filtered = []
    aexps_filtered = []
    for ielem in range(nelems):
        if rsieve_arr[ielem] == asieve_arr[ielem] == 1:
            smooths.append((ielem + lb_a, b))
            rexps_filtered.append(rexps[ielem])
            aexps_filtered.append(aexps[ielem])

    return smooths, rexps_filtered, aexps_filtered