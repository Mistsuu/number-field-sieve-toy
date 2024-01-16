from handleargs       import Args, get_args, check_args, sanitize_args
from polysearch       import polysearch
from rings            import ZZ
from find_bases       import find_algebraic_factor_bases, find_quadratic_character_bases, find_rational_factor_bases
from sortedcontainers import SortedDict
from sieve            import find_algebraic_and_rational_smooths
from calc             import algebraic_legendre_symbols

def factor(
    N: ZZ,
    d: ZZ,
    boundZ: int,
    boundA: int,
    boundQ: int,
    sieve_a_bound: tuple[int, int]
) -> int:
    # Find m, f such that
    # f(m) == 0 mod N
    # and deg(f) == d.
    f, m = polysearch(N, d)

    # Find 3 types of
    # bases
    rbases = find_rational_factor_bases(boundZ)
    abases = find_algebraic_factor_bases(f, boundA)
    qbases = find_quadratic_character_bases(f, boundA, boundQ)

    # Find pairs (a,b)
    # such that a+bm
    # is smooth in ZZ,
    # a+bO is smooth in
    # Z[O].
    sieve_b = 1
    smooth_candidates_info = SortedDict({})
    target_ncandidates = 1 + len(rbases) + len(abases) + len(qbases)
    
    while len(smooth_candidates_info) <= target_ncandidates:
        smooths, rexps, aexps = find_algebraic_and_rational_smooths(
                                    rbases,
                                    abases,
                                    sieve_a_bound,
                                    sieve_b,
                                    f, m
                                )

        for smooth, rexp, aexp in zip(smooths, rexps, aexps):
            a, b = smooth
            qchr = algebraic_legendre_symbols(smooth, qbases)
            smooth_candidates_info[smooth] = {
                'sr': a + b*m < 0,
                'er': rexp,
                'ea': aexp,
                'qc': qchr
            }

        # Generate sequence of
        # b = {1, -1, 2, -2, 3, -3, ...}
        sieve_b = -sieve_b
        if sieve_b > 0:
            sieve_b += 1

def main(
    args: Args
) -> None:
    N = ZZ(args.N)
    d = ZZ(args.d)
    boundZ = args.boundZ
    boundA = args.boundA
    boundQ = args.boundQ
    sieve_a_bound = (args.lbsieve_a, args.ubsieve_a)

    print(f'[i] Factoring {N=}...')
    p = factor(N, d, boundZ, boundA, boundQ, sieve_a_bound)
    print(f'[i] Found {p=}!')

if __name__ == '__main__':
    args = get_args()
    sanitize_args(args)
    check_args(args)
    main(args)