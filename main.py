from handleargs       import Args, get_args, check_args, sanitize_args
from polysearch       import polysearch
from rings            import ZZ, x, ZZx
from find_bases       import find_algebraic_factor_bases, find_quadratic_character_bases, find_rational_factor_bases
from sieve            import find_algebraic_and_rational_smooths
from calc             import algebraic_legendre_symbols
from sage.all         import GF, Matrix, isqrt

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

    # Find 3 types of bases:
    #   - Rational Factor Bases:
    #       Determines if a + bm is smooth.
    #   - Algebraic Factor Bases:
    #       Determines if a + bO is smooth.
    #   - Quadratic Character Bases:
    #       Determines if a + bO is a square.
    rbases = find_rational_factor_bases(boundZ)
    abases = find_algebraic_factor_bases(f, boundA)
    qbases = find_quadratic_character_bases(f, boundA, boundQ)

    # Find pairs (a,b)
    # such that a+bm
    # is smooth in ZZ,
    # a+bO is smooth in
    # Z[O].
    target_ncandidates = (1 + len(rbases) + len(abases) + len(qbases)) + 1
    smooths, rexps, aexps = \
        find_algebraic_and_rational_smooths(
            rbases,
            abases,
            sieve_a_bound,
            target_ncandidates,
            f, m
        )

    # Add more informations
    # surrounding smooth
    # elements, such as the 
    # sign of the integers
    # a + bm, or the quadratic
    # characters of a + bO
    # modulo some prime elements
    # in <qbases>
    smooth_candidates_info = {}
    for smooth, rexp, aexp in zip(smooths, rexps, aexps):
        a, b = smooth
        qchr = algebraic_legendre_symbols(smooth, qbases)
        smooth_candidates_info[smooth] = {
            'rsig': a + b*m < 0,
            'rexp': rexp,
            'aexp': aexp,
            'qchr': qchr
        }

    # Now we can build up
    # the matrix composed
    # of exponents and
    # quadratic characters
    # encoded as elements
    # modulo 2.
    M = []
    for smooth, info in smooth_candidates_info.items():
        row = []
        row.append(int(info['rsig']))
        for e in info['rexp']:
            row.append(e % 2)
        for e in info['aexp']:
            row.append(e % 2)
        for c in info['qchr']:
            row.append(int(c == -1))
        M.append(row)
        
        # We don't want more...
        if len(M) == target_ncandidates:
            break
    M_F2 = Matrix(GF(2), M)

    # The left kernel forms a vector
    # space that consists of vectors
    # in which the entries are the
    # "choose bits", in which if the
    # entry is 1, that corresponds to
    # a rational factor in a + bm &
    # an algebraic factor in a + bO
    # such that 
    # a + bm is square in Z &
    # a + bO is square in Z[O]
    print(M_F2.left_kernel())
    # for v in M_F2.left_kernel().basis():
    #     g2  = 1
    #     hO2 = 1
    #     for choosebit, (a, b) in zip(v, smooth_candidates_info):
    #         if choosebit:
    #             hO2 *= a + b*x
    #             hO2 %= f
    #             g2  *= a + b*m
    #     g = isqrt(g2)
    #     exit(-1)


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