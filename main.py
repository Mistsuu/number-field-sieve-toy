from handleargs       import Args, get_args, check_args, sanitize_args
from polysearch       import polysearch
from rings            import ZZ, x, ZZx
from find_bases       import find_algebraic_factor_bases, find_quadratic_character_bases, find_rational_factor_bases
from sieve            import find_algebraic_and_rational_smooths
from calc             import algebraic_legendre_symbols, ZO_sqrt
from sage.all         import GF, Matrix, gcd, isqrt

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
    # a rational element a + bm presented in prod(a + bm) &
    # an algebraic element a + bO presented in prod(a + bO)
    # such that 
    # g^2    := prod(a + bm) is square in Z &
    # h(O)^2 := prod(a + bO) is square in Z[O] (with high probability)
    for choose_bit_vec in M_F2.left_kernel().basis():
        g2 = 1
        hO2 = 1
        norm_hO = 1
        for choose_bit, smooth_elem in zip(choose_bit_vec, smooth_candidates_info):
            a, b = smooth_elem
            info = smooth_candidates_info[smooth_elem]
            if choose_bit:
                g2 *= (a + b*m)
                hO2 *= (a + b*x)
                hO2 %= f
                for (r, p), e in zip(abases, info['aexp']):
                    norm_hO *= p**e
        norm_hO = isqrt(norm_hO)
        print(f'[debug] {norm_hO = }')
        
        # Since g is just an integer
        # sqrt it is just an easy task.
        g = isqrt(g2)
        print(f'[i] {g = }')

        # Compute square root in Z[O]
        hO = ZO_sqrt(hO2, norm_hO, f)
        if not hO:
            print(f'[i] ------------------------------------------')
            continue
        print(f'[i] {hO = }')

        # Finally, we can map values
        # from Z[O] to Z using map:
        # O -> m.
        h = hO(m)
        print(f'[i] {h = }')

        # From then, we successfully
        # product g and h such that
        # g^2 = h^2 mod n, in which
        # case we can find factors in
        # 2/3 of the cases!
        if 1 < (p := int(gcd(g-h, N))) < N:
            return p
        if 1 < (p := int(gcd(g+h, N))) < N:
            return p

    raise ValueError("co cai nit, but run it again maybe ur lucky")

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