from handleargs       import Args, get_args, check_args, sanitize_args
from polysearch       import polysearch
from rings            import ZZ, x, ZZx
from find_bases       import find_algebraic_factor_bases, find_quadratic_character_bases, find_rational_factor_bases
from sieve            import find_algebraic_and_rational_smooths
from calc             import algebraic_legendre_symbols
from sage.all         import GF, Matrix, gcd
from recover          import recover_rational_square_then_sqrt_it_then_mod_N, recover_algebraic_square_then_sqrt_it_then_do_a_norm_map_then_mod_N

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
    print(f'[i] {f = }')
    print(f'[i] {m = }')

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
    print(f'[i] #rational_factor_bases = {len(rbases)}')
    print(f'[i] #algebraic_factor_bases = {len(abases)}')
    print(f'[i] #quadratic_character_bases = {len(qbases)}')

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
    print(f'[i] Building GF(2) matrix...')
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
    # a + bO is square in Z[O] (with high probability)
    print(f'[i] Solving {M_F2.dimensions()[0]}x{M_F2.dimensions()[1]} GF(2) matrix to derive g != h such that g^2 = h^2 mod n...')
    for choose_bit_vec in M_F2.left_kernel().basis():
        rchooses = []
        achooses = []
        rbaseexps = []
        abaseexps = []
        for choose_bit, smooth_elem in zip(choose_bit_vec, smooth_candidates_info):
            if choose_bit:
                smooth_info = smooth_candidates_info[smooth_elem]
                rchooses.append(smooth_elem[0] + m*smooth_elem[1])
                achooses.append(smooth_elem)
                rbaseexps.append(smooth_info['rexp'])
                abaseexps.append(smooth_info['aexp'])

        g = recover_rational_square_then_sqrt_it_then_mod_N(rchooses, rbases, rbaseexps, int(N))
        h = recover_algebraic_square_then_sqrt_it_then_do_a_norm_map_then_mod_N(achooses, abases, abaseexps, f, m, int(N))
        print(f' * g = {g}')
        print(f' * h = {h}')
        if h and 1 < (p := int(gcd(g-h, N))) < N:
            return p
        if h and 1 < (p := int(gcd(g+h, N))) < N:
            return p
        print(f' * -------------------')

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