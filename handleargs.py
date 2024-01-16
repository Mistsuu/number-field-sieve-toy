"""
    handleargs.py:
        Handle parsing arguments.
"""

import argparse
import dataclasses

# ============================================================================
#                       CLASS TO HOLD ARGUMENTS HERE.
# ============================================================================
@dataclasses.dataclass
class Args:
    N: int
    d: int
    boundZ: int
    boundA: int
    boundQ: int
    lbsieve_a: int
    ubsieve_a: int

# ============================================================================
#                       FUNCTION TO PARSE ARGUMENTS.
# ============================================================================
def get_args(

) -> Args:
    # Create default help page.
    parser = argparse.ArgumentParser(
                prog="number_field_sieve_toy",
                description="A toy implementation for the Number Field Sieve Algorithm.",
                epilog="Made by Mistsu ^.^ @01/2024"
             )

    # Arguments config here.
    parser.add_argument('N', type=int, help="Number to factor.")
    parser.add_argument('d', type=int, help="Degree of f(x) = 0 mod n with a ZZ-root m.")
    parser.add_argument('-bZ', '--upper-bound-Z', type=int, help="Upper-bound p for Rational Factor Bases {p}.", default=10000, dest="boundZ")
    parser.add_argument('-bA', '--upper-bound-A', type=int, help="Upper-bound p for Algebraic Factor Bases {(r, p)}.", default=10000, dest="boundA")
    parser.add_argument('-bQ', '--upper-bound-Q', type=int, help="Upper-bound q for Quadratic Character Bases {(s, q)}.", default=20000, dest="boundQ")
    parser.add_argument('-lsa', '--lbound-sieve-a', type=int, help="Lower-bound for a of (a, b) pairs.", default=-100000, dest="lbsieve_a")
    parser.add_argument('-usa', '--ubound-sieve-a', type=int, help="Upper-bound for a of (a, b) pairs.", default=100000, dest="ubsieve_a")

    # Get arguments!
    return Args(**vars(parser.parse_args()))


# ============================================================================
#                       FUNCTION TO CHECK ARGUMENTS.
# ============================================================================
def sanitize_args(
    args: Args
) -> None:
    args.N = abs(args.N)
    args.d = abs(args.d)
    args.boundZ = abs(args.boundZ)
    args.boundA = abs(args.boundA)
    args.boundQ = abs(args.boundQ)

def check_args(
    args: Args
) -> None:
    if args.N == 0 or args.N == 1:
        raise ValueError(f"N mustn't be 0 or 1! N = {args.N}")
    if args.N.bit_length() < args.d:
        raise ValueError(f"Can't find f(x) with d > {args.N.bit_length()}!")
    if args.boundA > args.boundQ:
        raise ValueError(f"Bound q of quadratic characters bases must be larger than bound p of algebraic factor bases (actually it doesn't matter, but this helps implementation looking much cleaner!)")
    