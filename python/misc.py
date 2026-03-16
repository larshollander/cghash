import gmpy2

default_rng = gmpy2.random_state()

def sample_prime_discriminant(bit_length, rng = default_rng):

    p = gmpy2.mpz(0)

    while not p.is_prime(60):
        p = 4 * gmpy2.mpz_urandomb(rng, bit_length - 2) + 3

    return -p
