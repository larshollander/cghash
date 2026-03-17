import gmpy2

as_bytes = lambda x: gmpy2.to_binary(x)[2:]
default_rng = gmpy2.random_state()

def sample_prime_discriminant(bit_length, rng = default_rng):

    p = gmpy2.mpz(0)

    while not p.is_prime(60):
        p = 4 * gmpy2.mpz_urandomb(rng, bit_length - 2) + 3

    return -p

def sample_message(bit_length, rng = default_rng):
    return as_bytes(gmpy2.mpz_urandomb(rng, bit_length))
