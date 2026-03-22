import gmpy2

as_bytes = lambda x: gmpy2.to_binary(x)[2:]
default_rng = gmpy2.random_state()

primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53)

def sample_prime_discriminant(bit_length, rng = default_rng, mr_reps = 30):

    p = gmpy2.mpz(0)

    while not p.is_prime(mr_reps):
        p = 4 * gmpy2.mpz_urandomb(rng, bit_length - 2) + 3

    return -p

def sample_message(bit_length, rng = default_rng):

    return as_bytes(gmpy2.mpz_urandomb(rng, bit_length))

def is_possibly_prime(n, n_primes = 16):

    return all((n % p for p in primes[:n_primes]))
