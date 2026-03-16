import gmpy2
import hashlib

from modsqrt import tonelli_shanks


hash_to_bytes = lambda x, length: hashlib.shake_256(x).digest(length)
as_bytes = lambda x: gmpy2.to_binary(x)[2:]


def wesolowski_hash(message, discriminant, byte_length = 0, miller_rabin_iterations = 30):
    
    if not byte_length:
        byte_length = ((gmpy2.isqrt(abs(discriminant) // 4)).bit_length() - 1) // 8

    counter = gmpy2.xmpz(0)

    while True:

        h = hash_to_bytes(message + as_bytes(counter), byte_length)
        counter += 1
        a = gmpy2.mpz.from_bytes(h)

        if a % 2 == 0:
            a += 1

        if gmpy2.legendre(discriminant, a) == 1 and a.is_prime(miller_rabin_iterations):
            break

    b = tonelli_shanks(discriminant, a)

    if b.is_even():
        b = a - b

    if a % 3 == 1:
        b = -b

    return (a,b)

