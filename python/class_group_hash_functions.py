import gmpy2
import hashlib
import sympy

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

        a += 1 - (a % 2)

        if gmpy2.legendre(discriminant, a) == 1 and a.is_prime(miller_rabin_iterations):
            break

    b = tonelli_shanks(discriminant, a)

    if b.is_even():
        b = a - b

    if a % 3 == 1:
        b = -b

    return (a,b)

def clr_hash(message, discriminant, byte_length = 0, n_small_primes = 2, miller_rabin_iterations = 30):

    if not byte_length:
        byte_length = ((gmpy2.isqrt(abs(discriminant) // 4)).bit_length() - 1) // 8

    small_prime_length = byte_length // (2*n_small_primes)
    large_prime_length = byte_length - n_small_primes * small_prime_length

    forms = []
    forms.append(wesolowski_hash(message, discriminant, large_prime_length))
    a = forms[0][0]

    counter = gmpy2.xmpz(0)

    for _ in range(n_small_primes):
        while True:

            a_new, b_new = wesolowski_hash(as_bytes(counter) + message, discriminant, small_prime_length)

            if all((a != a_new for (a, _) in forms)):
                break

            counter += 1

        forms.append((a_new, b_new))
        a *= a_new

    b = gmpy2.mpz(sympy.ntheory.modular.crt(*zip(*forms))[0])
    
    if b.is_even():
        b = a - b

    if a % 3 == 1:
        b = -b

    return (a,b)

def pseudoprime_hash_with_counter(message, discriminant, byte_length = 0):

    if not byte_length:
        byte_length = ((gmpy2.isqrt(abs(discriminant) // 4)).bit_length() - 1) // 8

    counter = gmpy2.xmpz(0)

    while True:

        h = hash_to_bytes(message + as_bytes(counter), byte_length)
        counter += 1
        a = gmpy2.mpz.from_bytes(h)

        a += 3 - (a % 4)

        b = gmpy2.powmod(discriminant % a, (a+1)//4, a)

        if (b**2 - discriminant) % a == 0:
            break

    if b.is_even():
        b = a - b

    if a % 3 == 1:
        b = -b

    return (a,b)

def pseudoprime_hash_with_gap(message, discriminant, byte_length = 0, gap = 0):

    if not gap:
        gap  = discriminant.bit_length()
        gap -= gap % 4

    if not byte_length:
        byte_length = ((gmpy2.isqrt(abs(discriminant) // 4)).bit_length() - 1) // 8

    h  = hash_to_bytes(message + as_bytes(counter), byte_length)
    h *= gap

    a = gmpy2.mpz.from_bytes(h) + 3
    b = gmpy2.powmod(discriminant % a, (a+1)//4, a)

    while True:

        if (b**2 - discriminant) % a == 0:
            break

        a += 4
        b *= discriminant % a
        b %= a

    if b.is_even():
        b = a - b

    if a % 3 == 1:
        b = -b

    return (a,b)

