import gmpy2
import hashlib
import sympy

from misc import is_possibly_prime, qnr_primes
from modsqrt import tonelli_shanks


hash_to_bytes = lambda x, length: hashlib.shake_256(x).digest(length)
as_bytes = lambda x: gmpy2.to_binary(x)[2:].ljust(4, b'\x00')
as_hex = lambda xs: "".join("{:02x}".format(x) for x in xs)


def wesolowski_hash(message, discriminant, byte_length = 0, miller_rabin_iterations = 30):
    
    if not byte_length:
        byte_length = (discriminant.bit_length()//2 - 1) // 8

    counter = gmpy2.xmpz(0)

    while True:

        h = hash_to_bytes(message + as_bytes(counter), byte_length)
        counter += 1

        a = gmpy2.mpz.from_bytes(h)
        swap_b = a % 2

        a += 3 - (a % 4)

        if gmpy2.legendre(discriminant, a) == 1 and a.is_prime(miller_rabin_iterations):
            break

    b = tonelli_shanks(discriminant, a)

    if b.is_even():
        b = a - b

    if swap_b:
        b = -b

    return ((a,b), counter)

def clr_hash(message, discriminant, byte_length = 0, n_small_primes = 2, miller_rabin_iterations = 30):

    if not byte_length:
        byte_length = (discriminant.bit_length()//2 - 1) // 8

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

    return ((a,b), counter)

def pseudoprime_hash_with_counter(message, discriminant, byte_length = 0):

    if not byte_length:
        byte_length = (discriminant.bit_length()//2 - 1) // 8

    counter = gmpy2.xmpz(0)

    while True:

        h = hash_to_bytes(message + as_bytes(counter), byte_length)
        counter += 1

        a = gmpy2.mpz.from_bytes(h)
        swap_b = a % 2

        a += 3 - (a % 4)

        if is_possibly_prime(a):

            b = gmpy2.powmod(discriminant % a, (a+1)//4, a)

            if (b**2 - discriminant) % a == 0:
                break

    if b.is_even():
        b = a - b

    if swap_b:
        b = -b

    return ((a,b), counter)

def pseudoprime_hash_with_gap(message, discriminant, byte_length = 0, gap = 0):

    if not gap:
        gap  = (discriminant * discriminant.bit_length()).bit_length()
        gap -= gap % 4

    if not byte_length:
        byte_length = (discriminant.bit_length()//2 - 1) // 8

    h  = hash_to_bytes(message, byte_length)

    a = gmpy2.mpz.from_bytes(h)
    swap_b = a % 2

    a += 3 - (a % 4)
    
    while True:

        if is_possibly_prime(a):

            b = gmpy2.powmod(discriminant % a, (a+1)//4, a)

            if (b**2 - discriminant) % a == 0:
                break

        a += 4

    if b.is_even():
        b = a - b

    if swap_b:
        b = -b

    return ((a,b), 0)

def hybrid_hash(message, discriminant, byte_length = 0, n_small_as = 2, miller_rabin_iterations = 30):

    if not byte_length:
        byte_length = (discriminant.bit_length()//2 - 1) // 8

    small_a_length = byte_length // (2*n_small_as)
    large_a_length = byte_length - n_small_as * small_a_length

    small_a_hash = wesolowski_hash if small_a_length >= 1024 else pseudoprime_hash_with_counter
    large_a_hash = wesolowski_hash if large_a_length >= 1024 else pseudoprime_hash_with_counter

    forms = []
    forms.append(large_a_hash(message, discriminant, large_a_length))
    a = forms[0][0]

    counter = gmpy2.xmpz(0)

    for _ in range(n_small_as):
        while True:

            a_new, b_new = small_a_hash(as_bytes(counter) + message, discriminant, small_a_length)

            if all((a != a_new for (a, _) in forms)):
                break

            counter += 1

        forms.append((a_new, b_new))
        a *= a_new

    b = gmpy2.mpz(sympy.ntheory.modular.crt(*zip(*forms))[0])
    
    if b.is_even():
        b = a - b

    return ((a,b), counter)

hash_functions = {
        "Wesolowski":wesolowski_hash,
        "CLR":clr_hash,
        "pseudoprime-counter":pseudoprime_hash_with_counter,
        "pseudoprime-gap":pseudoprime_hash_with_gap,
        "hybrid":hybrid_hash,
        }
