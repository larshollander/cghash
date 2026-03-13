import gmpy2

rng_state = gmpy2.random_state()

def tonelli_shanks(a, p):

    e, q = gmpy2.xmpz(0), p-1

    while q.is_even():
        q //= 2
        e += 1

    n_bits = gmpy2.bit_length(p) - 1

    while True:

        n = gmpy2.mpz_urandomb(rng_state, n_bits)

        if gmpy2.legendre(n, p) == -1:
            break

    z = pow(n, q, p)
    y = z
    r = e
    x = pow(a, (q-1)//2, p)
    b = a*(x**2 % p) % p
    x = a*x % p
    

    while b != 1:

        m = 1

        while pow(b, 2**m, p) != 1:
            m += 1

        if m == r:
            return None

        t = pow(y, 2**(r-m-1), p)
        y = pow(t, 2, p)
        r = m
        x = x*t % p
        b = b*y % p

    return x
