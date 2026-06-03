#ifndef GMP_INCLUDED
#define GMP_INCLUDED
#include <gmp.h>
#endif

#define SMALL_PRIME_PRODUCTS 30030, 25897, 26381, 29971
#define N_SMALL_PRIME_PRODUCTS 4
#define FIRST_PRIME_OMITTED 53

int isprime(unsigned long t)
{
  unsigned long q, r, d;

  assert (t >= 3 && (t & 1) != 0);
  d = 3;

  do {
      q = t / d;
      r = t - q * d;
      if (q < d)
	return 1;
      d += 2;
  } while (r != 0);

  return 0;
}

int mpz_possible_prime(mpz_t n)
{
  mpz_t r_;
  mpz_init(r_);

  unsigned long r;
  unsigned long spp[] = { SMALL_PRIME_PRODUCTS };

  for (unsigned long i = 0; i < N_SMALL_PRIME_PRODUCTS; i++)
  {
    r = mpz_fdiv_r_ui(r_, n, spp[i]);

    if (gcd(r, spp[i]) != 1) 
    { 
      return 0; 
    }
  }

  {
    unsigned long ln2;
    unsigned long p0, q, p;
    unsigned int primes[15];
    int nprimes;

    nprimes = 0;
    p = 1;
    ln2 = mpz_sizeinbase(n, 2);	
    /* individual primes go into q,
     * prime product into p,
     * primes dividing p into primes[]
     */

    for (q = FIRST_PRIME_OMITTED; q < ln2; q += 2)
      {
      if (isprime(q))
        {
          if (__builtin_umull_overflow(p, q, &p0)) 
            /* p0 overflowed,
             * prime product is now in p
             */
          {
            r = mpz_fdiv_r_ui(r_, n, p);

            while (--nprimes >= 0) { if (r % primes[nprimes] == 0) { return 0; } }

            p = q; // largest prime not checked yet
            nprimes = 0;
          }

          else { p = p0; }

          primes[nprimes++] = q;
        }
      }
  }

  return 1;
}

int mpz_probable_prime(mpz_t n)
{
  return (mpz_possible_prime(n) && mpz_millerrabin(n, 25));
}

int mod_eq_m1(mpz_t x, mpz_t n)
{
  mpz_add_ui(x, x, 1);
  int result = mpz_divisible_p(x, n);
  mpz_sub_ui(x, x, 1);

  return result;
}

int millerrabin_base_2 (mpz_t n)
{
  mpz_t x, exp;
  int k = 0;

  mpz_init_set_ui(x, 2);

  mpz_init(exp);
  mpz_sub_ui(exp, n, 1);

  while (mpz_even_p(exp))
  {
    mpz_tdiv_q_ui(exp, exp, 2);
    k++;
  }

  mpz_powm(x, x, exp, n);

  if (mpz_cmp_ui(x, 1) == 0 || mod_eq_m1(x, n))
  {
    mpz_clear(x);
    return 1;
  }

  for (int i = 1; i < k; i++)
    {
      mpz_powm_ui(x, x, 2, n);

      if (mod_eq_m1(x, n))
      {
        mpz_clear(x);
        return 1;
      }

      if (mpz_cmp_ui(x, 1) <= 0)
      {
        mpz_clear(x);
        return 0;
      }
    }

  mpz_clear(x);
  return 0;
}
