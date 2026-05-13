#include <assert.h>

#include "number_theory.h"

#define SMALL_PRIME_PRODUCTS 30030, 25897, 26381, 29971
#define N_SMALL_PRIME_PRODUCTS 4
#define FIRST_PRIME_OMITTED 53

int isprime (unsigned long t)
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
    ln2 = mpz_sizeinbase (n, 2);	
    /* individual primes go into q,
     * prime product into p,
     * primes dividing p into primes[]
     */

    for (q = FIRST_PRIME_OMITTED; q < ln2; q += 2)
      {
      if (isprime (q))
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
  if (mpz_possible_prime(n)) { return mpz_millerrabin(n, 25); }
}

int mr_base_2 (mpz_t n, mpz_t q, unsigned long k)
{
  mpz_t x, y;
  mpz_init_set_ui(x, 2);

  mpz_init(y);
  mpz_powm (y, x, q, n);

  if (mpz_cmp_ui (y, 1L) == 0 || mod_eq_m1 (y, n)) { return 1; }

  for (unsigned long i = 1; i < k; i++)
    {
      mpz_powm_ui (y, y, 2L, n);
      if (mod_eq_m1 (y, n)) { return 1; }
      if (mpz_cmp_ui (y, 1L) <= 0) { return 0; }
    }
  return 0;
}

void sample_prime_discriminant(mpz_t discriminant, gmp_randstate_t state, int bit_length)
{
  do
  {
    mpz_urandomb(discriminant, state, bit_length);
    mpz_setbit(discriminant, 0);
    mpz_setbit(discriminant, 1);
    mpz_setbit(discriminant, bit_length-1);
  }
  while (mpz_probab_prime_p(discriminant, 25) == 0);

  mpz_neg(discriminant, discriminant);
}
