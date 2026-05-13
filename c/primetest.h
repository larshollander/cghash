#include <assert.h>

#define SMALL_PRIME_PRODUCTS 30030, 25897, 26381, 29971
#define N_SMALL_PRIME_PRODUCTS 4
#define FIRST_PRIME_OMITTED 53

#define swap(A,B) A ^= B; B ^= A; A ^= B

unsigned long int gcd(unsigned long int a, unsigned long int b)
{
  while (b != 0)
  {
    a %= b;
    swap(a,b);
  }

  return a;
}

int isprime (unsigned long int t)
{
  unsigned long int q, r, d;

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

  unsigned long int r;
  unsigned long int spp[] = { SMALL_PRIME_PRODUCTS };

  for (unsigned long int i = 0; i < N_SMALL_PRIME_PRODUCTS; i++)
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

int mpz_mr_base_2(mpz_t n) 
{
  mpz_t base, n_minus_one, exp;

  mpz_init_set_ui(base, 2);

  mpz_init(n_minus_one);
  mpz_sub_ui(n_minus_one, n, 1);

  unsigned int s = 0;
  mpz_init_set(exp, n_minus_one);

  while (mpz_even_p(exp)) 
  {
    mpz_fdiv_q_2exp(exp, exp, 1);
    s++;
  }

  mpz_powm(base, base, exp, n);
  mpz_clear(exp);

  if (mpz_cmp_ui(base, 1) == 0)  {
    mpz_clear(base);
    mpz_clear(n_minus_one);
    return 1;
  }

  for(int i = 0; i < s-1; i++) {
      if (mpz_cmp(base, n_minus_one) == 0) {
        mpz_clear(base);
        mpz_clear(n_minus_one);
        return 1;
      }
      mpz_powm_ui (base, base, 2, n);
  }

  if (mpz_cmp(base, n_minus_one) == 0) {
    mpz_clear(base);
    mpz_clear(n_minus_one);
    return 1;
  }

  mpz_clear(base);
  mpz_clear(n_minus_one);

  return 0;
}

void sample_prime_discriminant(mpz_t discriminant, gmp_randstate_t state, int bit_length)
{
  do
  {
    mpz_urandomb(discriminant, state, bit_length);
    mpz_setbit(discriminant, 0);
    mpz_setbit(discriminant, 1);
  }
  while (mpz_probab_prime_p(discriminant, 25) == 0);

  mpz_neg(discriminant, discriminant);
}
