#define SMALL_PRIME_PRODUCTS 30030, 25897, 26381, 29971
#define N_SMALL_PRIME_PRODUCTS 4

#define swap(A,B) A ^= B; B ^= A; A ^= B

unsigned int gcd(unsigned int a, unsigned int b)
{
  while (b != 0)
  {
    a %= b;
    swap(a,b);
  }

  return a;
}

int is_possibly_prime(mpz_t a, mpz_t r_)
{
  unsigned int r;
  unsigned int spp[] = { SMALL_PRIME_PRODUCTS };

  for (unsigned int i = 0; i < N_SMALL_PRIME_PRODUCTS; i++)
  {
    r = mpz_fdiv_r_ui(r_, a, spp[i]);
    if (gcd(r, spp[i]) != 1) { return 0; }
  }

  return 1;
}
