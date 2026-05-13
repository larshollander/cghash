#define swap(A,B) A ^= B; B ^= A; A ^= B

unsigned long gcd(unsigned long a, unsigned long b)
{
  while (b != 0)
  {
    a %= b;
    swap(a,b);
  }

  return a;
}

void sqrt_mod_3_mod_4(mpz_t out, const mpz_t a, const mpz_t p)
{
  mpz_t base, exp;

  mpz_init_set(base, a);
  mpz_mod(base, base, p);

  mpz_init(exp);
  mpz_add_ui(exp, p, 1);
  mpz_tdiv_q_ui(exp, exp, 4);

  mpz_powm(out, base, exp, p);

  mpz_clear(base);
  mpz_clear(exp);
}

int all_coprime(mpz_t a, mpz_t* as, int n)
{
  mpz_t g, h;
  mpz_init(g);
  mpz_init_set_ui(h, 1);

  for (int i = 0; i < n; i++)
  {
    mpz_gcd(g, a, as[i]);
    
    if (mpz_cmp(g, h)) { return 0; }
  }

  return 1;
}

/**
 * Compute x such that x ≡ a[i] (mod m[i]) for all i
 * Assumes all m[i] are pairwise coprime
 * 
 * @param M Output mpz_t to store the product modulus
 * @param x Output mpz_t to store the result
 * @param a Array of residue values (mpz_t)
 * @param m Array of moduli (mpz_t)
 * @param n Number of congruences
 */
void crt(mpz_t M, mpz_t x, const mpz_t *m, const mpz_t *a, size_t n) 
{
    mpz_t Mi, yi, temp;

    // Initialize temporaries
    mpz_init(Mi);
    mpz_init(yi);
    mpz_init(temp);
    
    // Step 2: Apply CRT formula: x = Σ(a[i] * M[i] * y[i]) mod M
    
    for (size_t i = 0; i < n; i++) 
    {
        // M[i] = M / m[i]
        mpz_div(Mi, M, m[i]);
        
        // y[i] = M[i]^(-1) mod m[i] (modular inverse)
        mpz_invert(yi, Mi, m[i]);
        
        // Add term: a[i] * M[i] * y[i]
        mpz_mul(temp, a[i], Mi);
        mpz_mul(temp, temp, yi);
        mpz_add(x, x, temp);
    }
    
    // Final reduction: result = result mod M
    mpz_mod(x, x, M);
    
    // Clean up temporaries
    mpz_clear(Mi);
    mpz_clear(yi);
    mpz_clear(temp);
}
