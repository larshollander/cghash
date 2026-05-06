#include <gmp.h>

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
