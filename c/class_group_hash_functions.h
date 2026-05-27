#include <openssl/evp.h>
#include <gmp.h>
#include <stdio.h>
#include <assert.h>

#include "number_theory.h"
#include "primetest.h"

#define WESOLOWSKI           1
#define PSEUDOPRIME_COUNTER  2
#define CLR                  3

#define PRIME                5
#define PSEUDOPRIME          7
#define HYBRID              11

#define HALF_SIZE           13
#define NO_MINIMUM          17

#define MIN_BLOCK_SIZE 255

// #define DEBUG_HASH
// #define DEBUG_SHAKE
#define D_MOD_A

#define FORM a, b, discriminant
#define OUT out, out_len
#define MSG msg, msg_len
#define CTR ctr, ctr_len

#define S_FORM mpz_t a, mpz_t b, mpz_t discriminant
#define S_OUT unsigned char* out, int out_len
#define S_MSG unsigned char* msg, int msg_len
#define S_CTR unsigned char* ctr, int ctr_len

typedef struct
{
  mpz_t a;
  mpz_t b;
  mpz_t discriminant;
} Form;

typedef struct
{
  int algorithm;
  int coefficient_type;
  int large_coefficient;
  int n_small_primes;
  int coefficient_size;
} AlgSpec;

void shake256(S_OUT, S_MSG, S_CTR, FILE* fptr)
{
  OSSL_PARAM params[] = {
    { "xoflen", OSSL_PARAM_UNSIGNED_INTEGER, &out_len, sizeof(out_len), 0 },
    { NULL, 0, NULL, 0, 0 }
  };

  EVP_MD_CTX *ctx = EVP_MD_CTX_new();
  EVP_DigestInit_ex(ctx, EVP_shake256(), NULL);
  EVP_MD_CTX_set_params(ctx, params);

  EVP_DigestUpdate(ctx, MSG);
  EVP_DigestUpdate(ctx, CTR);
  EVP_DigestFinalXOF(ctx, OUT);

#ifdef DEBUG_SHAKE
  fprintf(fptr, "\n");
  for (int i = 0; i < msg_len; i++) { fprintf(fptr, "%02hhx", msg[i]); }
  fprintf(fptr, " ");
  for (int i = 0; i < ctr_len; i++) { fprintf(fptr, "%02hhx", ctr[i]); }
  fprintf(fptr, "\n");
  for (int i = 0; i < out_len; i++) { fprintf(fptr, "%02hhx", out[i]); }
#endif

  EVP_MD_CTX_free(ctx);
}

int pseudoprime_counter_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

#ifdef D_MOD_A
  mpz_t d_mod_a;
  mpz_init(d_mod_a);
#endif

  mpz_t cmp;
  mpz_init(cmp);

#ifdef DEBUG_SHAKE
  fprintf(fptr, "\nmsg length: %d\nctr length: %d", msg_len, counter_len);
#endif

  do
  {
    shake256(hash_value, coefficient_size/8, msg, msg_len, (unsigned char*)counter, counter_len, fptr);
    (*counter)++;
    
    mpz_import(a, coefficient_size/8, 1, 1, 0, 0, hash_value);

    swap_b = mpz_even_p(a);
    mpz_setbit(a, 0);
    mpz_setbit(a, 1);
#ifdef DEBUG_HASH
    gmp_fprintf(fptr, "\n%Zd", a);
#endif
#ifdef D_MOD_A
    mpz_mod(d_mod_a, discriminant, a);

    if (mpz_jacobi(d_mod_a, a) != 1 || mpz_possible_prime(a) == 0) { continue; }

    sqrt_mod_3_mod_4(b, d_mod_a, a);
    mpz_mul(cmp, b, b);
#else
    if (mpz_jacobi(discriminant, a) != 1 || mpz_possible_prime(a) == 0) { continue; }

    sqrt_mod_3_mod_4(b, discriminant, a);
    mpz_mul(cmp, b, b);
#endif
  }
#ifdef D_MOD_A
  while (!mpz_congruent_p(cmp, d_mod_a, a));

  mpz_clear(d_mod_a);
#else
  while (!mpz_congruent_p(cmp, discriminant, a));
#endif

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
  gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

  return *counter;
}

int wesolowski_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

#ifdef D_MOD_A
  mpz_t d_mod_a;
  mpz_init(d_mod_a);
#endif

#ifdef DEBUG_SHAKE
  fprintf(fptr, "\nmsg length: %d\nctr length: %d", msg_len, counter_len);
#endif

  do
  {
    shake256(hash_value, coefficient_size/8, msg, msg_len, (unsigned char*)counter, counter_len, fptr);
    (*counter)++;
    
    mpz_import(a, coefficient_size/8, 1, 1, 0, 0, hash_value);

    swap_b = mpz_even_p(a);
    mpz_setbit(a, 0);
    mpz_setbit(a, 1);
#ifdef DEBUG_HASH
    gmp_fprintf(fptr, "\n%Zd", a);
#endif
#ifdef D_MOD_A
    mpz_mod(d_mod_a, discriminant, a);
#endif
  }
#ifdef D_MOD_A
  while (mpz_jacobi(d_mod_a, a) != 1 || mpz_probable_prime(a) == 0);

  sqrt_mod_3_mod_4(b, d_mod_a, a);
  mpz_clear(d_mod_a);
#else
  while (mpz_jacobi(discriminant, a) != 1 || mpz_probable_prime(a) == 0);

  sqrt_mod_3_mod_4(b, discriminant, a);
#endif


  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
  gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

  return *counter;
}

int clr_pr_lc_hash(S_FORM, S_MSG, int n_small_primes, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int coefficient_size = mpz_sizeinbase(discriminant, 2)/2 - 1;
  int small_prime_size = coefficient_size / (2 * n_small_primes);
  int large_prime_size = coefficient_size - n_small_primes * small_prime_size;

  mpz_t as[n_small_primes + 1];
  mpz_t bs[n_small_primes + 1];

  for (int i = 0; i < n_small_primes; i++)
  {
    do
    {
      wesolowski_hash(FORM, MSG, small_prime_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  wesolowski_hash(FORM, MSG, large_prime_size, &counter, fptr);
  mpz_init_set(as[n_small_primes], a);
  mpz_init_set(bs[n_small_primes], b);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_primes, a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_primes, b);
#endif

  crt(a, b, as, bs, n_small_primes + 1);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }

  return counter;
}

int clr_pp_lc_hash(S_FORM, S_MSG, int n_small_primes, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int coefficient_size = mpz_sizeinbase(discriminant, 2)/2 - 1;
  int small_prime_size = coefficient_size / (2 * n_small_primes);
  int large_prime_size = coefficient_size - n_small_primes * small_prime_size;

  mpz_t as[n_small_primes + 1];
  mpz_t bs[n_small_primes + 1];

  for (int i = 0; i < n_small_primes; i++)
  {
    do
    {
      pseudoprime_counter_hash(FORM, MSG, small_prime_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  pseudoprime_counter_hash(FORM, MSG, large_prime_size, &counter, fptr);
  mpz_init_set(as[n_small_primes], a);
  mpz_init_set(bs[n_small_primes], b);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_primes, a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_primes, b);
#endif

  crt(a, b, as, bs, n_small_primes + 1);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }

  return counter;
}

int clr_hy_lc_hash(S_FORM, S_MSG, int n_small_primes, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int coefficient_size = mpz_sizeinbase(discriminant, 2)/2 - 1;
  int small_prime_size = coefficient_size / (2 * n_small_primes);
  int large_prime_size = coefficient_size - n_small_primes * small_prime_size;

  mpz_t as[n_small_primes + 1];
  mpz_t bs[n_small_primes + 1];

  for (int i = 0; i < n_small_primes; i++)
  {
    do
    {
      pseudoprime_counter_hash(FORM, MSG, small_prime_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  wesolowski_hash(FORM, MSG, large_prime_size, &counter, fptr);
  mpz_init_set(as[n_small_primes], a);
  mpz_init_set(bs[n_small_primes], b);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_primes, a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_primes, b);
#endif

  crt(a, b, as, bs, n_small_primes + 1);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }

  return counter;
}

int clr_pr_sc_hash(S_FORM, S_MSG, int n_small_primes, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int coefficient_size = mpz_sizeinbase(discriminant, 2)/2 - 1;
  int small_prime_size = coefficient_size / (n_small_primes + 1);

  if (coefficient_size < MIN_BLOCK_SIZE) { return 0; }
  if (small_prime_size < MIN_BLOCK_SIZE) 
  { 
    small_prime_size -= (MIN_BLOCK_SIZE - small_prime_size) / n_small_primes;
  } 

  int large_prime_size = coefficient_size - n_small_primes * small_prime_size;

  mpz_t as[n_small_primes + 1];
  mpz_t bs[n_small_primes + 1];

  for (int i = 0; i < n_small_primes; i++)
  {
    do
    {
      wesolowski_hash(FORM, MSG, small_prime_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  do 
  {
    wesolowski_hash(FORM, MSG, large_prime_size, &counter, fptr);
  } while (all_coprime(a, as, n_small_primes) == 0);

  mpz_init_set(as[n_small_primes], a);
  mpz_init_set(bs[n_small_primes], b);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_primes, a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_primes, b);
#endif

  crt(a, b, as, bs, n_small_primes + 1);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }

  return counter;
}

int clr_pp_sc_hash(S_FORM, S_MSG, int n_small_primes, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int coefficient_size = mpz_sizeinbase(discriminant, 2)/2 - 1;
  int small_prime_size = coefficient_size / (n_small_primes + 1);

  if (coefficient_size < MIN_BLOCK_SIZE) { return 0; }
  if (small_prime_size < MIN_BLOCK_SIZE) 
  { 
    small_prime_size -= (MIN_BLOCK_SIZE - small_prime_size) / n_small_primes;
  } 

  int large_prime_size = coefficient_size - n_small_primes * small_prime_size;

  mpz_t as[n_small_primes + 1];
  mpz_t bs[n_small_primes + 1];

  for (int i = 0; i < n_small_primes; i++)
  {
    do
    {
      pseudoprime_counter_hash(FORM, MSG, small_prime_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  do 
  {
    pseudoprime_counter_hash(FORM, MSG, large_prime_size, &counter, fptr);
  } while (all_coprime(a, as, n_small_primes) == 0);

  mpz_init_set(as[n_small_primes], a);
  mpz_init_set(bs[n_small_primes], b);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_primes, a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_primes, b);
#endif

  crt(a, b, as, bs, n_small_primes + 1);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }

  return counter;
}

int hash_to_form(AlgSpec spec, S_FORM, S_MSG, FILE* fptr)
{
  int counter = 0;

  switch (spec.algorithm * spec.coefficient_type * spec.large_coefficient)
  {
    case PSEUDOPRIME_COUNTER:
      return pseudoprime_counter_hash(FORM, MSG, spec.coefficient_size, &counter, fptr);

    case WESOLOWSKI:
      return wesolowski_hash(FORM, MSG, spec.coefficient_size, &counter, fptr);

    case CLR * PRIME * HALF_SIZE:
      return clr_pr_lc_hash(FORM, MSG, spec.n_small_primes, fptr);

    case CLR * PSEUDOPRIME * HALF_SIZE:
      return clr_pp_lc_hash(FORM, MSG, spec.n_small_primes, fptr);

    case CLR * HYBRID * HALF_SIZE:
      return clr_hy_lc_hash(FORM, MSG, spec.n_small_primes, fptr);

    case CLR * PRIME * NO_MINIMUM:
      return clr_pr_sc_hash(FORM, MSG, spec.n_small_primes, fptr);

    case CLR * PSEUDOPRIME * NO_MINIMUM:
      return clr_pp_sc_hash(FORM, MSG, spec.n_small_primes, fptr);

    default:
      fprintf(stderr, "unknown class group hash requested\n");
      return 0;
  }
}

