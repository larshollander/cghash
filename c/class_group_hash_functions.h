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

#define MIN_BLOCK_SIZE 384
#define WES_BEATS_PP 640

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

void shake256(unsigned char* out, int out_len, unsigned char* msg, int msg_len, unsigned char* ctr, int ctr_len, FILE* fptr)
{
  OSSL_PARAM params[] = {
    { "xoflen", OSSL_PARAM_UNSIGNED_INTEGER, &out_len, sizeof(out_len), 0 },
    { NULL, 0, NULL, 0, 0 }
  };

  EVP_MD_CTX *ctx = EVP_MD_CTX_new();
  EVP_DigestInit_ex(ctx, EVP_shake256(), NULL);
  EVP_MD_CTX_set_params(ctx, params);

  EVP_DigestUpdate(ctx, msg, msg_len);
  EVP_DigestUpdate(ctx, ctr, ctr_len);
  EVP_DigestFinalXOF(ctx, out, out_len);

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

int pseudoprime_counter_hash(Form f, unsigned char* msg, int msg_len, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

  mpz_t cmp;
  mpz_init_set_ui(cmp, 1);

  do
  {
    shake256(hash_value, coefficient_size/8, msg, msg_len, (unsigned char*)counter, counter_len, fptr);
    (*counter)++;
    
    mpz_import(f.a, coefficient_size/8, 1, 1, 0, 0, hash_value);

    swap_b = mpz_even_p(f.a);
    mpz_setbit(f.a, 0);
    mpz_setbit(f.a, 1);
#ifdef DEBUG_HASH
    gmp_fprintf(fptr, "\n%Zd", f.a);
#endif

    if (mpz_jacobi(f.discriminant, f.a) != 1 || mpz_possible_prime(f.a) == 0) { continue; }

    sqrt_mod_3_mod_4(f.b, f.discriminant, f.a);
    mpz_mul(cmp, f.b, f.b);
  }
  while (!mpz_congruent_p(cmp, f.discriminant, f.a));

  if (mpz_even_p(f.b)) { mpz_sub(f.b, f.a, f.b); }
  if (swap_b) { mpz_neg(f.b, f.b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
#endif

  return *counter;
}

int wesolowski_hash(Form f, unsigned char* msg, int msg_len, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

  do
  {
    shake256(hash_value, coefficient_size/8, msg, msg_len, (unsigned char*)counter, counter_len, fptr);
    (*counter)++;
    
    mpz_import(f.a, coefficient_size/8, 1, 1, 0, 0, hash_value);

    swap_b = mpz_even_p(f.a);
    mpz_setbit(f.a, 0);
    mpz_setbit(f.a, 1);
#ifdef DEBUG_HASH
    gmp_fprintf(fptr, "\n%Zd", f.a);
#endif
  }
  while (mpz_jacobi(f.discriminant, f.a) != 1 || mpz_probable_prime(f.a) == 0);

  sqrt_mod_3_mod_4(f.b, f.discriminant, f.a);

  if (mpz_even_p(f.b)) { mpz_sub(f.b, f.a, f.b); }
  if (swap_b) { mpz_neg(f.b, f.b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
#endif

  return *counter;
}

int fastest_small_hash(Form f, unsigned char* msg, int msg_len, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  if (coefficient_size >= WES_BEATS_PP)
  {
    return pseudoprime_counter_hash(f, msg, msg_len, coefficient_size, counter, fptr);
  }
  else
  {
    return wesolowski_hash(f, msg, msg_len, coefficient_size, counter, fptr);
  }
}

int clr_pr_lc_hash(Form f, unsigned char* msg, int msg_len, int n_small_primes, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int coefficient_size = mpz_sizeinbase(f.discriminant, 2)/2 - 1;
  int small_prime_size = coefficient_size / (2 * n_small_primes);
  int large_prime_size = coefficient_size - n_small_primes * small_prime_size;

  mpz_t as[n_small_primes + 1];
  mpz_t bs[n_small_primes + 1];

  for (int i = 0; i < n_small_primes; i++)
  {
    do
    {
      wesolowski_hash(f, msg, msg_len, small_prime_size, &counter, fptr);
    } while (all_coprime(f.a, as, i) == 0);

    mpz_init_set(as[i], f.a);
    mpz_init_set(bs[i], f.b);
  }

  wesolowski_hash(f, msg, msg_len, large_prime_size, &counter, fptr);
  mpz_init_set(as[n_small_primes], f.a);
  mpz_init_set(bs[n_small_primes], f.b);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_primes, f.a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_primes, f.b);
#endif

  crt(f.a, f.b, as, bs, n_small_primes + 1);

  if (mpz_even_p(f.b)) { mpz_sub(f.b, f.a, f.b); }

  return counter;
}

int clr_pp_lc_hash(Form f, unsigned char* msg, int msg_len, int n_small_primes, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int coefficient_size = mpz_sizeinbase(f.discriminant, 2)/2 - 1;
  int small_prime_size = coefficient_size / (2 * n_small_primes);
  int large_prime_size = coefficient_size - n_small_primes * small_prime_size;

  mpz_t as[n_small_primes + 1];
  mpz_t bs[n_small_primes + 1];

  for (int i = 0; i < n_small_primes; i++)
  {
    do
    {
      pseudoprime_counter_hash(f, msg, msg_len, small_prime_size, &counter, fptr);
    } while (all_coprime(f.a, as, i) == 0);

    mpz_init_set(as[i], f.a);
    mpz_init_set(bs[i], f.b);
  }

  pseudoprime_counter_hash(f, msg, msg_len, large_prime_size, &counter, fptr);
  mpz_init_set(as[n_small_primes], f.a);
  mpz_init_set(bs[n_small_primes], f.b);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_primes, f.a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_primes, f.b);
#endif

  crt(f.a, f.b, as, bs, n_small_primes + 1);

  if (mpz_even_p(f.b)) { mpz_sub(f.b, f.a, f.b); }

  return counter;
}

int clr_hy_lc_hash(Form f, unsigned char* msg, int msg_len, int n_small_primes, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int coefficient_size = mpz_sizeinbase(f.discriminant, 2)/2 - 1;
  int small_prime_size = coefficient_size / (2 * n_small_primes);
  int large_prime_size = coefficient_size - n_small_primes * small_prime_size;

  mpz_t as[n_small_primes + 1];
  mpz_t bs[n_small_primes + 1];

  for (int i = 0; i < n_small_primes; i++)
  {
    do
    {
      fastest_small_hash(f, msg, msg_len, small_prime_size, &counter, fptr);
    } while (all_coprime(f.a, as, i) == 0);

    mpz_init_set(as[i], f.a);
    mpz_init_set(bs[i], f.b);
  }

  fastest_small_hash(f, msg, msg_len, large_prime_size, &counter, fptr);
  mpz_init_set(as[n_small_primes], f.a);
  mpz_init_set(bs[n_small_primes], f.b);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_primes, f.a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_primes, f.b);
#endif

  crt(f.a, f.b, as, bs, n_small_primes + 1);

  if (mpz_even_p(f.b)) { mpz_sub(f.b, f.a, f.b); }

  return counter;
}

int clr_pr_sc_hash(Form f, unsigned char* msg, int msg_len, int n_small_primes, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int coefficient_size = mpz_sizeinbase(f.discriminant, 2)/2 - 1;
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
      wesolowski_hash(f, msg, msg_len, small_prime_size, &counter, fptr);
    } while (all_coprime(f.a, as, i) == 0);

    mpz_init_set(as[i], f.a);
    mpz_init_set(bs[i], f.b);
  }

  do 
  {
    wesolowski_hash(f, msg, msg_len, large_prime_size, &counter, fptr);
  } while (all_coprime(f.a, as, n_small_primes) == 0);

  mpz_init_set(as[n_small_primes], f.a);
  mpz_init_set(bs[n_small_primes], f.b);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_primes, f.a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_primes, f.b);
#endif

  crt(f.a, f.b, as, bs, n_small_primes + 1);

  if (mpz_even_p(f.b)) { mpz_sub(f.b, f.a, f.b); }

  return counter;
}

int clr_pp_sc_hash(Form f, unsigned char* msg, int msg_len, int n_small_primes, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int coefficient_size = mpz_sizeinbase(f.discriminant, 2)/2 - 1;
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
      pseudoprime_counter_hash(f, msg, msg_len, small_prime_size, &counter, fptr);
    } while (all_coprime(f.a, as, i) == 0);

    mpz_init_set(as[i], f.a);
    mpz_init_set(bs[i], f.b);
  }

  do 
  {
    pseudoprime_counter_hash(f, msg, msg_len, large_prime_size, &counter, fptr);
  } while (all_coprime(f.a, as, n_small_primes) == 0);

  mpz_init_set(as[n_small_primes], f.a);
  mpz_init_set(bs[n_small_primes], f.b);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_primes, f.a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_primes, f.b);
#endif

  crt(f.a, f.b, as, bs, n_small_primes + 1);

  if (mpz_even_p(f.b)) { mpz_sub(f.b, f.a, f.b); }

  return counter;
}

int clr_hy_sc_hash(Form f, unsigned char* msg, int msg_len, int n_small_primes, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int coefficient_size = mpz_sizeinbase(f.discriminant, 2)/2 - 1;
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
      fastest_small_hash(f, msg, msg_len, small_prime_size, &counter, fptr);
    } while (all_coprime(f.a, as, i) == 0);

    mpz_init_set(as[i], f.a);
    mpz_init_set(bs[i], f.b);
  }

  do 
  {
    fastest_small_hash(f, msg, msg_len, large_prime_size, &counter, fptr);
  } while (all_coprime(f.a, as, n_small_primes) == 0);

  mpz_init_set(as[n_small_primes], f.a);
  mpz_init_set(bs[n_small_primes], f.b);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_primes, f.a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_primes, f.b);
#endif

  crt(f.a, f.b, as, bs, n_small_primes + 1);

  if (mpz_even_p(f.b)) { mpz_sub(f.b, f.a, f.b); }

  return counter;
}

int hash_to_form(AlgSpec spec, Form f, unsigned char* msg, int msg_len, FILE* fptr)
{
  int counter = 0;

  switch (spec.algorithm * spec.coefficient_type * spec.large_coefficient)
  {
    case PSEUDOPRIME_COUNTER:
      return pseudoprime_counter_hash(f, msg, msg_len, spec.coefficient_size, &counter, fptr);

    case WESOLOWSKI:
      return wesolowski_hash(f, msg, msg_len, spec.coefficient_size, &counter, fptr);

    case CLR * PRIME * HALF_SIZE:
      return clr_pr_lc_hash(f, msg, msg_len, spec.n_small_primes, fptr);

    case CLR * PSEUDOPRIME * HALF_SIZE:
      return clr_hy_lc_hash(f, msg, msg_len, spec.n_small_primes, fptr);

    case CLR * HYBRID * HALF_SIZE:
      return clr_pp_sc_hash(f, msg, msg_len, spec.n_small_primes, fptr);

    case CLR * PRIME * NO_MINIMUM:
      return clr_pp_lc_hash(f, msg, msg_len, spec.n_small_primes, fptr);

    case CLR * PSEUDOPRIME * NO_MINIMUM:
      return clr_pr_sc_hash(f, msg, msg_len, spec.n_small_primes, fptr);

    case CLR * HYBRID * NO_MINIMUM:
      return clr_hy_sc_hash(f, msg, msg_len, spec.n_small_primes, fptr);

    default:
      fprintf(stderr, "unknown class group hash requested\n");
      return 0;
  }
}

