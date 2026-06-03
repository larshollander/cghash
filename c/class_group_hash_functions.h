#ifndef OPENSSL_EVP_INCLUDED
#define OPENSSL_EVP_INCLUDED
#include <openssl/evp.h>
#endif

#ifndef GMP_INCLUDED
#define GMP_INCLUDED
#include <gmp.h>
#endif

#ifndef STDIO_INCLUDED
#define STDIO_INCLUDED
#include <stdio.h>
#endif

#ifndef ASSERT_INCLUDED
#define ASSERT_INCLUDED
#include <assert.h>
#endif

#include "number_theory.h"
#include "primetest.h"

#define WESOLOWSKI       1
#define CLR              2
#define CLR_LF           3

#define PRIME            1
#define PSEUDOPRIME      5
#define PSEUDOPRIME_MR2  7

#define NO_REDUCE_D      1
#define REDUCE_D        11

#define MIN_BLOCK_SIZE 255

// #define DEBUG_HASH
// #define DEBUG_SHAKE

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
  int factor_type;
  int d_mod_a;
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

int wesolowski_pr_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

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
  }
  while (mpz_jacobi(discriminant, a) != 1 || mpz_probable_prime(a) == 0);

  sqrt_mod_3_mod_4(b, discriminant, a);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
  gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

  return *counter;
}

int wesolowski_pp_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

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
    if (mpz_jacobi(discriminant, a) != 1 || mpz_possible_prime(a) == 0) { continue; }

    sqrt_mod_3_mod_4(b, discriminant, a);
    mpz_mul(cmp, b, b);
  }
  while (!mpz_congruent_p(cmp, discriminant, a));

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
  gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

  return *counter;
}

int wesolowski_pp_mr2_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

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
    if (mpz_jacobi(discriminant, a) != 1 || mpz_possible_prime(a) == 0 || millerrabin_base_2(a) == 0) { continue; }

    sqrt_mod_3_mod_4(b, discriminant, a);
    mpz_mul(cmp, b, b);
  }
  while (!mpz_congruent_p(cmp, discriminant, a));

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
  gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

  return *counter;
}

int wesolowski_pr_dma_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

  mpz_t d_mod_a;
  mpz_init(d_mod_a);

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
    mpz_mod(d_mod_a, discriminant, a);
  }
  while (mpz_jacobi(d_mod_a, a) != 1 || mpz_probable_prime(a) == 0);

  sqrt_mod_3_mod_4(b, d_mod_a, a);
  mpz_clear(d_mod_a);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
  gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

  return *counter;
}

int wesolowski_pp_dma_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

  mpz_t d_mod_a;
  mpz_init(d_mod_a);

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
    mpz_mod(d_mod_a, discriminant, a);

    if (mpz_jacobi(d_mod_a, a) != 1 || mpz_possible_prime(a) == 0) { continue; }

    sqrt_mod_3_mod_4(b, d_mod_a, a);
    mpz_mul(cmp, b, b);
  }
  while (!mpz_congruent_p(cmp, d_mod_a, a));

  mpz_clear(d_mod_a);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
  gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

  return *counter;
}

int wesolowski_pp_mr2_dma_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

  mpz_t d_mod_a;
  mpz_init(d_mod_a);

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
    mpz_mod(d_mod_a, discriminant, a);

    if (mpz_jacobi(d_mod_a, a) != 1 || mpz_possible_prime(a) == 0 || millerrabin_base_2(a) == 0) { continue; }

    sqrt_mod_3_mod_4(b, d_mod_a, a);
    mpz_mul(cmp, b, b);
  }
  while (!mpz_congruent_p(cmp, d_mod_a, a));

  mpz_clear(d_mod_a);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
  gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

  return *counter;
}

int clr_lf_pr_hash(S_FORM, S_MSG, int n_small_primes, FILE* fptr)
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
      wesolowski_pr_hash(FORM, MSG, small_prime_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  wesolowski_pr_hash(FORM, MSG, large_prime_size, &counter, fptr);
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

int clr_lf_pp_hash(S_FORM, S_MSG, int n_small_primes, FILE* fptr)
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
      wesolowski_pp_hash(FORM, MSG, small_prime_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  wesolowski_pp_hash(FORM, MSG, large_prime_size, &counter, fptr);
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

int clr_lf_pp_mr2_hash(S_FORM, S_MSG, int n_small_primes, FILE* fptr)
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
      wesolowski_pp_mr2_hash(FORM, MSG, small_prime_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  wesolowski_pp_mr2_hash(FORM, MSG, large_prime_size, &counter, fptr);
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

int clr_sf_pr_hash(S_FORM, S_MSG, int n_small_primes, FILE* fptr)
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
      wesolowski_pr_hash(FORM, MSG, small_prime_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  do 
  {
    wesolowski_pr_hash(FORM, MSG, large_prime_size, &counter, fptr);
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

int clr_sf_pp_hash(S_FORM, S_MSG, int n_small_primes, FILE* fptr)
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
      wesolowski_pp_hash(FORM, MSG, small_prime_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  do 
  {
    wesolowski_pp_hash(FORM, MSG, large_prime_size, &counter, fptr);
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

int clr_sf_pp_mr2_hash(S_FORM, S_MSG, int n_small_primes, FILE* fptr)
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
      wesolowski_pp_mr2_hash(FORM, MSG, small_prime_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  do 
  {
    wesolowski_pp_mr2_hash(FORM, MSG, large_prime_size, &counter, fptr);
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

  switch (spec.algorithm * spec.factor_type * spec.d_mod_a)
  {
    case WESOLOWSKI * PRIME           * NO_REDUCE_D:
      return wesolowski_pr_hash(FORM, MSG, spec.coefficient_size, &counter, fptr);

    case WESOLOWSKI * PRIME           * REDUCE_D   :
      return wesolowski_pr_dma_hash(FORM, MSG, spec.coefficient_size, &counter, fptr);

    case WESOLOWSKI * PSEUDOPRIME     * NO_REDUCE_D:
      return wesolowski_pp_hash(FORM, MSG, spec.coefficient_size, &counter, fptr);

    case WESOLOWSKI * PSEUDOPRIME     * REDUCE_D   :
      return wesolowski_pp_dma_hash(FORM, MSG, spec.coefficient_size, &counter, fptr);

    case WESOLOWSKI * PSEUDOPRIME_MR2 * NO_REDUCE_D:
      return wesolowski_pp_mr2_hash(FORM, MSG, spec.coefficient_size, &counter, fptr);

    case WESOLOWSKI * PSEUDOPRIME_MR2 * REDUCE_D   :
      return wesolowski_pp_mr2_dma_hash(FORM, MSG, spec.coefficient_size, &counter, fptr);

   // case CLR     * PRIME           * NO_REDUCE_D:
   //   return clr_sf_pr_hash(FORM, MSG, spec.n_small_primes, fptr);
   //
   // case CLR     * PRIME           * REDUCE_D   :
   //   return clr_sf_pr_dma_hash(FORM, MSG, spec.n_small_primes, fptr);
   //
   // case CLR     * PSEUDOPRIME     * NO_REDUCE_D:
   //   return clr_sf_pp_hash(FORM, MSG, spec.n_small_primes, fptr);
   //
   // case CLR     * PSEUDOPRIME     * REDUCE_D   :
   //   return clr_sf_pp_dma_hash(FORM, MSG, spec.n_small_primes, fptr);
   //
   // case CLR     * PSEUDOPRIME_MR2 * NO_REDUCE_D:
   //   return clr_sf_pp_mr2_hash(FORM, MSG, spec.n_small_primes, fptr);
   //
   // case CLR     * PSEUDOPRIME_MR2 * REDUCE_D   :
   //   return clr_sf_pp_dma_mr2_hash(FORM, MSG, spec.n_small_primes, fptr);
   //
   // case CLR_LF     * PRIME           * NO_REDUCE_D:
   //   return clr_lf_pr_hash(FORM, MSG, spec.n_small_primes, fptr);
   //
   // case CLR_LF     * PRIME           * REDUCE_D   :
   //   return clr_lf_pr_dma_hash(FORM, MSG, spec.n_small_primes, fptr);
   //
   // case CLR_LF     * PSEUDOPRIME     * NO_REDUCE_D:
   //   return clr_lf_pp_hash(FORM, MSG, spec.n_small_primes, fptr);
   //
   // case CLR_LF     * PSEUDOPRIME     * REDUCE_D   :
   //   return clr_lf_pp_dma_hash(FORM, MSG, spec.n_small_primes, fptr);
   //
   // case CLR_LF     * PSEUDOPRIME_MR2 * NO_REDUCE_D:
   //   return clr_lf_pp_mr2_hash(FORM, MSG, spec.n_small_primes, fptr);
   //
   // case CLR_LF     * PSEUDOPRIME_MR2 * REDUCE_D   :
   //   return clr_lf_pp_dma_mr2_hash(FORM, MSG, spec.n_small_primes, fptr);

    default:
      fprintf(stderr, "unknown class group hash requested\n");
      exit(1);
      return 0;
  }
}

