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

#define SINGLE_FACTOR    1
#define CLR_LF           2
#define CLR_SF           3

#define PRIME            1
#define PSEUDOPRIME      5
#define PSEUDOPRIME_MR2  7

#define NO_REDUCE_D      1
#define REDUCE_D        11

#define MIN_LF_SIZE     16
#define MIN_FACTOR_SIZE 16

#define DEBUG_HASH
// #define DEBUG_SHAKE

#define FORM a, b, discriminant
#define OUT out, out_len
#define MSG msg, msg_len
#define CTR ctr, ctr_len

#define S_FORM mpz_t a, mpz_t b, mpz_t discriminant
#define S_OUT unsigned char* out, int out_len
#define S_MSG unsigned char* msg, int msg_len
#define S_CTR unsigned char* ctr, int ctr_len

#define GET_A \
  shake256(hash_value, coefficient_size/8 + 1, msg, msg_len, (unsigned char*)counter, counter_len, fptr);\
  (*counter)++;\
  mpz_import(a, coefficient_size/8 + 1, 1, 1, 0, 0, hash_value);\
  for (int i = coefficient_size; i < coefficient_size + 8; i++) { mpz_clrbit(a, i); }\
  swap_b = mpz_even_p(a);\
  mpz_setbit(a, 0);\
  mpz_setbit(a, 1);

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
  int n_small_factors;
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

  EVP_DigestUpdate(ctx, (unsigned char*) &msg_len, sizeof(msg_len));
  EVP_DigestUpdate(ctx, MSG);
  EVP_DigestUpdate(ctx, CTR);
  EVP_DigestFinalXOF(ctx, OUT);

#ifdef DEBUG_SHAKE
  for (int i = 0; i < msg_len; i++) { fprintf(fptr, "%02hhx", msg[i]); }
  fprintf(fptr, " ");
  for (int i = 0; i < ctr_len; i++) { fprintf(fptr, "%02hhx", ctr[i]); }
  fprintf(fptr, "\n");
  for (int i = 0; i < out_len; i++) { fprintf(fptr, "%02hhx", out[i]); }
  fprintf(fptr, "\n");
#endif

  EVP_MD_CTX_free(ctx);
}

int pr_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size/8 + 1];

  int swap_b = 0;

#ifdef DEBUG_SHAKE
  fprintf(fptr, "\nmsg length: %d\nctr length: %d\n", msg_len, counter_len);
#endif

  do
  {
    GET_A
#ifdef DEBUG_HASH
    gmp_fprintf(fptr, "\n%Zd", a);
#endif
  }
  while (mpz_possible_prime(a) == 0 || mpz_jacobi(discriminant, a) != 1 || mpz_millerrabin(a, 25) == 0);

  sqrt_mod_3_mod_4(b, discriminant, a);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
  gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

  return *counter;
}

int pp_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size/8 + 1];

  int swap_b = 0;

  mpz_t cmp;
  mpz_init(cmp);

#ifdef DEBUG_SHAKE
  fprintf(fptr, "\nmsg length: %d\nctr length: %d\n", msg_len, counter_len);
#endif

  do
  {
    GET_A
#ifdef DEBUG_HASH
    gmp_fprintf(fptr, "\n%Zd", a);
#endif
    if (mpz_possible_prime(a) == 0 || mpz_jacobi(discriminant, a) != 1) { continue; }

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

int pp_mr2_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size/8 + 1];

  int swap_b = 0;

  mpz_t cmp;
  mpz_init(cmp);

#ifdef DEBUG_SHAKE
  fprintf(fptr, "\nmsg length: %d\nctr length: %d\n", msg_len, counter_len);
#endif

  do
  {
    GET_A
#ifdef DEBUG_HASH
    gmp_fprintf(fptr, "\n%Zd", a);
#endif
    if (mpz_possible_prime(a) == 0 || mpz_jacobi(discriminant, a) != 1 || millerrabin_base_2(a) == 0) { continue; }

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

int pr_dma_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size/8 + 1];

  int continue_used;
  int swap_b = 0;

  mpz_t d_mod_a;
  mpz_init(d_mod_a);

#ifdef DEBUG_SHAKE
  fprintf(fptr, "\nmsg length: %d\nctr length: %d\n", msg_len, counter_len);
#endif

  do
  {
    GET_A
#ifdef DEBUG_HASH
    gmp_fprintf(fptr, "\n%Zd", a);
#endif

    continue_used = 0;

    if (mpz_possible_prime(a) == 0) { continue_used = 1; continue; }

    mpz_mod(d_mod_a, discriminant, a);

    if (mpz_jacobi(d_mod_a, a) != 1) { continue_used = 1; continue; }
  }
  while (continue_used || mpz_millerrabin(a, 25) == 0);

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

int pp_dma_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size/8 + 1];

  int continue_used;
  int swap_b = 0;

  mpz_t d_mod_a;
  mpz_init(d_mod_a);

  mpz_t cmp;
  mpz_init(cmp);

#ifdef DEBUG_SHAKE
  fprintf(fptr, "\nmsg length: %d\nctr length: %d\n", msg_len, counter_len);
#endif

  do
  {
    GET_A
#ifdef DEBUG_HASH
    gmp_fprintf(fptr, "\n%Zd", a);
#endif

    continue_used = 0;

    if (mpz_possible_prime(a) == 0) { continue_used = 1; continue; }

    mpz_mod(d_mod_a, discriminant, a); 

    if (mpz_jacobi(d_mod_a, a) != 1) { continue_used = 1; continue; }

    sqrt_mod_3_mod_4(b, d_mod_a, a);
    mpz_mul(cmp, b, b);
  }
  while (continue_used || !mpz_congruent_p(cmp, d_mod_a, a));

  mpz_clear(d_mod_a);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
  gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

  return *counter;
}

int pp_mr2_dma_hash(S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  int counter_len = sizeof(*counter);

  unsigned char hash_value[coefficient_size/8 + 1];

  int continue_used;
  int swap_b = 0;

  mpz_t d_mod_a;
  mpz_init(d_mod_a);

  mpz_t cmp;
  mpz_init(cmp);

#ifdef DEBUG_SHAKE
  fprintf(fptr, "\nmsg length: %d\nctr length: %d\n", msg_len, counter_len);
#endif

  do
  {
    GET_A
#ifdef DEBUG_HASH
    gmp_fprintf(fptr, "\n%Zd", a);
#endif

    continue_used = 0;

    if (mpz_possible_prime(a) == 0) { continue_used = 1; continue; } 

    mpz_mod(d_mod_a, discriminant, a); 

    if (mpz_jacobi(d_mod_a, a) != 1) { continue_used = 1; continue; }
    if (millerrabin_base_2(a) == 0) { continue_used = 1; continue; }

    sqrt_mod_3_mod_4(b, d_mod_a, a);
    mpz_mul(cmp, b, b);
  }
  while (continue_used || !mpz_congruent_p(cmp, d_mod_a, a));

  mpz_clear(d_mod_a);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

#ifdef DEBUG_HASH
  fprintf(fptr, "\nhash done\n");
  gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

  return *counter;
}

int single_factor_hash(AlgSpec spec, S_FORM, S_MSG, int coefficient_size, unsigned int* counter, FILE* fptr)
{
  switch (spec.factor_type * spec.d_mod_a)
  {
    case PRIME                      :
      return pr_hash(FORM, MSG, coefficient_size, counter, fptr);

    case PRIME * REDUCE_D           :
      return pr_dma_hash(FORM, MSG, coefficient_size, counter, fptr);

    case PSEUDOPRIME                :
      return pp_hash(FORM, MSG, coefficient_size, counter, fptr);

    case PSEUDOPRIME * REDUCE_D     :
      return pp_dma_hash(FORM, MSG, coefficient_size, counter, fptr);

    case PSEUDOPRIME_MR2            :
      return pp_mr2_hash(FORM, MSG, coefficient_size, counter, fptr);

    case PSEUDOPRIME_MR2 * REDUCE_D :
      return pp_dma_hash(FORM, MSG, coefficient_size, counter, fptr);

    default:
      fprintf(stderr, "unknown class group hash requested\n");
      exit(1);
  }
}

int clr_lf_hash(AlgSpec spec, S_FORM, S_MSG, int n_small_factors, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int small_factor_size = spec.coefficient_size / (2 * n_small_factors);
  int large_factor_size = spec.coefficient_size - n_small_factors * small_factor_size;

  if (small_factor_size < MIN_FACTOR_SIZE) { return 0; }

  mpz_t as[n_small_factors + 1];
  mpz_t bs[n_small_factors + 1];

  for (int i = 0; i < n_small_factors; i++)
  {
    do
    {
      single_factor_hash(spec, FORM, MSG, small_factor_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  single_factor_hash(spec, FORM, MSG, large_factor_size, &counter, fptr);

  mpz_init_set(as[n_small_factors], a);
  mpz_init_set(bs[n_small_factors], b);

  for (int i = 0; i < n_small_factors; i++) { mpz_mul(a, a, as[i]); }

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "  as[%d] = %Zd\n", n_small_factors, a);
  gmp_fprintf(fptr, "  bs[%d] = %Zd\n", n_small_factors, b);
#endif

  crt(a, b, as, bs, n_small_factors + 1);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }

  mpz_mul_si(b, b, mpz_sgn(bs[0]));

  return counter;
}

int clr_sf_hash(AlgSpec spec, S_FORM, S_MSG, int n_small_factors, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int small_factor_size = spec.coefficient_size / (n_small_factors + 1);

  if (small_factor_size < MIN_LF_SIZE) 
  { 
    small_factor_size -= (MIN_LF_SIZE - small_factor_size) / n_small_factors;
  } 
  if (small_factor_size < MIN_FACTOR_SIZE) { return 0; }

  int large_factor_size = spec.coefficient_size - n_small_factors * small_factor_size;

  mpz_t as[n_small_factors + 1];
  mpz_t bs[n_small_factors + 1];

  for (int i = 0; i < n_small_factors; i++)
  {
    do
    {
      single_factor_hash(spec, FORM, MSG, small_factor_size, &counter, fptr);
    } while (all_coprime(a, as, i) == 0);

    mpz_init_set(as[i], a);
    mpz_init_set(bs[i], b);
  }

  do 
  {
    single_factor_hash(spec, FORM, MSG, large_factor_size, &counter, fptr);
  } while (all_coprime(a, as, n_small_factors) == 0);

  mpz_init_set(as[n_small_factors], a);
  mpz_init_set(bs[n_small_factors], b);

  for (int i = 0; i < n_small_factors; i++) { mpz_mul(a, a, as[i]); }

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "\n");

  for (int i = 0; i <= n_small_factors; i++)
  {
    gmp_fprintf(fptr, "  as[%d] = %Zd\n", i, as[i]);
    gmp_fprintf(fptr, "  bs[%d] = %Zd\n", i, bs[i]);
  }
#endif

  crt(a, b, as, bs, n_small_factors + 1);

#ifdef DEBUG_HASH
  gmp_fprintf(fptr, "\nhash done\n");

  gmp_fprintf(fptr, "  a = %Zd\n", a);
  gmp_fprintf(fptr, "  b = %Zd\n", b);
#endif

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }

  mpz_mul_si(b, b, mpz_sgn(bs[0]));

  return counter;
}

int hash_to_form(AlgSpec spec, S_FORM, S_MSG, FILE* fptr)
{
  int counter = 0;

  if (spec.algorithm == SINGLE_FACTOR)
  {
    return single_factor_hash(spec, FORM, MSG, spec.coefficient_size, &counter, fptr);
  }
  else if (spec.algorithm == CLR_LF)
  {
    return clr_lf_hash(spec, FORM, MSG, spec.n_small_factors, fptr);
  }
  else if (spec.algorithm == CLR_SF)
  {
    return clr_sf_hash(spec, FORM, MSG, spec.n_small_factors, fptr);
  }
  else
  {
    fprintf(stderr, "unknown class group hash requested\n");
    exit(1);
  }
}

