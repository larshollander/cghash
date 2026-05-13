#include <openssl/evp.h>
#include <gmp.h>
#include <stdio.h>

#include "primetest.h"
#include "modsqrt.h"

#define EVP_MAX_MD_SIZE 512

#define PSEUDOPRIME_COUNTER 0
#define WESOLOWSKI 1
#define PSEUDOPRIME_COUNTER_DEBUG 2
#define WESOLOWSKI_DEBUG 3

void shake256(unsigned char* out, int out_len, unsigned char* msg, int msg_len, unsigned char* counter, int counter_len)
{
  OSSL_PARAM params[] = {
    { "xoflen", OSSL_PARAM_UNSIGNED_INTEGER, &out_len, sizeof(out_len), 0 },
    { NULL, 0, NULL, 0, 0 }
  };

  EVP_MD_CTX *ctx = EVP_MD_CTX_new();
  EVP_DigestInit_ex(ctx, EVP_shake256(), NULL);
  EVP_MD_CTX_set_params(ctx, params);

  EVP_DigestUpdate(ctx, msg, msg_len);
  EVP_DigestUpdate(ctx, counter, counter_len);
  EVP_DigestFinal_ex(ctx, out, &out_len);

  EVP_MD_CTX_free(ctx);
}

int pseudoprime_counter_hash(mpz_t a, mpz_t b, unsigned char* msg, int msg_len, mpz_t discriminant, int coefficient_size)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

  mpz_t cmp;
  mpz_init_set_ui(cmp, 1);

  do
  {
    shake256(hash_value, coefficient_size/8, msg, msg_len, (unsigned char*)&counter, counter_len);
    counter++;
    
    mpz_import(a, coefficient_size/8, 1, 1, 0, 0, hash_value);

    swap_b = mpz_even_p(a);
    mpz_setbit(a, 0);
    mpz_setbit(a, 1);

    if (mpz_jacobi(discriminant, a) != 1 || mpz_possible_prime(a) == 0) { continue; }

    sqrt_mod_3_mod_4(b, discriminant, a);
    mpz_mul(cmp, b, b);
  }
  while (!mpz_congruent_p(cmp, discriminant, a));

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

  return counter;
}

int wesolowski_hash(mpz_t a, mpz_t b, unsigned char* msg, int msg_len, mpz_t discriminant, int coefficient_size)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

  do
  {
    shake256(hash_value, coefficient_size/8, msg, msg_len, (unsigned char*)&counter, counter_len);
    counter++;
    
    mpz_import(a, coefficient_size/8, 1, 1, 0, 0, hash_value);

    swap_b = mpz_even_p(a);
    mpz_setbit(a, 0);
    mpz_setbit(a, 1);
  }
  while (mpz_jacobi(discriminant, a) != 1 || mpz_probable_prime(a) == 0);

  sqrt_mod_3_mod_4(b, discriminant, a);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

  return counter;
}

void shake256_debug(unsigned char* out, int out_len, unsigned char* msg, int msg_len, unsigned char* ctr, int ctr_len, FILE* fptr)
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
  EVP_DigestFinal_ex(ctx, out, &out_len);

  EVP_MD_CTX_free(ctx);

  fprintf(fptr, "\n");
  for (int i = 0; i < msg_len; i++) { fprintf(fptr, "%02hhx", msg[i]); }
  fprintf(fptr, " ");
  for (int i = 0; i < ctr_len; i++) { fprintf(fptr, "%02hhx", ctr[i]); }
  fprintf(fptr, "\n");
  for (int i = 0; i < out_len; i++) { fprintf(fptr, "%02hhx", out[i]); }
}

int pseudoprime_counter_hash_debug(mpz_t a, mpz_t b, unsigned char* msg, int msg_len, mpz_t discriminant, int coefficient_size, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

  mpz_t cmp;
  mpz_init_set_ui(cmp, 1);

  do
  {
    shake256_debug(hash_value, coefficient_size/8, msg, msg_len, (unsigned char*)&counter, counter_len, fptr);
    counter++;
    
    mpz_import(a, coefficient_size/8, 1, 1, 0, 0, hash_value);

    swap_b = mpz_even_p(a);
    mpz_setbit(a, 0);
    mpz_setbit(a, 1);
    gmp_fprintf(fptr, "\n%Zd", a);

    if (mpz_jacobi(discriminant, a) != 1 || mpz_possible_prime(a) == 0) { continue; }

    sqrt_mod_3_mod_4(b, discriminant, a);
    mpz_mul(cmp, b, b);
  }
  while (!mpz_congruent_p(cmp, discriminant, a));

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

  fprintf(fptr, "\nhash done\n");

  return counter;
}

int wesolowski_hash_debug(mpz_t a, mpz_t b, unsigned char* msg, int msg_len, mpz_t discriminant, int coefficient_size, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  unsigned char hash_value[coefficient_size];

  int swap_b = 0;

  do
  {
    shake256_debug(hash_value, coefficient_size/8, msg, msg_len, (unsigned char*)&counter, counter_len, fptr);
    counter++;
    
    mpz_import(a, coefficient_size/8, 1, 1, 0, 0, hash_value);

    swap_b = mpz_even_p(a);
    mpz_setbit(a, 0);
    mpz_setbit(a, 1);
    gmp_fprintf(fptr, "\n%Zd", a);
  }
  while (mpz_jacobi(discriminant, a) != 1 || mpz_probable_prime(a) == 0);

  sqrt_mod_3_mod_4(b, discriminant, a);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

  fprintf(fptr, "\nhash done\n");

  return counter;
}

int hash_to_form(int cg_hash_algorithm, mpz_t a, mpz_t b, unsigned char* msg, int msg_len, mpz_t discriminant, FILE* fptr)
{
  int discriminant_size = mpz_sizeinbase(discriminant, 2);
  int coefficient_size;

  switch (cg_hash_algorithm)
  {
    case PSEUDOPRIME_COUNTER:
      coefficient_size = discriminant_size/2 - 1;
      return pseudoprime_counter_hash(a, b, msg, msg_len, discriminant, coefficient_size);
    case WESOLOWSKI:
      coefficient_size = discriminant_size/2 - 1;
      return wesolowski_hash(a, b, msg, msg_len, discriminant, coefficient_size);
    case PSEUDOPRIME_COUNTER_DEBUG:
      coefficient_size = discriminant_size/2 - 1;
      return pseudoprime_counter_hash_debug(a, b, msg, msg_len, discriminant, coefficient_size, fptr);
    case WESOLOWSKI_DEBUG:
      coefficient_size = discriminant_size/2 - 1;
      return wesolowski_hash_debug(a, b, msg, msg_len, discriminant, coefficient_size, fptr);
    default:
      fprintf(stderr, "unknown class group hash requested\n");
      return 0;
  }
}

