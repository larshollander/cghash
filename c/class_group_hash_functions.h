#include <openssl/evp.h>
#include <gmp.h>

#include "primetest.h"
#include "modsqrt.h"

#define EVP_MAX_MD_SIZE 512 //bytes

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

int pseudoprime_counter_hash(mpz_t a, mpz_t b, unsigned char* msg, int msg_len, mpz_t discriminant)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int discriminant_size = mpz_sizeinbase(discriminant, 2);
  int coefficient_size = discriminant_size/2 - 1;
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

    if (mpz_jacobi(discriminant, a) != 1 || !is_possibly_prime(a, b)) { continue; }

    sqrt_mod_3_mod_4(b, discriminant, a);
    mpz_mul(cmp, b, b);
  }
  while (!mpz_congruent_p(cmp, discriminant, a));

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

  return counter;
}

int wesolowski_hash(mpz_t a, mpz_t b, unsigned char* msg, int msg_len, mpz_t discriminant)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int discriminant_size = mpz_sizeinbase(discriminant, 2);
  int coefficient_size = discriminant_size/2 - 1;
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
  while (mpz_jacobi(discriminant, a) != 1 || mpz_probab_prime_p(a, 25) == 0);

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

int pseudoprime_counter_hash_debug(mpz_t a, mpz_t b, unsigned char* msg, int msg_len, mpz_t discriminant, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int discriminant_size = mpz_sizeinbase(discriminant, 2);
  int coefficient_size = discriminant_size/2 - 1;
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

    if (mpz_jacobi(discriminant, a) != 1 || !is_possibly_prime(a, b)) { continue; }

    fprintf(fptr, "*");

    sqrt_mod_3_mod_4(b, discriminant, a);
    mpz_mul(cmp, b, b);
  }
  while (!mpz_congruent_p(cmp, discriminant, a));

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

  fprintf(fptr, "\nhash done\n");
  //for (int i = 0; i < coefficient_size/8; i++) { fprintf(fptr, "%02hhx", hash_value[i]); }

  return counter;
}

int wesolowski_hash_debug(mpz_t a, mpz_t b, unsigned char* msg, int msg_len, mpz_t discriminant, FILE* fptr)
{
  unsigned int counter = 0;
  int counter_len = sizeof(counter);

  int discriminant_size = mpz_sizeinbase(discriminant, 2);
  int coefficient_size = discriminant_size/2 - 1;
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
  while (mpz_jacobi(discriminant, a) != 1 || mpz_probab_prime_p(a, 25) == 0);

  sqrt_mod_3_mod_4(b, discriminant, a);

  if (mpz_even_p(b)) { mpz_sub(b, a, b); }
  if (swap_b) { mpz_neg(b, b); }

  fprintf(fptr, "\n");
  //for (int i = 0; i < coefficient_size/8; i++) { fprintf(fptr, "%02hhx", hash_value[i]); }

  return counter;
}
