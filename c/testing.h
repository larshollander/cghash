#ifndef GMP_INCLUDED
#define GMP_INCLUDED
#include <gmp.h>
#endif

#ifndef STDIO_INCLUDED
#define STDIO_INCLUDED
#include <stdio.h>
#endif

#ifndef TIME_INCLUDED
#define TIME_INCLUDED
#include <time.h>
#endif

#include "class_group_hash_functions.h"

// #define LOG_COUNTERS
// #define LOG_COEFFICIENTS

void sample_prime_discriminant(mpz_t discriminant, gmp_randstate_t state, int bit_length)
{
  do
  {
    mpz_urandomb(discriminant, state, bit_length);
    mpz_setbit(discriminant, 0);
    mpz_setbit(discriminant, 1);
    mpz_setbit(discriminant, bit_length-1);
  }
  while (mpz_probab_prime_p(discriminant, 25) == 0);

  mpz_neg(discriminant, discriminant);
}

void sample_discriminants(mpz_t* discriminants, int n_discriminants, int discriminant_size)
{
  gmp_randstate_t state;
  gmp_randinit_default(state);

  printf("\rdiscriminants sampled: 0/%d\n\033[A", n_discriminants);

  for (int i = 0; i < n_discriminants; i++)
  {
    mpz_init(discriminants[i]);
    sample_prime_discriminant(discriminants[i], state, discriminant_size);
    printf("\rdiscriminants sampled: %d/%d\n\033[A", i + 1, n_discriminants);
  }

  printf("\n");
}

void test_discriminant(AlgSpec spec, mpz_t discriminant, int max_input, long* sum_counters, double* sum_times, int n_hashes, int n_done, FILE* fptr)
{
  mpz_t a, b;
  mpz_init(a);
  mpz_init(b);

  gmp_fprintf(fptr, "D = %Zd\n", discriminant);

  double begin = (double)clock();

  for (unsigned int input = 0; input < max_input; input++)
  {
    *sum_counters += hash_to_form(spec, FORM, (unsigned char*)&input, sizeof(input), fptr);

#ifdef LOG_COUNTERS
    fprintf(fptr, "counter: %ld\n", *sum_counters);
#endif

#ifdef LOG_COEFFICIENTS
    gmp_fprintf(fptr, "a = %Zd\nb = %Zd\n", a, b);
#endif

    printf("\rvalues hashed: %d/%d\n\033[A", n_done + input + 1, n_hashes);
  }

  double end = (double)clock();
  *sum_times += (end - begin) / CLOCKS_PER_SEC;

  mpz_clear(a);
  mpz_clear(b);
}

void set_names(char* filename, char* hashname, AlgSpec* spec)
{
  switch (spec->algorithm)
  {
    case SINGLE_FACTOR:
      strncpy(filename, "logs/", 5);
      strncpy(hashname, "", 0);
      filename += 5;
      hashname += 0;
      break;

    case CLR_LF:
      filename += snprintf(filename, 15, "logs/clr-%d-lf-", spec->n_small_factors);
      hashname += snprintf(hashname, 10, "clr-%d-lf-", spec->n_small_factors);
      break;

    case CLR_SF:
      filename += snprintf(filename, 15, "logs/clr-%d-sf-", spec->n_small_factors);
      hashname += snprintf(hashname, 10, "clr-%d-sf-", spec->n_small_factors);
      break;

    default:
      break;
  }

  switch (spec->factor_type)
  {
    case PRIME:
      strncpy(filename, "pr", 2);
      strncpy(hashname, "pr", 2);
      filename += 2;
      hashname += 2;
      break;

    case PSEUDOPRIME:
      strncpy(filename, "pp", 2);
      strncpy(hashname, "pp", 2);
      filename += 2;
      hashname += 2;
      break;

    case PSEUDOPRIME_MR2:
      strncpy(filename, "pp-mr2", 6);
      strncpy(hashname, "pp-mr2", 6);
      filename += 6;
      hashname += 6;
      break;

    default:
      break;
  }

  switch (spec->d_mod_a)
  {
    case REDUCE_D:
      strncpy(filename, "-dma.txt", 8);
      strncpy(hashname, "-dma", 4);
      filename += 8;
      hashname += 4;
      break;

    case NO_REDUCE_D:
      strncpy(filename, ".txt", 4);
      filename += 4;
      break;

    default:
      break;
  }

  strncpy(filename, "\x00", 1);
  strncpy(hashname, "\x00", 1);
}

void test_algorithm(AlgSpec spec, mpz_t* discriminants, int n_discriminants, int discriminant_size, int max_input)
{
  int i, j;

  int n_hashes = max_input * n_discriminants;
  long sum_counters;
  double sum_times;

  char filename[30];
  char hashname[21];

  set_names(filename, hashname, &spec);

  FILE* fptr = fopen(filename, "w");

  sum_counters = 0;
  sum_times = 0.0;
    
  for (i = 0; i < n_discriminants; i++)
    {
      test_discriminant(spec, discriminants[i], max_input, &sum_counters, &sum_times, n_hashes, i * max_input, fptr);
    }

  printf("%20s | %6d | %6d | ", hashname, discriminant_size, spec.coefficient_size);
  printf("%11ld | %13.2e\n", sum_counters/n_hashes, sum_times/n_hashes);

  fclose(fptr);
}
