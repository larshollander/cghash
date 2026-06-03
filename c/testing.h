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

#define SMALL_PRIME_SIZES 255, 383, 511, 767, 1023, 1535, 2047, 3071
#define N_SPS 8

#define SMALL_PRIME_COUNTS 1, 2, 3, 4, 5, 6, 7, 8
#define N_SPC 8

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
  switch (spec->factor_type * spec->d_mod_a)
  {
    case PRIME * NO_REDUCE_D:
      strncpy(filename, "logs/xxx_pr.txt", 15);
      strncpy(hashname, "xxx_pr", 6);
      break;

    case PRIME * REDUCE_D:
      strncpy(filename, "logs/xxx_pr_dma.txt", 19);
      strncpy(hashname, "xxx_pr_dma", 10);
      break;

    case PSEUDOPRIME * NO_REDUCE_D:
      strncpy(filename, "logs/xxx_pp.txt", 15);
      strncpy(hashname, "xxx_pp", 6);
      break;

    case PSEUDOPRIME * REDUCE_D:
      strncpy(filename, "logs/xxx_pp_dma.txt", 19);
      strncpy(hashname, "xxx_pp", 10);
      break;

    case PSEUDOPRIME_MR2 * NO_REDUCE_D:
      strncpy(filename, "logs/xxx_pp_mr2.txt", 19);
      strncpy(hashname, "xxx_pp_mr2", 10);
      break;

    case PSEUDOPRIME_MR2 * REDUCE_D:
      strncpy(filename, "logs/xxx_pp_dma_mr2.txt", 23);
      strncpy(hashname, "xxx_pp_dma_mr2", 14);
      break;

    default:
      break;
  }

  switch (spec->algorithm)
  {
    case WESOLOWSKI:
      strncpy(filename, "logs/wes", 8);
      strncpy(hashname, "wes", 3);
      break;

    case CLR:
      strncpy(filename, "logs/clr", 8);
      strncpy(hashname, "clr", 3);
      break;

    default:
      break;
  }
}

void test_algorithm(AlgSpec spec, mpz_t* discriminants, int n_discriminants, int discriminant_size, int max_input)
{
  int i, j;

  int n_hashes = max_input * n_discriminants;
  int small_prime_sizes[N_SPS] = { SMALL_PRIME_SIZES };
  int small_prime_counts[N_SPC] = { SMALL_PRIME_COUNTS };

  long sum_counters;
  double sum_times;

  char filename[24];
  char hashname[15];

  for (i = 0; i < 20; i++) { filename[i] = '\x00'; }
  for (i = 0; i < 11; i++) { hashname[i] = '\x00'; }

  set_names(filename, hashname, &spec);

  FILE* fptr = fopen(filename, "w");

  if (spec.algorithm == WESOLOWSKI)
  {
    for (j = 0; j < N_SPS; j++)
    {
      int small_prime_size = small_prime_sizes[j];

      if (small_prime_size > discriminant_size/2 - 1) { break; }

      fprintf(fptr, "\nk = %d\n", small_prime_size);
      printf("\n  (%s-%d-%d)\n", hashname, discriminant_size, small_prime_size);

      spec.coefficient_size = small_prime_size;

      sum_counters = 0;
      sum_times = 0.0;

      for (i = 0; i < n_discriminants; i++)
      {
        test_discriminant(spec, discriminants[i], max_input, &sum_counters, &sum_times, n_hashes, i * max_input, fptr);
      }

      printf("\navg. steps: %ld\navg. time (s): %.2e\n", sum_counters/n_hashes, sum_times/n_hashes);
    }
  }

  if (spec.algorithm == CLR)
  {
    for (j = 0; j < N_SPC; j++)
    {
      int n_small_primes = small_prime_counts[j];

      fprintf(fptr, "\nk = %d\n", n_small_primes);
      printf("\n  (%s-%d-%d)\n", hashname, n_small_primes, discriminant_size);

      spec.n_small_primes = n_small_primes;

      sum_counters = 0;
      sum_times = 0.0;

      for (i = 0; i < n_discriminants; i++)
      {
        test_discriminant(spec, discriminants[i], max_input, &sum_counters, &sum_times, n_hashes, i * max_input, fptr);
      }

      printf("\navg. steps: %ld\navg. time (s): %.2e\n", sum_counters/n_hashes, sum_times/n_hashes);
    }
  }

  fclose(fptr);
}
