#ifndef STDIO_INCLUDED
#define STDIO_INCLUDED
#include <stdio.h>
#endif

#ifndef STRING_INCLUDED
#define STRING_INCLUDED
#include <string.h>
#endif

#include "testing.h"

#define D_LENGTH 4096
#define D_BASE 10

#define MAX_ARGS_SKIPPED 4

typedef struct
{
  int n_specs;
  int n_discriminants;
  int discriminant_size;
  int max_input;
  int prompt_discriminant;
  int prompt_max_input;
  int prompt_algorithm;
} Args;

void debug_spec(AlgSpec* spec)
{
  printf("algorithm: %d\n", spec->algorithm);
  printf("factor type: %d\n", spec->factor_type);
  printf("use D%%=a optimisation: %d\n", spec->d_mod_a);
  printf("# of small factors: %d\n", spec->n_small_factors);
  printf("coefficient size: %d\n", spec->coefficient_size);
}

void prompt_parameters(mpz_t discriminant, AlgSpec* spec, Args* args)
{
  if (args->prompt_discriminant)
  {
    unsigned char buf[D_LENGTH];
    for (int i = 0; i < D_LENGTH; i++) { buf[i] = 0; }

    printf("discriminant:");
    scanf("%s", buf);

    mpz_init_set_str(discriminant, buf, D_BASE);
    args->discriminant_size = (int)mpz_sizeinbase(discriminant, 2);
  }

  if (args->prompt_max_input)
  {
    printf("number of inputs: ");
    scanf("%d", &args->max_input);
  }

  if (args->prompt_algorithm)
  {
    unsigned char buf[1];
    printf("algorithm ([w]esolowski, [c]lr): ");
    scanf("%s", buf);

    switch (buf[0])
    {
      case 'w':
        spec->algorithm = WESOLOWSKI;
        break;

      case 'c':
        spec->algorithm = CLR_SF;
        break;

      default:
        break;
    }

    printf("factor type (p[r]ime, p[s]eudoprime): ");
    scanf("%s", buf);

    switch (buf[0])
    {
      case 'r':
        spec->factor_type = PRIME;
        break;

      case 's':
        spec->factor_type = PSEUDOPRIME;
        break;

      default:
        break;
    }

    if (spec->factor_type == PSEUDOPRIME)
    {
      printf("use base 2 miller-rabin filter ([y]es, [n]o): ");
      scanf("%s", buf);

      switch (buf[0])
      {
        case 'y':
          spec->factor_type = PSEUDOPRIME_MR2;
          break;

        default:
          break;
      }
    }

    if (spec->algorithm == CLR_SF)
    {
      unsigned char buf[1];
      printf("use large factor ([y]es, [n]o): ");
      scanf("%s", buf);

      if (buf[0] == 'y') { spec->algorithm = CLR_LF; }
    }

    printf("apply D %%= a optimisation ([y]es, [n]o): ");
    scanf("%s", buf);

    switch (buf[0])
    {
      case 'y':
        spec->d_mod_a = REDUCE_D;
        break;

      case 'n':
        spec->d_mod_a = NO_REDUCE_D;
        break;

      default:
        break;
    }
  }
}

Args parse_args(int argc, char** argv)
{
  Args args = { 0, 0, 0, 0, 0, 0, 0 };

  args.n_specs = 1;

  switch (argc)
  {
    case 1:
      args.n_discriminants = 1;
      args.prompt_discriminant = 1;
      args.prompt_max_input = 1;
      args.prompt_algorithm = 1;
      break;

    case 2:
      args.n_discriminants = 1;
      args.discriminant_size = atoi(argv[1]);
      args.prompt_max_input = 1;
      args.prompt_algorithm = 1;
      break;

    case 3:
      args.n_discriminants = atoi(argv[1]);
      args.discriminant_size = atoi(argv[2]);
      args.prompt_max_input = 1;
      args.prompt_algorithm = 1;
      break;

    case 4:
      args.max_input = atoi(argv[1]);
      args.n_discriminants = atoi(argv[2]);
      args.discriminant_size = atoi(argv[3]);
      args.prompt_algorithm = 1;
      break;

    default:
      args.max_input = atoi(argv[1]);
      args.n_discriminants = atoi(argv[2]);
      args.discriminant_size = atoi(argv[3]);
      args.n_specs = argc - 4;
      break;
  }

  return args;
}

AlgSpec parse_spec(char* buf)
{
  int i = 0;
  AlgSpec spec = { 0, 1, 1, 1, 0 };

  // algorithm (always): 
  //   'w' for wesolowski, 
  //   'c' for clr,
  //   else returns invalid spec
  switch (buf[i++])
  {
    case 'w':
      spec.algorithm = WESOLOWSKI;
      break;

    case 'c':
      spec.algorithm = CLR_SF;
      break;

    default:
      return spec;
  }

  // factor type (always): 
  //   'r' for prime, 
  //   's' for pseudoprime,
  //   else returns invalid spec
  switch (buf[i++])
  {
    case 'r':
      spec.factor_type = PRIME;
      break;

    case 's':
      spec.factor_type = PSEUDOPRIME;
      break;

    default:
      return spec;
  }

  // miller-rabin base 2 (optional): 
  //   'm', 
  //   requires pseudoprime factors
  if (buf[i] == 'm')
  {
    spec.factor_type = PSEUDOPRIME_MR2;
    i++;
  }

  // large factor (optional): 
  //   'l', 
  //   requires clr hash
  if (buf[i] == 'l')
  {
    spec.algorithm = CLR_LF;
    i++;
  }

  // D %= a (optional):
  //   'd',
  //   for any spec
  if (buf[i] == 'd') 
  { 
    spec.d_mod_a = REDUCE_D; 
    i++;
  }

  // number of small factors (needed for clr):
  //   decimal integer
  spec.n_small_factors = atoi(buf + i);

  return spec;
}

void parse_specs(AlgSpec* specs, char** bufs, int n_specs)
{
  int skipped = 0;

  for (int i = 0; i < n_specs && skipped <= MAX_ARGS_SKIPPED; )
  {
    specs[i] = parse_spec(bufs[skipped + i]);

    if (specs[i].algorithm != 0) { i++; }
    else { skipped++; }
  }
}

int prompt_continue()
{
  char buf[1];
  printf("\n[q]uit or [c]ontinue: ");
  scanf("%s", buf);
  
  switch (buf[0])
  {
    case 'q':
      return 0;

    default:
      return 1;
  }
}

int main(int argc, char** argv)
{
  int i;

  Args args = parse_args(argc, argv);

  if (args.n_discriminants == 0) { return 1; }

  mpz_t discriminants[args.n_discriminants];
  AlgSpec specs[args.n_specs];

  if (!args.prompt_discriminant) 
  { 
    sample_discriminants(discriminants, args.n_discriminants, args.discriminant_size); 
  }

  if (!args.prompt_algorithm)
  {
    parse_specs(specs, argv, args.n_specs);
  }

  do
  {
    prompt_parameters(discriminants[0], &specs[0], &args);

    for (i = 0; i < args.n_specs; i++)
    {
      test_algorithm(specs[i], discriminants, args.n_discriminants, args.discriminant_size, args.max_input);
    }
  } while (args.prompt_algorithm && prompt_continue());

  for (int i = 0; i < args.n_discriminants; i++)
  {
    mpz_clear(discriminants[i]);
  }

  return 0;
}
