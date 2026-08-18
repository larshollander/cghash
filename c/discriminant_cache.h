#ifndef CTYPE_INCLUDED
#define CTYPE_INCLUDED
#include <ctype.h>
#endif

#ifndef GMP_INCLUDED
#define GMP_INCLUDED
#include <gmp.h>
#endif

#ifndef STDIO_INCLUDED
#define STDIO_INCLUDED
#include <stdio.h>
#endif

#ifndef STDLIB_INCLUDED
#define STDLIB_INCLUDED
#include <stdlib.h>
#endif

#ifndef STRING_INCLUDED
#define STRING_INCLUDED
#include <string.h>
#endif

// #define SHOW_PROGRESS

/**
 * Saves an array of mpz_t values to a text file with a given key.
 * Overwrites any existing entry with the same key.
 */
int save_mpz_array_text(const char *filename, int key, const mpz_t *arr, size_t count) {
  FILE *fp = fopen(filename, "r");
  FILE *tmp_fp = fopen("gmp_tmp.txt", "w");
  
  if (!tmp_fp) {
    if (fp) fclose(fp);
    return -1; // Cannot create temp file
  }

  int found_key = 0;

  // If file exists, copy lines that DO NOT match the key we are updating
  if (fp) {
    char line[20480]; // Buffer for reading lines
    while (fgets(line, sizeof(line), fp)) {
      int file_key;
      // Peek at the first number to see if it's our key
      if (sscanf(line, "%u", &file_key) == 1) {
        if (file_key == key) {
          found_key = 1;
          // Skip this line (we will write the new one later)
          continue; 
        }
      }
      // Copy the line to temp file
      fputs(line, tmp_fp);
    }
    fclose(fp);
  }

  // Write the new/updated entry
  fprintf(tmp_fp, "%u %zu", key, count);
  for (size_t i = 0; i < count; ++i) {
    // gmp_printf supports %Zd for decimal output
    fprintf(tmp_fp, " ");
    mpz_out_str(tmp_fp, 10, arr[i]);
  }
  fprintf(tmp_fp, "\n");

  fclose(tmp_fp);

  // Replace original file with temp file
  remove(filename);
  rename("gmp_tmp.txt", filename);

  return 0;
}

/**
 * Loads an array of mpz_t values from a text file by key.
 * Returns the number of elements loaded.
 */
size_t load_mpz_array_text(const char *filename, int key, mpz_t *arr, size_t arr_size) {
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    return 0;
  }

  char line[20480];

  while (fgets(line, sizeof(line), fp)) {
    int file_key;
    size_t file_count;

    // Parse the header: KEY COUNT
    // We use sscanf to peek at the start of the line
    if (sscanf(line, "%u %zu", &file_key, &file_count) != 2) {
      continue; // Malformed line, skip
    }

    if (file_key == key) {
      // Found the target!

      // Read at most arr_size values
      file_count = arr_size < file_count ? arr_size : file_count;

      // Initialize arrary entries
      for (size_t i = 0; i < file_count; i++) {
        mpz_init(arr[i]);
      }

      // Parse the rest of the line
      // We need to skip the first two numbers we already read
      char *ptr = line;
      // Advance past key and count
      while (*ptr && !isspace(*ptr)) ptr++; // skip key
      while (*ptr && isspace(*ptr)) ptr++;  // skip space
      while (*ptr && !isspace(*ptr)) ptr++; // skip count
      while (*ptr && isspace(*ptr)) ptr++;  // skip space

      for (size_t i = 0; i < file_count; ++i) {
        // Find end of current token
        char *start = ptr;
        while (*ptr && !isspace(*ptr)) ptr++;
        *ptr = '\0'; // Null terminate the token
        
        if (strlen(start) > 0) {
          mpz_set_str(arr[i], start, 10);
        } else {
          mpz_set_ui(arr[i], 0);
        }
        ptr++; // Move past the null terminator for next iteration
      }

      fclose(fp);
      return file_count;
    }
  }

  fclose(fp);
  return 0; // Key not found
}

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

void sample_discriminants(mpz_t* discriminants, int n_read, int n_discriminants, int discriminant_size)
{
  gmp_randstate_t state;
  gmp_randinit_default(state);

#ifdef SHOW_PROGRESS
  printf("\rdiscriminants sampled: %d/%d\n\033[A", n_read, n_discriminants);
#endif

  for (int i = n_read; i < n_discriminants; i++)
  {
    mpz_init(discriminants[i]);
    sample_prime_discriminant(discriminants[i], state, discriminant_size);
#ifdef SHOW_PROGRESS
    printf("\rdiscriminants sampled: %d/%d\n\033[A", i + 1, n_discriminants);
#endif
  }

  printf("\n");
}

