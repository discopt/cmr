#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <stdint.h>

#include <cmr/matrix.h>
#include <cmr/named.h>

typedef enum
{
  NONE,
  R_10,
  R_12,
  F_7,
  K_3_3,
  K_5,
  I_k,
} NamedMatrix;

typedef enum
{
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

static
CMR_ERROR recognize(
  const char* inputMatrixFileName,
  FileFormat inputFormat,
  double epsilon)
{
  clock_t readClock = clock();
  FILE* inputMatrixFile = strcmp(inputMatrixFileName, "-") ? fopen(inputMatrixFileName, "r") : stdin;
  if (!inputMatrixFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_DBLMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRdblmatCreateFromDenseStream(cmr, inputMatrixFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRdblmatCreateFromSparseStream(cmr, inputMatrixFile, &matrix) );
  if (inputMatrixFile != stdin)
    fclose(inputMatrixFile);
  fprintf(stderr, "Read %zux%zu matrix with %zu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  CMR_CHRMAT* mat = NULL;
  CMR_ERROR error = CMRdblmatToChr(cmr, matrix, epsilon, &mat);
  if (error != CMR_ERROR_INPUT && error != CMR_ERROR_OVERFLOW)
  {
    CMR_CALL( error );

    size_t isR10;
    CMR_CALL( CMRisR10Matrix(cmr, mat, &isR10) );
    if (isR10)
      printf("R_10.%zu\n", isR10);

    size_t isR12;
    CMR_CALL( CMRisR12Matrix(cmr, mat, &isR12) );
    if (isR12)
      printf("R_12.%zu\n", isR10);

    size_t isIdentity;
    CMR_CALL( CMRisIdentityMatrix(cmr, mat, &isIdentity) );
    if (isIdentity)
      printf("I_%zu\n", isIdentity);

    CMR_CALL( CMRchrmatFree(cmr, &mat) );
  }

  CMR_CALL( CMRdblmatFree(cmr, &matrix) );

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

static
CMR_ERROR generate(
  const char* outputMatrixFileName,
  FileFormat outputFormat,
  NamedMatrix generateNamed,
  bool generateDual,
  size_t generateIndex
)
{
  assert(outputMatrixFileName);
  assert(generateNamed != NONE);

  char name[32];
  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;

  if (generateIndex == 0)
    generateIndex = 1;

  switch(generateNamed)
  {
  case R_10:
    CMR_CALL( CMRcreateR10Matrix(cmr, generateIndex, &matrix) );
    sprintf(name, "R_10");
  break;
  case R_12:
    CMR_CALL( CMRcreateR12Matrix(cmr, generateIndex, &matrix) );
    sprintf(name, "R_12");
  break;
  case K_5:
    CMR_CALL( CMRcreateK5Matrix(cmr, generateIndex, &matrix) );
    sprintf(name, "K_5");
  break;
  case K_3_3:
    CMR_CALL( CMRcreateK33Matrix(cmr, generateIndex, &matrix) );
    sprintf(name, "K_3,3");
  break;
  case I_k:
    CMR_CALL( CMRcreateIdentityMatrix(cmr, generateIndex, &matrix) );
    sprintf(name, "I_%zu", generateIndex);
  break;
  default:
    assert(!"Named matroid generation not implemented.");
    return CMR_ERROR_INPUT;
  }

  if (generateDual)
  {
    CMR_CHRMAT* dual = NULL;
    CMR_CALL( CMRchrmatTranspose(cmr, matrix, &dual) );
    strcat(name, "*");
    CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    matrix = dual;
  }

  /* Write to file. */

  bool outputMatrixToFile = strcmp(outputMatrixFileName, "-");
  FILE* outputMatrixFile = outputMatrixToFile ? fopen(outputMatrixFileName, "w") : stdout;
  fprintf(stderr, "Writing %s matrix to %s%s%s in %s format.\n", name, outputMatrixToFile ? "file <" : "",
    outputMatrixToFile ? outputMatrixFileName : "stdout", outputMatrixToFile ? ">" : "",
    outputFormat == FILEFORMAT_MATRIX_DENSE ? "dense" : "sparse");
  if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatPrintDense(cmr, matrix, outputMatrixFile, '0', false) );
  else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatPrintSparse(cmr, matrix, outputMatrixFile) );
  else
    assert(false);
  if (outputMatrixToFile)
    fclose(outputMatrixFile);

  /* Cleanup. */

  CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

/**
 * \brief Prints the usage of the \p program to stdout.
 *
 * \returns \c EXIT_FAILURE.
 */

int printUsage(const char* program)
{
  fputs("Usage:\n", stderr);

  fprintf(stderr, "%s IN-MAT [OPTION]...\n", program);
  fputs("  (1) determines whether the matrix given in file IN-MAT is (up to row/column permutations)\n", stderr);
  fputs("      one of the named matrices below.\n", stderr);
  fputs("\n", stderr);

  fprintf(stderr, "%s OUT-MAT -g NAME [OPTION]...\n", program);
  fputs("  (2) generates one of the matrices below.\n", stderr);
  fputs("\n", stderr);

  fputs("Representation matrices of:\n", stderr);
  fputs("  R_10       regular matroid R_10.\n", stderr);
  // fputs("  R_12       regular matroid R_12.\n", stderr);
  // fputs("  F_7        irregular Fano-matroid F_7.\n", stderr);
  // fputs("  K_3,3      graphic matroid of complete bipartite graph with 3+3 vertices.\n", stderr);
  // fputs("  K_5        graphic matroid of complete graph with 5 vertices.\n", stderr);
  fputs("Variants (can be combined):\n", stderr);
  fputs("  <NAME>*    refers to the dual matroid:\n", stderr);
  fputs("  <NAME>.<k> refers to a specific representation matrix, k = 1,2,...\n", stderr);
  fputs("\n", stderr);
  fputs("Other matrices:\n", stderr);
  fputs("  I_<SIZE>   Identity matrix of order SIZE.\n", stderr);
  fputs("\n", stderr);

  fputs("Options specific to (1):\n", stderr);
  fputs("  -i         Format of file IN-MAT; default: dense.\n", stderr);
  fputs("  -e EPSILON Allows rounding of numbers up to tolerance EPSILON; default: 1.0e-9.\n", stderr);
  fputs("\n", stderr);

  fputs("Options specific to (2):\n", stderr);
  fputs("  -o         Format of file OUT-MAT; default: dense.\n", stderr);
  fputs("\n", stderr);

  fputs("Formats for matrices: dense, sparse\n", stderr);
  fputs("If IN-MAT is `-' then the input matrix is read from stdin.\n", stderr);
  fputs("If OUT-MAT is `-' then the output matrix is written to stdout.\n", stderr);

  return EXIT_FAILURE;
}

static
bool parseNamed(
  const char* name,
  bool* pgenerateDual,
  char indexSeparator,
  size_t numIndices,
  size_t* pgenerateIndex
)
{
  if (pgenerateDual)
    *pgenerateDual = false;
  *pgenerateIndex = 0;
  if (*name == '\0')
    return true;

  if (pgenerateDual && *name == '*')
  {
    *pgenerateDual = true;
    ++name;
  }

  if (*name == indexSeparator)
  {
    ++name;
    char* remaining;
    int k = strtol(name, &remaining, 10);
    if (remaining == name || k <= 0 || (size_t)k > numIndices)
      return false;

    *pgenerateIndex = k;

    if (pgenerateDual && !*pgenerateDual && (*remaining == '*'))
    {
      *pgenerateDual = true;
      ++remaining;
    }

    return (*remaining == '\0');
  }
  else
    return (*name == '\0');
}

int main(int argc, char** argv)
{
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  FileFormat outputFormat = FILEFORMAT_MATRIX_DENSE;
  double epsilon = 1.0e-9;
  NamedMatrix generateNamed = NONE;
  bool generateDual = false;
  size_t generateIndex = 0;
  char* matrixFileName = NULL;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-g") && a+1 < argc)
    {
      if (!strncmp("R_10", argv[a+1], 4))
      {
        if (!parseNamed(&argv[a+1][4], &generateDual, '.', 2, &generateIndex))
        {
          fprintf(stderr, "Error: Invalid suffix for named R_10 matroid <%s>.\n\n", argv[a+1]);
          return printUsage(argv[0]);
        }
        generateNamed = R_10;
      }
      else if (!strncmp("R_12", argv[a+1], 4))
      {
        if (!parseNamed(&argv[a+1][4], &generateDual, '.', 2, &generateIndex))
        {
          fprintf(stderr, "Error: Invalid suffix for named R_12 matroid <%s>.\n\n", argv[a+1]);
          return printUsage(argv[0]);
        }
        generateNamed = R_12;
      }
      else if (!strncmp("K_3,3", argv[a+1], 5))
      {
        if (!parseNamed(&argv[a+1][5], &generateDual, '.', 2, &generateIndex))
        {
          fprintf(stderr, "Error: Invalid suffix for named K_3,3 matroid <%s>.\n\n", argv[a+1]);
          return printUsage(argv[0]);
        }
        generateNamed = K_3_3;
      }
      else if (!strncmp("K_5", argv[a+1], 3))
      {
        if (!parseNamed(&argv[a+1][3], &generateDual, '.', 1, &generateIndex))
        {
          fprintf(stderr, "Error: Invalid suffix for named K_5 matroid <%s>.\n\n", argv[a+1]);
          return printUsage(argv[0]);
        }
        generateNamed = K_5;
      }
      else if (!strncmp("F_7", argv[a+1], 3))
      {
        if (!parseNamed(&argv[a+1][3], &generateDual, '.', 2, &generateIndex))
        {
          fprintf(stderr, "Error: Invalid suffix for named F_7 matroid <%s>.\n\n", argv[a+1]);
          return printUsage(argv[0]);
        }
        generateNamed = F_7;
      }
      else if (!strncmp("I", argv[a+1], 1))
      {
        if (!parseNamed(&argv[a+1][1], NULL, '_', SIZE_MAX, &generateIndex))
        {
          fprintf(stderr, "Error: Invalid suffix for identity matroid <%s>.\n\n", argv[a+1]);
          return printUsage(argv[0]);
        }
        generateNamed = I_k;
      }
      else
      {
        fprintf(stderr, "Error: Invalid named matroid <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-i") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        inputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        inputFormat = FILEFORMAT_MATRIX_SPARSE;
      else
      {
        fprintf(stderr, "Error: Unknown input file format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-o") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        outputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        outputFormat = FILEFORMAT_MATRIX_SPARSE;
      else
      {
        fprintf(stderr, "Error: Unknown output file format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-e") && a+1 < argc)
    {
      char* p;
      epsilon = strtod(argv[a+1], &p);
      if (*p != '\0' || epsilon < 0.0 || epsilon > 0.5)
      {
        fprintf(stderr, "Error: Invalid tolerance <%s>", argv[a+1]);
        return printUsage(argv[0]);
      }
      a++;
    }
    else if (!matrixFileName)
      matrixFileName = argv[a];
    else
    {
      printf("Error: Two input/output matrix files <%s> and <%s> specified.\n\n", matrixFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (!matrixFileName)
  {
    fputs("Error: No input/output matrix specified.\n\n", stderr);
    return printUsage(argv[0]);
  }

  CMR_ERROR error;
  if (generateNamed == NONE)
  {
    error = recognize(matrixFileName, inputFormat, epsilon);
  }
  else
  {
    error = generate(matrixFileName, outputFormat, generateNamed, generateDual, generateIndex);
  }

  switch (error)
  {
  case CMR_ERROR_INPUT:
    puts("Input error.");
    return EXIT_FAILURE;
  case CMR_ERROR_MEMORY:
    puts("Memory error.");
    return EXIT_FAILURE;
  default:
    return EXIT_SUCCESS;
  }
}

