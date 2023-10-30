#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include <cmr/matrix.h>
#include <cmr/tu.h>

typedef enum
{
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

typedef enum
{
  ALGORITHM_DEFAULT = 0,          /**< Default algorithm. */
  ALGORITHM_ENUMERATE = 1,        /**< Enumeration algorithm. */
  ALGORITHM_GRAPH = 2,            /**< Graph-based algorithm. */
} Algorithm;

/**
 * \brief Tests matrix from a file for total unimodularity.
 */

static
CMR_ERROR testBalanced(
  const char* inputMatrixFileName,      /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,               /**< Format of the input matrix. */
  const char* outputSubmatrixFileName,  /**< File name of output file for non-balanced submatrix. */
  bool printStats,                      /**< Whether to print statistics to stderr. */
  Algorithm algorithm,                  /**< Algorithm to use. */
  double timeLimit                      /**< Time limit to impose. */
)
{
  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_CHRMAT* matrix = NULL;
  clock_t readClock = clock();
  CMR_ERROR error = CMR_OKAY;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    error = CMRchrmatCreateFromDenseFile(cmr, inputMatrixFileName, "-", &matrix);
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    error = CMRchrmatCreateFromSparseFile(cmr, inputMatrixFileName, "-", &matrix);
  else
    CMR_CALL(CMR_ERROR_INVALID);

  if (error)
  {
    fprintf(stderr, "Input error: %s\n", CMRgetErrorMessage(cmr));
    CMR_CALL( CMRfreeEnvironment(&cmr) );
    return CMR_ERROR_INPUT;
  }

  fprintf(stderr, "Read %lux%lu matrix with %lu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Actual test. */

  bool isBalanced = false; 
  
  
  printf("Matrix %sbalanced.\n", isBalanced ? "IS " : "IS NOT ");

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
  fprintf(stderr, "%s IN-MAT [OPTION]...\n\n", program);
  fputs("  determines whether the matrix given in file IN-MAT is balanced.\n\n", stderr);
  fputs("Options:\n", stderr);
  fputs("  -i FORMAT  Format of file IN-MAT, among `dense' and `sparse'; default: dense.\n", stderr);
  fputs("  -N NON-SUB Write a minimal non-balanced submatrix to file NON-SUB; default: skip computation.\n", stderr);
  fputs("  -s         Print statistics about the computation to stderr.\n\n", stderr);
  fputs("Advanced options:\n", stderr);
  fputs("  --time-limit LIMIT Allow at most LIMIT seconds for the computation.\n", stderr);
  fputs("  --algorithm ALGO   Algorithm to use, among `enumerate` and `graph`; default: choose best.\n\n", stderr);
  fputs("If IN-MAT is `-' then the matrix is read from stdin.\n", stderr);
  fputs("If NON-SUB is `-' then the submatrix is written to stdout.\n", stderr);

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  char* inputMatrixFileName = NULL;
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  char* outputSubmatrix = NULL;
  bool printStats = false;
  Algorithm algorithm = ALGORITHM_DEFAULT;
  double timeLimit = DBL_MAX;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-i") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        inputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        inputFormat = FILEFORMAT_MATRIX_SPARSE;
      else
      {
        fprintf(stderr, "Error: unknown input file format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-N") && a+1 < argc)
      outputSubmatrix = argv[++a];
    else if (!strcmp(argv[a], "-s"))
      printStats = true;
    else if (!strcmp(argv[a], "--algorithm") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "enumerate"))
        algorithm = ALGORITHM_ENUMERATE;
      else if (!strcmp(argv[a+1], "graph"))
        algorithm = ALGORITHM_GRAPH;
      else
      {
        fprintf(stderr, "Error: Invalid algorithm <%s> specified.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "--time-limit") && (a+1 < argc))
    {
      if (sscanf(argv[a+1], "%lf", &timeLimit) == 0 || timeLimit <= 0)
      {
        fprintf(stderr, "Error: Invalid time limit <%s> specified.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!inputMatrixFileName)
      inputMatrixFileName = argv[a];
    else
    {
      fprintf(stderr, "Error: Two input files <%s> and <%s> specified.\n\n", inputMatrixFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (!inputMatrixFileName)
  {
    fputs("Error: No input file specified.\n\n", stderr);
    return printUsage(argv[0]);
  }

  CMR_ERROR error;
  error = testBalanced(inputMatrixFileName, inputFormat, outputSubmatrix, printStats, algorithm, timeLimit);

  switch (error)
  {
  case CMR_ERROR_INPUT:
    /* The actual function will have reported the details. */
    return EXIT_FAILURE;
  case CMR_ERROR_MEMORY:
    puts("Memory error.");
    return EXIT_FAILURE;
  default:
    return EXIT_SUCCESS;
  }
}
