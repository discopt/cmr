#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <float.h>

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/equimodular.h>

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,       /**< Whether the file format of input/output was defined by the user. */
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

/**
 * \brief Tests matrix from a file for total unimodularity.
 */

static
CMR_ERROR testEquimodularity(
  const char* inputMatrixFileName,  /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,           /**< Format of the input matrix. */
  bool transpose,                   /**< Whether to test the transpose matrix. */
  bool strong,                      /**< Whether to test for strong equimodularity. */
  bool unimodular,                  /**< Whether to only test for unimodularity. */
  bool printStats,                  /**< Whether to print statistics to stderr. */
  double timeLimit                  /**< Time limit to impose. */
)
{
  assert(!transpose || !strong);

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_INTMAT* matrix = NULL;
  clock_t readClock = clock();
  CMR_ERROR error = CMR_OKAY;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    error = CMRintmatCreateFromDenseFile(cmr, inputMatrixFileName, "-", &matrix);
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    error = CMRintmatCreateFromSparseFile(cmr, inputMatrixFileName, "-", &matrix);
  else
    CMR_CALL(CMR_ERROR_INVALID);

  if (error)
  {
    fprintf(stderr, "Input error: %s\n", CMRgetErrorMessage(cmr));
    CMR_CALL( CMRfreeEnvironment(&cmr) );
    return CMR_ERROR_INPUT;
  }

  fprintf(stderr, "Read %zux%zu matrix with %zu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Actual test. */

  bool propertyOriginal, propertyTranspose;
  int64_t kOriginal, kTranspose;
  bool checkOriginal = !transpose || strong;
  bool checkTranspose = transpose || strong;
  time_t startClock;

  CMR_EQUIMODULAR_PARAMS params;
  CMR_CALL( CMRequimodularParamsInit(&params) );
  CMR_EQUIMODULAR_STATS stats;
  CMR_CALL( CMRequimodularStatsInit(&stats));

  if (checkOriginal)
  {
    startClock = clock();
    if (unimodular)
    {
      CMR_CALL( CMRunimodularTest(cmr, matrix, &propertyOriginal, &params, &stats, timeLimit) );
      fprintf(stderr, "Determined in %f seconds that it is %sunimodular.\n",
        (clock() - startClock) * 1.0 / CLOCKS_PER_SEC, propertyOriginal ? "" : "NOT ");
    }
    else
    {
      kOriginal = 0;
      CMR_CALL( CMRequimodularTest(cmr, matrix, &propertyOriginal, &kOriginal, &params, &stats, timeLimit) );
      fprintf(stderr, "Determined in %f seconds that it is ", (clock() - startClock) * 1.0 / CLOCKS_PER_SEC);
      if (propertyOriginal)
        fprintf(stderr, "equimodular with determinant gcd %" PRId64 ".\n", kOriginal);
      else
        fprintf(stderr, "NOT equimodular.\n");
    }
  }

  if (checkTranspose)
  {
    startClock = clock();
    CMR_INTMAT* transposed = NULL;
    CMR_CALL( CMRintmatTranspose(cmr, matrix, &transposed) );

    if (unimodular)
    {
      CMR_CALL( CMRunimodularTest(cmr, transposed, &propertyTranspose, &params, &stats, timeLimit - stats.totalTime) );
      fprintf(stderr, "Determined in %f seconds that its transpose is %sunimodular.\n",
        (clock() - startClock) * 1.0 / CLOCKS_PER_SEC, propertyTranspose ? "" : "NOT ");
    }
    else
    {
      kTranspose = 0;
      CMR_CALL( CMRequimodularTest(cmr, transposed, &propertyOriginal, &kTranspose, &params, &stats,
        timeLimit - stats.totalTime) );
      fprintf(stderr, "Determined in %f seconds that its transpose is ", (clock() - startClock) * 1.0 / CLOCKS_PER_SEC);
      if (propertyOriginal)
        fprintf(stderr, "equimodular with determinant gcd %" PRId64 ".\n", kTranspose);
      else
        fprintf(stderr, "NOT equimodular.\n");
    }

    CMR_CALL( CMRintmatFree(cmr, &transposed));
  }

  if (strong)
  {
    if (unimodular)
    {
      fprintf(stderr, "The matrix is %sstrongly unimodular.\n", propertyOriginal && propertyTranspose ? "" : "NOT ");
    }
    else
    {
      if (propertyOriginal && propertyTranspose)
      {
        assert(kOriginal == kTranspose);
        fprintf(stderr, "The matrix is strongly equimodular with k = %" PRId64 ".\n", kOriginal);
      }
      else
        fprintf(stderr, "The matrix is NOT strongly equimodular.\n");
    }
  }

  if (printStats)
    CMR_CALL( CMRequimodularStatsPrint(stderr, &stats, "") );

  /* Cleanup. */

  CMR_CALL( CMRintmatFree(cmr, &matrix) );
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
  fputs("  determines whether the matrix given in file IN-MAT is (strongly) equimodular for determinant gcd k.\n",
    stderr);
  fputs("\n", stderr);

  fputs("Options:\n", stderr);
  fputs("  -i FORMAT  Format of file IN-MAT; default: dense.\n", stderr);
  fputs("  -t         Test the transpose matrix instead.\n", stderr);
  fputs("  -s         Test for strong equimodularity.\n", stderr);
  fputs("  -u         Test only for unimodularity, i.e., k = 1.\n", stderr);
  fputs("\n", stderr);

  fputs("Advanced options:\n", stderr);
  fputs("  --stats              Print statistics about the computation to stderr.\n", stderr);
  fputs("  --time-limit LIMIT   Allow at most LIMIT seconds for the computation.\n", stderr);
  fputs("\n", stderr);

  fputs("Formats for matrices: dense, sparse\n", stderr);
  fputs("If IN-MAT is `-', then the input will be read from stdin.\n", stderr);

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  FileFormat inputFormat = FILEFORMAT_UNDEFINED;
  bool transpose = false;
  bool strong = false;
  bool unimodular = false;
  bool printStats = false;
  char* instanceFileName = NULL;
  double timeLimit = DBL_MAX;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-t"))
      transpose = true;
    else if (!strcmp(argv[a], "-s"))
      strong = true;
    else if (!strcmp(argv[a], "-u"))
      unimodular = true;
    else if (!strcmp(argv[a], "--stats"))
      printStats = true;
    else if (!strcmp(argv[a], "-i") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        inputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        inputFormat = FILEFORMAT_MATRIX_SPARSE;
      else
      {
        printf("Error: unknown input file format <%s>.\n\n", argv[a+1]);
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
    else if (!instanceFileName)
      instanceFileName = argv[a];
    else
    {
      printf("Error: Two input files <%s> and <%s> specified.\n\n", instanceFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (!instanceFileName)
  {
    puts("No input file specified.\n");
    return printUsage(argv[0]);
  }
  else if (transpose && strong)
  {
    puts("Asked for transpose and strong equimodularity at once.\n");
    return printUsage(argv[0]);
  }

  if (inputFormat == FILEFORMAT_UNDEFINED)
    inputFormat = FILEFORMAT_MATRIX_DENSE;

  CMR_ERROR error;
  error = testEquimodularity(instanceFileName, inputFormat, transpose, strong, unimodular, printStats, timeLimit);

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

