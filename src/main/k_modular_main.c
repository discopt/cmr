#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include <cmr/matrix.h>
#include <cmr/k_modular.h>

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,       /**< Whether the file format of input/output was defined by the user. */
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

/**
 * \brief Prints the usage of the \p program to stdout.
 * 
 * \returns \c EXIT_FAILURE.
 */

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... FILE\n\n", program);
  puts("Tests matrix in FILE for (strong) k-modularity (and determines k).");
  puts("Options:");
  puts("  -i FORMAT  Format of input FILE; default: `dense'.");
  puts("  -t         Test the transpose matrix instead.");
  puts("  -s         Test for strong k-modularity.");
  puts("  -u         Test only for unimodularity, i.e., 1-modularity.");
  puts("Formats for matrices: dense, sparse");
  puts("If FILE is `-', then the input will be read from stdin.");

  return EXIT_FAILURE;
}

/**
 * \brief Tests matrix from a file for total unimodularity.
 */

static
CMR_ERROR testKmodularity(
  const char* instanceFileName, /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,       /**< Format of the input matrix. */
  bool transpose,               /**< Whether to test the transpose matrix. */
  bool strong,                  /**< Whether to test for strong k-modularity. */
  bool unimodular               /**< Whether to only test for unimodularity. */
)
{
  assert(!transpose || !strong);

  clock_t startClock, endTime;
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  startClock = clock();
  CMR_CHRMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, instanceFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, instanceFile, &matrix) );
  if (instanceFile != stdin)
    fclose(instanceFile);
  fprintf(stderr, "Read %lux%lu matrix with %lu nonzeros.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros);

  /* Actual test. */

  bool propertyOriginal, propertyTranspose;
  size_t kOriginal, kTranspose;
  bool checkOriginal = !transpose;
  bool checkTranspose = transpose || strong;

  if (checkOriginal)
  {
    startClock = clock();
    if (unimodular)
    {
      CMR_CALL( CMRtestUnimodularity(cmr, matrix, &propertyOriginal) );
      fprintf(stderr, "Determined in %f seconds that it is %sunimodular.\n",
        (clock() - startClock) * 1.0 / CLOCKS_PER_SEC, propertyOriginal ? "" : "NOT ");
    }
    else
    {
      CMR_CALL( CMRtestKmodularity(cmr, matrix, &propertyOriginal, &kOriginal) );
      fprintf(stderr, "Determined in %f seconds that it is ", (clock() - startClock) * 1.0 / CLOCKS_PER_SEC);
      if (propertyOriginal)
        fprintf(stderr, "%lu-modular.\n", kOriginal);
      else
        fprintf(stderr, "NOT k-modular.\n");
    }
  }

  if (checkTranspose)
  {
    startClock = clock();
    CMR_CHRMAT* transposed = NULL;
    CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transposed) );

    if (unimodular)
    {
      CMR_CALL( CMRtestUnimodularity(cmr, transposed, &propertyTranspose) );
      fprintf(stderr, "Determined in %f seconds that its transpose is %sunimodular.\n",
        (clock() - startClock) * 1.0 / CLOCKS_PER_SEC, propertyTranspose ? "" : "NOT ");
    }
    else
    {
      CMR_CALL( CMRtestKmodularity(cmr, transposed, &propertyOriginal, &kTranspose) );
      fprintf(stderr, "Determined in %f seconds that its transpose is ", (clock() - startClock) * 1.0 / CLOCKS_PER_SEC);
      if (propertyOriginal)
        fprintf(stderr, "%lu-modular.\n", kTranspose);
      else
        fprintf(stderr, "NOT k-modular.\n");
    }

    CMR_CALL( CMRchrmatFree(cmr, &transposed));
  }

  if (strong)
  {
    fprintf(stderr, "The matrix is %s%smodular.\n", propertyOriginal && propertyTranspose ? "" : "NOT ",
      unimodular ? "uni" : "k-");
  }

  /* Cleanup. */

  CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

int main(int argc, char** argv)
{
  FileFormat inputFormat = FILEFORMAT_UNDEFINED;
  bool transpose = false;
  bool strong = false;
  bool unimodular = false;
  char* instanceFileName = NULL;
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
    puts("Asked for transpose and strong k-modularity at once.\n");
    return printUsage(argv[0]);
  }

  if (inputFormat == FILEFORMAT_UNDEFINED)
    inputFormat = FILEFORMAT_MATRIX_DENSE;

  CMR_ERROR error;
  error = testKmodularity(instanceFileName, inputFormat, transpose, strong, unimodular);

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
