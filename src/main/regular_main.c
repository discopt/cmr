#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <cmr/matrix.h>
#include <cmr/regular.h>

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
CMR_ERROR testRegularity(
  const char* instanceFileName, /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,       /**< Format of the input matrix. */
  FileFormat outputFormat,      /**< Format of the output submatrix. */
  bool directGraphicness,       /**< Whether to use fast graphicness routines. */
  bool seriesParallel,          /**< Whether to allow series-parallel operations in the decomposition tree. */
  bool completeTree,            /**< Whether to compute the complete decomposition tree. */
  bool planarityCheck,          /**< Whether to test for planarity. */
  bool printTree,               /**< Whether to print the decomposition tree. */
  bool outputSubmatrixElements, /**< Whether to print the elements of a minimal submatrix containing an irregular minor. */
  bool outputSubmatrix,         /**< Whether to print a minimal submatrix containing an irregular minor. */
  bool printStats               /**< Whether to print statistics to stderr. */
)
{
  clock_t readClock = clock();
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_CHRMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, instanceFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, instanceFile, &matrix) );
  if (instanceFile != stdin)
    fclose(instanceFile);
  fprintf(stderr, "Read %lux%lu matrix with %lu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Actual test. */

  bool isRegular;
  CMR_DEC* decomposition = NULL;
  CMR_MINOR* minor = NULL;
  CMR_REGULAR_PARAMETERS params;
  CMR_CALL( CMRparamsRegularInit(&params) );
  params.directGraphicness = directGraphicness;
  params.seriesParallel = seriesParallel;
  params.completeTree = completeTree;
  params.planarityCheck = planarityCheck;
  params.matrices = printTree ? CMR_DEC_CONSTRUCT_ALL : CMR_DEC_CONSTRUCT_NONE;
  CMR_REGULAR_STATISTICS stats;
  CMR_CALL( CMRstatsRegularInit(&stats) );
  CMR_CALL( CMRtestBinaryRegular(cmr, matrix, &isRegular, printTree ? &decomposition : NULL,
    (outputSubmatrix || outputSubmatrixElements) ? &minor : NULL, &params, &stats) );

  fprintf(stderr, "Matrix %sregular.\n", isRegular ? "IS " : "IS NOT ");
  if (printStats)
    CMR_CALL( CMRstatsRegularPrint(stderr, &stats, NULL) );

  if (decomposition)
  {
    if (printTree)
      CMR_CALL( CMRdecPrint(cmr, decomposition, stdout, 2, true, true, true) );
    CMR_CALL( CMRdecFree(cmr, &decomposition) );
  }

  if (minor)
  {
    CMR_CALL( CMRminorFree(cmr, &minor) );
  }

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
  printf("Usage: %s [OPTION]... FILE\n\n", program);
  puts("Tests matrix in FILE for total unimodularity.");
  puts("Options:");
  puts("  -i FORMAT  Format of input FILE; default: `dense'.");
  puts("  -o FORMAT  Format of output matrices; default: `dense'.");
  puts("  -p         Test graphic matrices also for cographicness (i.e., planarity).");
  puts("  -d         Compute the complete decomposition tree of the matroid.");
  puts("  -D         Output the computed decomposition tree of the matroid.");
  puts("  -n         Output the elements of a minimal non-regular submatrix.");
  puts("  -N         Output a minimal non-regular submatrix.");
  puts("  -s         Print statistics about the computation to stderr.");
  puts("Parameter options:");
  puts("  --no-direct-graphic   Check only 3-connected matrices for regularity.");
  puts("  --no-series-parallel  Do not allow series-parallel operations in decomposition tree.");
  puts("Formats for matrices: dense, sparse");
  puts("If FILE is `-', then the input will be read from stdin.");

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  FileFormat inputFormat = FILEFORMAT_UNDEFINED;
  FileFormat outputFormat = FILEFORMAT_MATRIX_DENSE;
  bool directGraphicness = true;
  bool seriesParallel = true;
  bool completeTree = false;
  bool planarityCheck = false;
  bool printTree = false;
  bool outputSubmatrixElements = false;
  bool outputSubmatrix = false;
  bool printStats = false;
  char* instanceFileName = NULL;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "--no-direct-graphic"))
      directGraphicness = false;
    else if (!strcmp(argv[a], "--no-series-parallel"))
      seriesParallel = false;
    else if (!strcmp(argv[a], "-d"))
      completeTree = true;
    else if (!strcmp(argv[a], "-D"))
      printTree = true;
    else if (!strcmp(argv[a], "-p"))
      planarityCheck = true;
    else if (!strcmp(argv[a], "-n"))
      outputSubmatrixElements = true;
    else if (!strcmp(argv[a], "-N"))
      outputSubmatrix = true;
    else if (!strcmp(argv[a], "-s"))
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
    else if (!strcmp(argv[a], "-o") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        outputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        outputFormat = FILEFORMAT_MATRIX_SPARSE;
      else
      {
        printf("Error: unknown output format <%s>.\n\n", argv[a+1]);
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

  if (inputFormat == FILEFORMAT_UNDEFINED)
    inputFormat = FILEFORMAT_MATRIX_DENSE;

  CMR_ERROR error;
  error = testRegularity(instanceFileName, inputFormat, outputFormat, directGraphicness, seriesParallel, completeTree,
    planarityCheck, printTree, outputSubmatrixElements, outputSubmatrix, printStats);

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
