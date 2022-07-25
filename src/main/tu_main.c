#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <cmr/matrix.h>
#include <cmr/tu.h>

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
CMR_ERROR testTotalUnimodularity(
  const char* instanceFileName, /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,       /**< Format of the input matrix. */
  FileFormat outputFormat,      /**< Format of the output submatrix. */
  bool directGraphicness,       /**< Whether to use fast graphicness routines. */
  bool seriesParallel,          /**< Whether to allow series-parallel operations in the decomposition tree. */
  bool printTree,               /**< Whether to print the decomposition tree. */
  bool outputNonTUElements,     /**< Whether to print the elements of a non-TU submatrix. */
  bool outputNonTUMatrix,       /**< Whether to print a non-TU submatrix. */
  bool printStats               /**< Whether to print statistics to stderr. */
)
{
  clock_t readClock = clock();
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
  {
    fprintf(stderr, "Unable to open file <%s>\n", instanceFileName);
    return CMR_ERROR_INPUT;
  }

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_CHRMAT* matrix = NULL;
  CMR_ERROR error;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
  {
    error = CMRchrmatCreateFromDenseStream(cmr, instanceFile, &matrix);
    if (error == CMR_ERROR_INPUT)
      fprintf(stderr, "Error when reading dense matrix from <%s>: %s\n", instanceFileName, CMRgetErrorMessage(cmr));
  }
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
  {
    error = CMRchrmatCreateFromSparseStream(cmr, instanceFile, &matrix);
    if (error == CMR_ERROR_INPUT)
      fprintf(stderr, "Error when reading dense matrix from <%s>: %s\n", instanceFileName, CMRgetErrorMessage(cmr));
  }
  if (instanceFile != stdin)
    fclose(instanceFile);
  CMR_CALL(error);
  fprintf(stderr, "Read %lux%lu matrix with %lu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Actual test. */

  bool isTU;
  CMR_DEC* decomposition = NULL;
  CMR_SUBMAT* submatrix = NULL;
  CMR_TU_PARAMETERS params;
  CMR_CALL( CMRparamsTotalUnimodularityInit(&params) );
  params.regular.directGraphicness = directGraphicness;
  params.regular.seriesParallel = seriesParallel;
  CMR_TU_STATISTICS stats;
  CMR_CALL(CMRstatsTotalUnimodularityInit(&stats));
  CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, printTree ? &decomposition : NULL,
    (outputNonTUMatrix || outputNonTUElements) ? &submatrix : NULL, &params, &stats) );

  fprintf(stderr, "Matrix %stotally unimodular.\n", isTU ? "IS " : "IS NOT ");
  if (printStats)
    CMR_CALL( CMRstatsTotalUnimodularityPrint(stderr, &stats, NULL) );

  if (decomposition)
  {
    CMR_CALL( CMRdecFree(cmr, &decomposition) );
  }

  if (submatrix)
  {
    if (outputNonTUElements)
    {
      fprintf(stderr, "\nNon-TU submatrix consists of these elements:\n");
      printf("%ld rows:", submatrix->numRows);
      for (size_t r = 0; r < submatrix->numRows; ++r)
        printf(" %ld", submatrix->rows[r]+1);
      printf("\n%ld columns: ", submatrix->numColumns);
      for (size_t c = 0; c < submatrix->numColumns; ++c)
        printf(" %ld", submatrix->columns[c]+1);
      printf("\n");
    }

    if (outputNonTUMatrix)
    {
      CMR_CHRMAT* violatorMatrix = NULL;
      CMR_CALL( CMRchrmatZoomSubmat(cmr, matrix, submatrix, &violatorMatrix) );
      fprintf(stderr, "\nExtracted %lux%lu non-TU submatrix with %lu nonzeros.\n", violatorMatrix->numRows,
        violatorMatrix->numColumns, violatorMatrix->numNonzeros);
      if (outputFormat == FILEFORMAT_MATRIX_DENSE)
        CMR_CALL( CMRchrmatPrintDense(cmr, violatorMatrix, stdout, '0', false) );
      else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
        CMR_CALL( CMRchrmatPrintSparse(cmr, violatorMatrix, stdout) );
      CMR_CALL( CMRchrmatFree(cmr, &violatorMatrix) );
    }

    CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
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
  puts("  -d         Output the decomposition tree of the underlying regular matroid.");
  puts("  -n         Output the elements of a minimal non-totally-unimodular submatrix.");
  puts("  -N         Output a minimal non-totally-unimodular submatrix.");
  puts("  -s         Print statistics about the computation to stderr.");
  puts("Parameter options:");
  puts("  --no-direct-graphic  Check only 3-connected matrices for regularity.");
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
  bool printTree = false;
  bool outputNonTUElements = false;
  bool outputNonTUMatrix = false;
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
      printTree = true;
    else if (!strcmp(argv[a], "-n"))
      outputNonTUElements = true;
    else if (!strcmp(argv[a], "-N"))
      outputNonTUMatrix = true;
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
  error = testTotalUnimodularity(instanceFileName, inputFormat, outputFormat, directGraphicness, seriesParallel,
    printTree, outputNonTUElements, outputNonTUMatrix, printStats);

  switch (error)
  {
  case CMR_ERROR_INPUT:
    return EXIT_FAILURE;
  case CMR_ERROR_MEMORY:
    puts("Memory error.");
    return EXIT_FAILURE;
  default:
    return EXIT_SUCCESS;
  }
}
