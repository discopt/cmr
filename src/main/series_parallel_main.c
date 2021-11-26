#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <cmr/matrix.h>
#include <cmr/series_parallel.h>

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,
  FILEFORMAT_MATRIX_DENSE = 1,
  FILEFORMAT_MATRIX_SPARSE = 2
} FileFormat;

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... FILE\n\n", program);
  puts("Applies all possible series-parallel reductions to the ternary or binary matrix in FILE.");
  puts("Options:");
  puts("  -i FORMAT  Format of input FILE; default: `dense'.");
  puts("  -o FORMAT  Format of output matrices; default: `dense'.");
  puts("  -sp        Output the list of series-parallel reductions.");
  puts("  -r         Output the elements of the reduced matrix.");
  puts("  -R         Output the reduced matrix.");
  puts("  -n         Output the elements of a minimal non-series-parallel submatrix.");
  puts("  -N         Output a minimal non-series-parallel submatrix.");
  puts("  -s         Print statistics about the computation to stderr.");
  puts("Matrix formats: dense, sparse");
  puts("If FILE is `-', then the input will be read from stdin.");
  return EXIT_FAILURE;
}

CMR_ERROR matrixSeriesParallel2Sums(
  const char* instanceFileName, /**< File name of instance. */
  FileFormat inputFormat,       /**< Format of input matrix. */
  FileFormat outputFormat,      /**< Format of output matrices. */
  bool outputReductions,        /**< Whether to output the list of series-parallel reductions. */
  bool outputReducedElements,   /**< Whether to output the elements of the reduced matrix. */
  bool outputReducedMatrix,     /**< Whether to output the reduced matrix. */
  bool outputNonSPElements,     /**< Whether to output the elements of a non-SP submatrix if not series-parallel. */
  bool outputNonSPMatrix,       /**< Whether to output a non-SP submatrix if not series-parallel. */
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

  /* Run the search. */

  CMR_SP_REDUCTION* reductions = NULL;
  size_t numReductions = 0;
  CMR_CALL( CMRallocBlockArray(cmr, &reductions, matrix->numRows + matrix->numColumns) );
  CMR_SUBMAT* reducedSubmatrix = NULL;
  CMR_SUBMAT* violatorSubmatrix = NULL;

  CMR_SP_STATISTICS stats;
  CMR_CALL( CMRstatsSeriesParallelInit(&stats) );
  CMR_CALL( CMRtestTernarySeriesParallel(cmr, matrix, NULL, reductions, &numReductions,
    (outputReducedElements || outputReducedMatrix) ? &reducedSubmatrix : NULL,
    (outputNonSPElements || outputNonSPMatrix) ? &violatorSubmatrix : NULL, &stats) );

  fprintf(stderr, "Matrix %sseries-parallel. %ld reductions can be applied.\n",
    numReductions == matrix->numRows + matrix->numColumns ? "IS " : "is NOT ", numReductions);
  if (printStats)
    CMR_CALL( CMRstatsSeriesParallelPrint(stderr, &stats, NULL) );

  if (outputReductions)
  {
    fprintf(stderr, "Printing %ld series-parallel reductions.\n", numReductions);
    printf("%ld\n", numReductions);
    for (size_t i = 0; i < numReductions; ++i)
      printf("%s\n", CMRspReductionString(reductions[i], NULL));
  }

  if (outputReducedElements)
  {
    fprintf(stderr, "\nReduced submatrix consists of these elements:\n");
    printf("%ld rows:", reducedSubmatrix->numRows);
    for (size_t r = 0; r < reducedSubmatrix->numRows; ++r)
      printf(" %ld", reducedSubmatrix->rows[r]+1);
    printf("\n%ld columns: ", reducedSubmatrix->numColumns);
    for (size_t c = 0; c < reducedSubmatrix->numColumns; ++c)
      printf(" %ld", reducedSubmatrix->columns[c]+1);
    printf("\n");
  }

  if (outputReducedMatrix)
  {
    CMR_CHRMAT* reducedMatrix = NULL;
    CMR_CALL( CMRchrmatZoomSubmat(cmr, matrix, reducedSubmatrix, &reducedMatrix) );
    fprintf(stderr, "\nExtracted reduced %lux%lu matrix with %lu nonzeros.\n", reducedMatrix->numRows,
      reducedMatrix->numColumns, reducedMatrix->numNonzeros);
    if (outputFormat == FILEFORMAT_MATRIX_DENSE)
      CMR_CALL( CMRchrmatPrintDense(cmr, reducedMatrix, stdout, '0', false) );
    else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
      CMR_CALL( CMRchrmatPrintSparse(cmr, reducedMatrix, stdout) );
    CMR_CALL( CMRchrmatFree(cmr, &reducedMatrix) );
  }

  if (violatorSubmatrix && outputNonSPElements)
  {
    fprintf(stderr, "\nMinimal non-series-parallel submatrix consists of these elements of the input matrix:\n");
    printf("%ld rows:", violatorSubmatrix->numRows);
    for (size_t r = 0; r < violatorSubmatrix->numRows; ++r)
      printf(" %ld", violatorSubmatrix->rows[r]+1);
    printf("\n%ld columns: ", violatorSubmatrix->numColumns);
    for (size_t c = 0; c < violatorSubmatrix->numColumns; ++c)
      printf(" %ld", violatorSubmatrix->columns[c]+1);
    printf("\n");
  }

  if (violatorSubmatrix && outputNonSPMatrix)
  {
    CMR_CHRMAT* violatorMatrix = NULL;
    CMR_CALL( CMRchrmatZoomSubmat(cmr, matrix, violatorSubmatrix, &violatorMatrix) );
    fprintf(stderr, "\nMinimal %lux%lu non-series parallel matrix with %lu nonzeros.\n", violatorMatrix->numRows,
      violatorMatrix->numColumns, violatorMatrix->numNonzeros);
    if (outputFormat == FILEFORMAT_MATRIX_DENSE)
      CMR_CALL( CMRchrmatPrintDense(cmr, violatorMatrix, stdout, '0', false) );
    else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
      CMR_CALL( CMRchrmatPrintSparse(cmr, violatorMatrix, stdout) );
    CMR_CALL( CMRchrmatFree(cmr, &violatorMatrix) );
  }
  
  /* Cleanup. */

  CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
  CMR_CALL( CMRsubmatFree(cmr, &reducedSubmatrix) );
  CMR_CALL( CMRfreeBlockArray(cmr, &reductions) );
  CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

int main(int argc, char** argv)
{
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  FileFormat outputFormat = FILEFORMAT_MATRIX_DENSE;
  char* instanceFileName = NULL;
  bool outputReductions = false;
  bool outputReducedElements = false;
  bool outputReducedMatrix = false;
  bool outputNonSPElements = false;
  bool outputNonSPMatrix = false;
  bool printStats = false;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-sp"))
      outputReductions = true;
    else if (!strcmp(argv[a], "-r"))
      outputReducedElements = true;
    else if (!strcmp(argv[a], "-R"))
      outputReducedMatrix = true;
    else if (!strcmp(argv[a], "-n"))
      outputNonSPElements = true;
    else if (!strcmp(argv[a], "-N"))
      outputNonSPMatrix = true;
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

  CMR_ERROR error = matrixSeriesParallel2Sums(instanceFileName, inputFormat, outputFormat, outputReductions,
    outputReducedElements, outputReducedMatrix, outputNonSPElements, outputNonSPMatrix, printStats);
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
