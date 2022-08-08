#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <cmr/matrix.h>
#include <cmr/series_parallel.h>

typedef enum
{
  FILEFORMAT_MATRIX_DENSE = 1,
  FILEFORMAT_MATRIX_SPARSE = 2
} FileFormat;

CMR_ERROR recognizeSeriesParallel(
  const char* inputMatrixFileName,      /**< File name of input matrix (may be `-` for stdin). */
  FileFormat inputFormat,               /**< Format of input matrix. */
  const char* outputReductionsFileName, /**< File name for output reductions (may be `-` for stdout). */
  const char* outputReducedFileName,    /**< File name for output of reduced matrix (may be `-` for stdout). */
  const char* outputSubmatrixFileName,  /**< File name for minimal non-series-parallel submatrix (may be `-` for stdout). */
  bool printStats                       /**< Whether to print statistics to stderr. */
)
{
  clock_t readClock = clock();
  FILE* inputMatrixFile = strcmp(inputMatrixFileName, "-") ? fopen(inputMatrixFileName, "r") : stdin;
  if (!inputMatrixFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_CHRMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, inputMatrixFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, inputMatrixFile, &matrix) );
  if (inputMatrixFile != stdin)
    fclose(inputMatrixFile);
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
    outputReducedFileName ? &reducedSubmatrix : NULL, outputSubmatrixFileName ? &violatorSubmatrix : NULL, &stats) );

  fprintf(stderr, "Matrix %sseries-parallel. %ld reductions can be applied.\n",
    numReductions == matrix->numRows + matrix->numColumns ? "IS " : "is NOT ", numReductions);
  if (printStats)
    CMR_CALL( CMRstatsSeriesParallelPrint(stderr, &stats, NULL) );

  if (outputReductionsFileName)
  {
    bool outputReductionsToFile = strcmp(outputReductionsFileName, "-");
    FILE* outputReductionsFile = outputReductionsToFile ? fopen(outputReductionsFileName, "w") : stdout;
    fprintf(stderr, "Writing %ld series-parallel reductions to %s%s%s.\n", numReductions,
      outputReductionsToFile ? "file <" : "", outputReductionsToFile ? outputReductionsFileName : "stdout",
      outputReductionsToFile ? ">" : "");    

    fprintf(outputReductionsFile, "%ld\n", numReductions);
    for (size_t i = 0; i < numReductions; ++i)
      fprintf(outputReductionsFile, "%s\n", CMRspReductionString(reductions[i], NULL));

    if (outputReductionsToFile)
      fclose(outputReductionsFile);
  }

  if (outputReducedFileName)
  {
    bool outputReducedToFile = strcmp(outputReducedFileName, "-");
    fprintf(stderr, "Writing reduced submatrix to %s%s%s.\n", outputReducedToFile ? "file <" : "",
      outputReducedToFile ? outputReducedFileName : "stdout", outputReducedToFile ? ">" : "");    

    CMR_CALL( CMRsubmatWriteToFile(cmr, reducedSubmatrix, matrix->numRows, matrix->numColumns, outputReducedFileName) );
  }

  if (violatorSubmatrix && outputSubmatrixFileName)
  {
    bool outputSubmatrixToFile = strcmp(outputSubmatrixFileName, "-");
    fprintf(stderr, "Writing minimal non-series-parallel submatrix to %s%s%s.\n", outputSubmatrixToFile ? "file <" : "",
      outputSubmatrixToFile ? outputSubmatrixFileName : "stdout", outputSubmatrixToFile ? ">" : "");    

    CMR_CALL( CMRsubmatWriteToFile(cmr, violatorSubmatrix, matrix->numRows, matrix->numColumns,
      outputSubmatrixFileName) );
  }

  /* Cleanup. */

  CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
  CMR_CALL( CMRsubmatFree(cmr, &reducedSubmatrix) );
  CMR_CALL( CMRfreeBlockArray(cmr, &reductions) );
  CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

int printUsage(const char* program)
{
  fputs("Usage:\n", stderr);
  fprintf(stderr, "%s IN-MAT [OPTION]...\n\n", program);
  fputs("  determines whether the matrix given in file IN-MAT is series-parallel.\n\n", stderr);
  fputs("Options:\n", stderr);
  fputs("  -i FORMAT       Format of file IN-MAT, among `dense' and `sparse'; default: dense.\n", stderr);
  fputs("  -S OUT-SP       Write the list of series-parallel reductions to file OUT-SP; default: skip computation.\n", stderr);
  fputs("  -R OUT-REDUCED  Write the reduced submatrix to file `OUT-REDUCED`; default: skip computation.\n", stderr);
  fputs("  -N NON-SUB      Write a minimal non-series-parallel submatrix to file `NON-SUB`; default: skip computation.\n", stderr);
  fputs("  -s`             Print statistics about the computation to stderr.\n\n", stderr);
  fputs("If IN-MAT is `-' then the matrix is read from stdin.\n", stderr);
  fputs("If OUT-SP, OUT-REDUCED or NON-SUB is `-' then the list of reductions (resp. the submatrix) is written to stdout.\n", stderr);

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  char* inputMatrixFileName = NULL;
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  char* outputReductionsFileName = NULL;
  char* outputReducedFileName = NULL;
  char* outputSubmatrixFileName = NULL;
  bool printStats = false;
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
        printf("Error: unknown input file format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-S") && a+1 < argc)
      outputReductionsFileName = argv[++a];
    else if (!strcmp(argv[a], "-R") && a+1 < argc)
      outputReducedFileName = argv[++a];
    else if (!strcmp(argv[a], "-N") && a+1 < argc)
      outputSubmatrixFileName = argv[++a];
    else if (!strcmp(argv[a], "-s"))
      printStats = true;
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
    fputs("No input file specified.\n\n", stderr);
    return printUsage(argv[0]);
  }

  CMR_ERROR error = recognizeSeriesParallel(inputMatrixFileName, inputFormat, outputReductionsFileName,
    outputReducedFileName, outputSubmatrixFileName, printStats);

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
