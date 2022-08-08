#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <cmr/matrix.h>
#include <cmr/tu.h>

typedef enum
{
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

/**
 * \brief Tests matrix from a file for total unimodularity.
 */

static
CMR_ERROR testTotalUnimodularity(
  const char* inputMatrixFileName,      /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,               /**< Format of the input matrix. */
  const char* outputTreeFileName,       /**< File name to print decomposition tree to, or \c NULL. */
  const char* outputSubmatrixFileName,  /**< File name to print non-TU submatrix to, or \c NULL. */
  bool printStats,                      /**< Whether to print statistics to stderr. */
  bool directGraphicness,               /**< Whether to use fast graphicness routines. */
  bool seriesParallel                   /**< Whether to allow series-parallel operations in the decomposition tree. */
)
{
  clock_t readClock = clock();
  FILE* inputMatrixFile = strcmp(inputMatrixFileName, "-") ? fopen(inputMatrixFileName, "r") : stdin;
  if (!inputMatrixFile)
  {
    fprintf(stderr, "Unable to open file <%s>\n", inputMatrixFileName);
    return CMR_ERROR_INPUT;
  }

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_CHRMAT* matrix = NULL;
  CMR_ERROR error;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
  {
    error = CMRchrmatCreateFromDenseStream(cmr, inputMatrixFile, &matrix);
    if (error == CMR_ERROR_INPUT)
      fprintf(stderr, "Error when reading dense matrix from <%s>: %s\n", inputMatrixFileName, CMRgetErrorMessage(cmr));
  }
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
  {
    error = CMRchrmatCreateFromSparseStream(cmr, inputMatrixFile, &matrix);
    if (error == CMR_ERROR_INPUT)
      fprintf(stderr, "Error when reading dense matrix from <%s>: %s\n", inputMatrixFileName, CMRgetErrorMessage(cmr));
  }
  else
    assert(false);

  if (inputMatrixFile != stdin)
    fclose(inputMatrixFile);

  CMR_CALL(error);

  fprintf(stderr, "Read %lux%lu matrix with %lu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Actual test. */

  bool isTU;
  CMR_DEC* decomposition = NULL;
  CMR_SUBMAT* submatrix = NULL;
  CMR_TU_PARAMETERS params;
  CMR_CALL( CMRparamsTotalUnimodularityInit(&params) );
  params.regular.completeTree = outputTreeFileName;
  params.regular.matrices = outputTreeFileName ? CMR_DEC_CONSTRUCT_ALL : CMR_DEC_CONSTRUCT_NONE;
  params.regular.directGraphicness = directGraphicness;
  params.regular.seriesParallel = seriesParallel;
  CMR_TU_STATISTICS stats;
  CMR_CALL( CMRstatsTotalUnimodularityInit(&stats));
  CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, outputTreeFileName ? &decomposition : NULL,
    outputSubmatrixFileName ? &submatrix : NULL, &params, &stats) );

  printf("Matrix %stotally unimodular.\n", isTU ? "IS " : "IS NOT ");
  if (printStats)
    CMR_CALL( CMRstatsTotalUnimodularityPrint(stderr, &stats, NULL) );

  if (decomposition)
  {
    // TODO: Write decomposition.
    assert(!"NOT IMPLEMENTED");
    exit(EXIT_FAILURE);
  }

  if (submatrix && outputSubmatrixFileName)
  {
    bool outputSubmatrixToFile = strcmp(outputSubmatrixFileName, "-");
    fprintf(stderr, "Writing minimal non-totally-unimodular submatrix to %s%s%s.\n", outputSubmatrixToFile ? "file <" : "",
      outputSubmatrixToFile ? outputSubmatrixFileName : "stdout", outputSubmatrixToFile ? ">" : "");    

    CMR_CALL( CMRsubmatWriteToFile(cmr, submatrix, matrix->numRows, matrix->numColumns, outputSubmatrixFileName) );
  }

  /* Cleanup. */

  CMR_CALL( CMRdecFree(cmr, &decomposition) );
  CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
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
  fputs("  determines whether the matrix given in file IN-MAT is totally unimodular.\n\n", stderr);
  fputs("Options:\n", stderr);
  fputs("  -i FORMAT  Format of file IN-MAT, among `dense' and `sparse'; default: dense.\n", stderr);
  fputs("  -D OUT-DEC Write a decomposition tree of the underlying regular matroid to file OUT-DEC; default: skip computation.\n", stderr);
  fputs("  -N NON-SUB Write a minimal non-totally-unimodular submatrix to file NON-SUB; default: skip computation.\n", stderr);
  fputs("  -s         Print statistics about the computation to stderr.\n\n", stderr);
  fputs("Advanced options:\n", stderr);
  fputs("  --no-direct-graphic  Check only 3-connected matrices for regularity.\n", stderr);
  fputs("  --no-series-parallel Do not allow series-parallel operations in decomposition tree.\n\n", stderr);
  fputs("If IN-MAT is `-' then the matrix is read from stdin.\n", stderr);
  fputs("If OUT-DEC or NON-SUB is `-' then the decomposition tree (resp. the submatrix) is written to stdout.\n", stderr);

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  char* inputMatrixFileName = NULL;
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  char* outputTree = NULL;
  char* outputSubmatrix = NULL;
  bool printStats = false;
  bool directGraphicness = true;
  bool seriesParallel = true;
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
    else if (!strcmp(argv[a], "-D") && a+1 < argc)
      outputTree = argv[++a];
    else if (!strcmp(argv[a], "-N") && a+1 < argc)
      outputSubmatrix = argv[++a];
    else if (!strcmp(argv[a], "-s"))
      printStats = true;
    else if (!strcmp(argv[a], "--no-direct-graphic"))
      directGraphicness = false;
    else if (!strcmp(argv[a], "--no-series-parallel"))
      seriesParallel = false;
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
  error = testTotalUnimodularity(inputMatrixFileName, inputFormat, outputTree, outputSubmatrix, printStats,
    directGraphicness, seriesParallel);

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
