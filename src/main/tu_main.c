#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/tu.h>
#include <cmr/linear_algebra.h>

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
  bool planarityCheck,                  /**< Whether to check for planarity. */
  bool seriesParallel,                  /**< Whether to allow series-parallel operations in the decomposition tree. */
  int decomposeStrategy,                /**< Which strategy to use for 3-separations. */
  bool naiveSubmatrix,                  /**< Use naive bad submatrix heuristic instead of greedy algorithm. */
  CMR_TU_ALGORITHM algorithm,           /**< Algorithm to use for TU test. */
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

  fprintf(stderr, "Read %zux%zu matrix with %zu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Actual test. */

  bool isTU;
  CMR_SEYMOUR_NODE* decomposition = NULL;
  CMR_SUBMAT* submatrix = NULL;
  CMR_TU_PARAMS params;
  CMR_CALL( CMRtuParamsInit(&params) );
  params.algorithm = algorithm;
  params.seymour.stopWhenIrregular = !outputTreeFileName;
  params.seymour.directGraphicness = directGraphicness;
  params.seymour.seriesParallel = seriesParallel;
  params.seymour.planarityCheck = planarityCheck;
  params.seymour.decomposeStrategy = decomposeStrategy;
  params.naiveSubmatrix = naiveSubmatrix;
  CMR_TU_STATS stats;
  CMR_CALL( CMRtuStatsInit(&stats));
  error = CMRtuTest(cmr, matrix, &isTU, outputTreeFileName ? &decomposition : NULL,
    outputSubmatrixFileName ? &submatrix : NULL, &params, &stats, timeLimit);
  if (error == CMR_ERROR_TIMEOUT)
  {
    fprintf(stderr, "Time limit exceeded!\n");
    CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    if (printStats)
      CMR_CALL( CMRtuStatsPrint(stderr, &stats, NULL) );

    CMR_CALL( CMRfreeEnvironment(&cmr) );
    return error;
  }
  CMR_CALL( error );

  printf("Matrix %stotally unimodular.\n", isTU ? "IS " : "IS NOT ");
  if (printStats)
    CMR_CALL( CMRtuStatsPrint(stderr, &stats, NULL) );

  if (decomposition)
  {
    FILE* outputTreeFile;
    if (strcmp(outputTreeFileName, "-"))
      outputTreeFile = fopen(outputTreeFileName, "w");
    else
      outputTreeFile = stdout;
    CMR_CALL( CMRseymourPrint(cmr, decomposition, outputTreeFile, true, true, true, true, true, true) );
    if (outputTreeFile != stdout)
      fclose(outputTreeFile);
  }

  if (submatrix && outputSubmatrixFileName)
  {
    /* Extract submatrix to compute its determinant. */
    CMR_CHRMAT* violator = NULL;
    int64_t determinant = 0;
    CMR_CALL( CMRchrmatSlice(cmr, matrix, submatrix, &violator) );
    CMR_CALL( CMRchrmatDeterminant(cmr, violator, &determinant) );
    CMR_CALL( CMRchrmatFree(cmr, &violator) );

    bool outputSubmatrixToFile = strcmp(outputSubmatrixFileName, "-");
    fprintf(stderr, "Writing minimal non-totally-unimodular submatrix with absolute determinant %ld to %s%s%s.\n", determinant,
      outputSubmatrixToFile ? "file <" : "", outputSubmatrixToFile ? outputSubmatrixFileName : "stdout",
      outputSubmatrixToFile ? ">" : "");

    CMR_CALL( CMRsubmatWriteToFile(cmr, submatrix, matrix->numRows, matrix->numColumns, outputSubmatrixFileName) );
  }

  /* Cleanup. */

  if (decomposition)
    CMR_CALL( CMRseymourRelease(cmr, &decomposition) );
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

  fprintf(stderr, "%s IN-MAT [OPTION]...\n", program);
  fputs("  determines whether the matrix given in file IN-MAT is totally unimodular.\n", stderr);
  fputs("\n", stderr);

  fputs("Options:\n", stderr);
  fputs("  -i FORMAT  Format of file IN-MAT; default: dense.\n", stderr);
  fputs("  -D OUT-DEC Write a decomposition tree of the underlying regular matroid to file OUT-DEC; "
    "default: skip computation.\n", stderr);
  fputs("  -N NON-SUB Write a minimal non-totally-unimodular submatrix to file NON-SUB; default: skip computation.\n",
    stderr);
  fputs("\n", stderr);

  fputs("Advanced options:\n", stderr);
  fputs("  --stats              Print statistics about the computation to stderr.\n", stderr);
  fputs("  --time-limit LIMIT   Allow at most LIMIT seconds for the computation.\n", stderr);
  fputs("  --decompose          Strategy for decomposing; among {seymour, truemper, pivotless}.\n", stderr);
  fputs("  --no-direct-graphic  Check only 3-connected matrices for regularity.\n", stderr);
  fputs("  --no-planarity       Do not test for planarity.\n", stderr);
  fputs("  --no-series-parallel Do not allow series-parallel operations in decomposition tree.\n", stderr);
  fputs("  --naive-submatrix    Use naive bad submatrix algorithm instead of greedy heuristic.\n", stderr);
  fputs("  --algo ALGO          Use algorithm from {decomposition, eulerian, partition}; default: decomposition.\n",
    stderr);
  fputs("\n", stderr);

  fputs("Formats for matrices: dense, sparse\n", stderr);
  fputs("If IN-MAT is `-' then the matrix is read from stdin.\n", stderr);
  fputs("If OUT-DEC or NON-SUB is `-' then the decomposition tree (resp. the submatrix) is written to stdout.\n",
    stderr);

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
  bool planarityCheck = true;
  bool seriesParallel = true;
  bool naiveSubmatrix = true;
  int decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR;
  double timeLimit = DBL_MAX;
  CMR_TU_ALGORITHM algorithm = CMR_TU_ALGORITHM_DECOMPOSITION;
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
    else if (!strcmp(argv[a], "--stats"))
      printStats = true;
    else if (!strcmp(argv[a], "--no-direct-graphic"))
      directGraphicness = false;
    else if (!strcmp(argv[a], "--no-planarity"))
      planarityCheck = false;
    else if (!strcmp(argv[a], "--no-series-parallel"))
      seriesParallel = false;
    else if (!strcmp(argv[a], "--naive-submatrix"))
      naiveSubmatrix = true;
    else if (!strcmp(argv[a], "--decompose") && a+1 < argc)
    {
      ++a;
      if (!strcasecmp(argv[a], "seymour"))
        decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR;
      else if (!strcasecmp(argv[a], "truemper"))
        decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_TRUEMPER;
      else if (!strcasecmp(argv[a], "pivotless"))
        decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_PIVOTLESS;
      else
      {
        fprintf(stderr, "Error: Invalid threesum-strategy <%s> specified.\n\n", argv[a]);
        return printUsage(argv[0]);
      }
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
    else if (!strcmp(argv[a], "--algo") && (a+1 < argc))
    {
      if (!strcmp(argv[a+1], "decomposition"))
        algorithm = CMR_TU_ALGORITHM_DECOMPOSITION;
      else if (!strcmp(argv[a+1], "eulerian"))
        algorithm = CMR_TU_ALGORITHM_EULERIAN;
      else if (!strcmp(argv[a+1], "partition"))
        algorithm = CMR_TU_ALGORITHM_PARTITION;
      else
      {
        fprintf(stderr, "Error: Invalid algorithm <%s> specified.\n\n", argv[a+1]);
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
  error = testTotalUnimodularity(inputMatrixFileName, inputFormat, outputTree, outputSubmatrix, printStats,
    directGraphicness, planarityCheck, seriesParallel, decomposeStrategy, naiveSubmatrix, algorithm, timeLimit);

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
