#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/regular.h>

typedef enum
{
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

/**
 * \brief Tests matrix from a file for regularity.
 */

static
CMR_ERROR testRegularity(
  const char* inputMatrixFileName,  /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,           /**< Format of the input matrix. */
  const char* outputTreeFileName,   /**< File name to print decomposition tree to, or \c NULL. */
  const char* outputMinorFileName,  /**< File name to print non-regular minor to, or \c NULL. */
  bool printStats,                  /**< Whether to print statistics to stderr. */
  bool directGraphicness,           /**< Whether to use fast graphicness routines. */
  bool seriesParallel,              /**< Whether to allow series-parallel operations in the decomposition tree. */
  bool simpleThreeSeparations,      /**< Whether to test for simple 3-separations. */
  int decomposeStrategy,            /**< Which strategy to use for 3-separations. */
  double timeLimit                  /**< Time limit to impose. */
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

  bool isRegular;
  CMR_SEYMOUR_NODE* decomposition = NULL;
  CMR_MINOR* minor = NULL;
  CMR_REGULAR_PARAMS params;
  CMR_CALL( CMRregularParamsInit(&params) );
  params.seymour.stopWhenIrregular = !outputTreeFileName;
  params.seymour.directGraphicness = directGraphicness;
  params.seymour.seriesParallel = seriesParallel;
  params.seymour.simpleThreeSeparations = simpleThreeSeparations;
  params.seymour.decomposeStrategy = decomposeStrategy;
  CMR_REGULAR_STATS stats;
  CMR_CALL( CMRregularStatsInit(&stats) );
  CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, outputTreeFileName ? &decomposition : NULL,
    outputMinorFileName ? &minor : NULL, &params, &stats, timeLimit) );

  fprintf(stderr, "Matrix %sregular.\n", isRegular ? "IS " : "IS NOT ");
  if (printStats)
    CMR_CALL( CMRregularStatsPrint(stderr, &stats, NULL) );

  if (decomposition)
    CMR_CALL( CMRseymourPrint(cmr, decomposition, stderr, true, true, true, true, true, true) );

  if (minor && outputMinorFileName)
  {
    bool outputMinorToFile = strcmp(outputMinorFileName, "-");
    fprintf(stderr, "Writing minimal non-regular submatrix to %s%s%s.\n", outputMinorToFile ? "file <" : "",
      outputMinorToFile ? outputMinorFileName : "stdout", outputMinorToFile ? ">" : "");    

    CMR_CALL( CMRminorWriteToFile(cmr, minor, matrix->numRows, matrix->numColumns, outputMinorFileName) );
  }

  /* Cleanup. */

  if (decomposition)
    CMR_CALL( CMRseymourRelease(cmr, &decomposition) );
  CMR_CALL( CMRminorFree(cmr, &minor) );
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
  fputs("  determines whether the matrix given in file IN-MAT is regular.\n", stderr);
  fputs("\n", stderr);

  fputs("Options:\n", stderr);
  fputs("  -i FORMAT    Format of file IN-MAT; default: dense.\n", stderr);
  fputs("  -D OUT-DEC   Write a decomposition tree of the regular matroid to file OUT-DEC; default: skip"
    " computation.\n", stderr);
  fputs("  -N NON-MINOR Write a minimal non-regular minor to file NON-MINOR; default: skip computation.\n", stderr);
  fputs("\n", stderr);

  fputs("Advanced options:\n", stderr);
  fputs("  --stats              Print statistics about the computation to stderr.\n", stderr);
  fputs("  --time-limit LIMIT   Allow at most LIMIT seconds for the computation.\n", stderr);
  fputs("  --decompose STRATEGY Strategy for decomposing among {DP, YP, P3, D3, Y3}; default: D3.\n", stderr);
  fputs("  --no-direct-graphic  Check only 3-connected matrices for regularity.\n", stderr);
  fputs("  --no-series-parallel Do not allow series-parallel operations in decomposition tree.\n", stderr);
  fputs("  --no-simple-3-sepa   Do not allow testing for simple 3-separations.\n", stderr);
  fputs("\n", stderr);

  fputs("Decomposition strategies: 1st letter for distributed, 2nd for concentrated rank(s).\n", stderr);
  fputs("  D Delta-sum (distributed ranks)\n", stderr);
  fputs("  Y Y-sum (distributed ranks)\n", stderr);
  fputs("  3 3-sum (concentrated rank)\n", stderr);
  fputs("  P pivot (changes rank type)\n", stderr);
  fputs("Note that D3 and Y3 do not produce pivots.\n\n", stderr);

  fputs("Formats for matrices: dense, sparse\n", stderr);
  fputs("If IN-MAT is `-' then the matrix is read from stdin.\n", stderr);
  fputs("If OUT-DEC or NON-MINOR is `-' then the decomposition tree (resp. the minor) is written to stdout.\n", stderr);

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  char* inputMatrixFileName = NULL;
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  char* outputTree = NULL;
  char* outputMinor = NULL;
  bool printStats = false;
  bool directGraphicness = true;
  bool seriesParallel = true;
  bool simpleThreeSeparations = true;
  double timeLimit = DBL_MAX;
  int decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_DELTASUM
    | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_THREESUM;
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
    else if (!strcmp(argv[a], "-D") && a+1 < argc)
      outputTree = argv[++a];
    else if (!strcmp(argv[a], "-N") && a+1 < argc)
      outputMinor = argv[++a];
    else if (!strcmp(argv[a], "--stats"))
      printStats = true;
    else if (!strcmp(argv[a], "--no-direct-graphic"))
      directGraphicness = false;
    else if (!strcmp(argv[a], "--no-series-parallel"))
      seriesParallel = false;
    else if (!strcmp(argv[a], "--no-simple-3-sepa"))
      simpleThreeSeparations = false;
    else if (!strcmp(argv[a], "--decompose") && a+1 < argc)
    {
      ++a;
      if (!strcasecmp(argv[a], "DP"))
      {
        decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_DELTASUM
          | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_PIVOT;
      }
      else if (!strcasecmp(argv[a], "YP"))
      {
        decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_YSUM
          | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_PIVOT;
      }
      else if (!strcasecmp(argv[a], "P3"))
      {
        decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_PIVOT
          | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_THREESUM;
      }
      else if (!strcasecmp(argv[a], "D3"))
      {
        decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_DELTASUM
          | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_THREESUM;
      }
      else if (!strcasecmp(argv[a], "Y3"))
      {
        decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_YSUM
          | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_THREESUM;
      }
      else
      {
        fprintf(stderr, "Error: Invalid decomposition strategy <%s> specified.\n\n", argv[a]);
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
    else if (!inputMatrixFileName)
      inputMatrixFileName = argv[a];
    else
    {
      printf("Error: Two input files <%s> and <%s> specified.\n\n", inputMatrixFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (!inputMatrixFileName)
  {
    puts("Error: No input file specified.\n");
    return printUsage(argv[0]);
  }

  CMR_ERROR error;
  error = testRegularity(inputMatrixFileName, inputFormat, outputTree, outputMinor, printStats, directGraphicness,
    seriesParallel, simpleThreeSeparations, decomposeStrategy, timeLimit);

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
