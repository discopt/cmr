#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

#include <cmr/matrix.h>
#include <cmr/ctu.h>

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
  puts("Determines whether matrix in FILE is complement totally unimodular or applies row- or column-complement operations (or both at the same time).");
  puts("Options:");
  puts("  -i FORMAT  Format of input FILE; default: `dense'.");
  puts("  -o FORMAT  Format of output matrices; default: `dense'.");
  puts("  -r ROW     Perform a row complement operation; do not test for complement total unimodularity.");
  puts("  -c COLUMN  Perform a row complement operation; do not test for complement total unimodularity.");
  puts("  -n         Output a complement operations that leads to a non-totally-unimodular matrix.");
  puts("  -N         Output a complemented matrix that is non-totally-unimodular.");
  puts("  -s         Print statistics about the computation to stderr.");
  puts("Formats for matrices: dense, sparse");
  puts("If FILE is `-', then the input will be read from stdin.");

  return EXIT_FAILURE;
}

static
CMR_ERROR complementMatrix(
  const char* instanceFileName, /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,       /**< Format of the input matrix. */
  FileFormat outputFormat,      /**< Format of the output submatrix. */
  size_t complementRow,         /**< Complement row. */
  size_t complementColumn       /**< Complement column. */
)
{
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

  /* Carry out complementing. */
  CMR_CHRMAT* complemented = NULL;
  CMR_CALL( CMRcomplementRowColumn(cmr, matrix, complementRow, complementColumn, &complemented) );
  if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatPrintDense(cmr, complemented, stdout, '0', false) );
  else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatPrintSparse(cmr, complemented, stdout) );

  CMR_CALL( CMRchrmatFree(cmr, &complemented) );
  
  /* Cleanup. */

  CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

/**
 * \brief Tests matrix from a file for total unimodularity.
 */

static
CMR_ERROR testComplementTotalUnimodularity(
  const char* instanceFileName, /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,       /**< Format of the input matrix. */
  FileFormat outputFormat,      /**< Format of the output submatrix. */
  bool outputComplements,       /**< Whether to print the complement operations leading to a non-TU matrix. */
  bool outputComplemented,      /**< Whether to print the complemented non-TU matrix. */
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

  bool isCTU;
  size_t complementRow = SIZE_MAX;
  size_t complementColumn = SIZE_MAX;
  CMR_CTU_STATISTICS stats;
  CMR_CALL( CMRstatsComplementTotalUnimodularityInit(&stats) );
  CMR_CALL( CMRtestComplementTotalUnimodularity(cmr, matrix, &isCTU, &complementRow, &complementColumn, &stats) );

  fprintf(stderr, "Matrix %scomplement totally unimodular.\n", isCTU ? "IS " : "IS NOT ");
  if (printStats)
    CMR_CALL( CMRstatsComplementTotalUnimodularityPrint(stderr, &stats, NULL) );

  if (complementRow < SIZE_MAX || complementColumn < SIZE_MAX)
  {
    if (outputComplements)
    {
      fprintf(stderr, "\nNon-TU matrix is obtained by");
      if (complementRow < SIZE_MAX)
      {
        fprintf(stderr, " complementing row %lu", complementRow + 1);
        if (complementColumn < SIZE_MAX)
          fprintf(stderr, " and complementing column %lu", complementColumn + 1);
      }
      else
        fprintf(stderr, " complementing column %lu", complementColumn + 1);
      fprintf(stderr, ".\n");
    }

    if (outputComplemented)
    {
      CMR_CHRMAT* complemented = NULL;
      CMR_CALL( CMRcomplementRowColumn(cmr, matrix, complementRow, complementColumn, &complemented) );

      fprintf(stderr, "\nComplemented non-TU matrix with %lu nonzeros:\n", complemented->numNonzeros);
      if (outputFormat == FILEFORMAT_MATRIX_DENSE)
        CMR_CALL( CMRchrmatPrintDense(cmr, complemented, stdout, '0', false) );
      else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
        CMR_CALL( CMRchrmatPrintSparse(cmr, complemented, stdout) );

      CMR_CALL( CMRchrmatFree(cmr, &complemented) );
    }
  }

  /* Cleanup. */

  CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

int main(int argc, char** argv)
{
  FileFormat inputFormat = FILEFORMAT_UNDEFINED;
  FileFormat outputFormat = FILEFORMAT_MATRIX_DENSE;
  size_t complementRow = SIZE_MAX;
  size_t complementColumn = SIZE_MAX;
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
    else if (!strcmp(argv[a], "-r") && a+1 < argc)
    {
      char* p;
      complementRow = strtoull(argv[a+1], &p, 10);
      if (*p != '\0' || complementRow <= 0)
      {
        fprintf(stderr, "Error: invalid complement row <%s>", argv[a+1]);
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a++;
    }
    else if (!strcmp(argv[a], "-c") && a+1 < argc)
    {
      char* p;
      complementColumn = strtoull(argv[a+1], &p, 10);
      if (*p != '\0' || complementColumn <= 0)
      {
        fprintf(stderr, "Error: invalid complement column <%s>", argv[a+1]);
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a++;
    }
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
  if ((outputSubmatrix || outputSubmatrixElements) && (complementRow < SIZE_MAX || complementColumn < SIZE_MAX))
  {
    puts("Error: row/column complementing and testing for complement total unimodularity at the same time is not allowed.");
    return printUsage(argv[0]);
  }

  if (inputFormat == FILEFORMAT_UNDEFINED)
    inputFormat = FILEFORMAT_MATRIX_DENSE;

  CMR_ERROR error;
  if (complementRow < SIZE_MAX || complementColumn < SIZE_MAX)
  {
    error = complementMatrix(instanceFileName, inputFormat, outputFormat, complementRow, complementColumn);
  }
  else
  {
    error = testComplementTotalUnimodularity(instanceFileName, inputFormat, outputFormat, outputSubmatrixElements,
      outputSubmatrix, printStats);
  }

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
