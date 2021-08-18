#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <cmr/matrix.h>
#include <cmr/camion.h>

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
  puts("Checks whether matrix in FILE is Camion-signed.");
  puts("Options:");
  puts("  -i FORMAT  Format of input FILE; default: `dense'.");
  puts("  -o FORMAT  Format of output; default: `dense'.");
  puts("  -s         Output the elements of a minimal non-camion submatrix.");
  puts("  -S         Output a minimal non-camion submatrix.");
  puts("Formats for matrices: dense, sparse");
  puts("If FILE is `-', then the input will be read from stdin.");
  return EXIT_FAILURE;
}

/**
 * \brief Tests matrix from a file for total unimodularity.
 */

static
CMR_ERROR testCamionSigned(
  const char* instanceFileName, /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,       /**< Format of the input matrix. */
  FileFormat outputFormat,      /**< Format of the output submatrix. */
  bool printSubmatrixElements,  /**< Whether to print the elements of a non-camion submatrix. */
  bool printSubmatrix           /**< Whether to print a non-camion submatrix. */
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
    CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, &matrix, instanceFile) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, &matrix, instanceFile) );
  if (instanceFile != stdin)
    fclose(instanceFile);
  fprintf(stderr, "Read %dx%d matrix with %d nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - startClock) * 1.0 / CLOCKS_PER_SEC);

  /* Actual test. */

  bool isCamion;
  CMR_SUBMAT* submatrix = NULL;
  startClock = clock();
  CMR_CALL( CMRtestCamionSigned(cmr, matrix, &isCamion,
    (printSubmatrix || printSubmatrixElements) ? &submatrix : NULL) );

  fprintf(stderr, "Determined in %f seconds that it is %sCamion-signed.\n",
    (clock() - startClock) * 1.0 / CLOCKS_PER_SEC, isCamion ? "" : "NOT ");

  if (submatrix)
  {
    if (printSubmatrixElements)
    {
      fprintf(stderr, "\nNon-camion submatrix consists of these elements:\n");
      printf("%ld rows:", submatrix->numRows);
      for (size_t r = 0; r < submatrix->numRows; ++r)
        printf(" %ld", submatrix->rows[r]+1);
      printf("\n%ld columns: ", submatrix->numColumns);
      for (size_t c = 0; c < submatrix->numColumns; ++c)
        printf(" %ld", submatrix->columns[c]+1);
      printf("\n");
    }

    if (printSubmatrix)
    {
      startClock = clock();
      CMR_CHRMAT* violatorMatrix = NULL;
      CMR_CALL( CMRchrmatFilterSubmat(cmr, matrix, submatrix, &violatorMatrix) );
      endTime = clock();
      fprintf(stderr, "\nExtracted %dx%d non-camion submatrix with %d nonzeros in %f seconds.\n", violatorMatrix->numRows,
        violatorMatrix->numColumns, violatorMatrix->numNonzeros, (endTime - startClock) * 1.0 / CLOCKS_PER_SEC );
      if (outputFormat == FILEFORMAT_MATRIX_DENSE)
        CMR_CALL( CMRchrmatPrintDense(cmr, stdout, violatorMatrix, '0', false) );
      else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
        CMR_CALL( CMRchrmatPrintSparse(stdout, violatorMatrix) );
      CMR_CALL( CMRchrmatFree(cmr, &violatorMatrix) );
    }

    CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
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
  bool printSubmatrixElements = false;
  bool printSubmatrix = false;
  char* instanceFileName = NULL;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-s"))
      printSubmatrixElements = true;
    else if (!strcmp(argv[a], "-S"))
      printSubmatrix = true;
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
  error = testCamionSigned(instanceFileName, inputFormat, outputFormat, printSubmatrixElements, printSubmatrix);

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
