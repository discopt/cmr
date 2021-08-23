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
  puts("  -s         Output the elements of a minimal non-totally-unimodular submatrix.");
  puts("  -S         Output a minimal non-totally-unimodular submatrix.");
  puts("Formats for matrices: dense, sparse");
  puts("If FILE is `-', then the input will be read from stdin.");

  return EXIT_FAILURE;
}

/**
 * \brief Tests matrix from a file for total unimodularity.
 */

static
CMR_ERROR testTotalUnimodularity(
  const char* instanceFileName, /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,       /**< Format of the input matrix. */
  FileFormat outputFormat,      /**< Format of the output submatrix. */
  bool printTree,               /**< Whether to print the decomposition tree. */
  bool printSubmatrixElements,  /**< Whether to print the elements of a non-TU submatrix. */
  bool printSubmatrix           /**< Whether to print a non-TU submatrix. */
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

  /* Actual test. */

  bool isTU;
  CMR_TU_DEC* decomposition = NULL;
  CMR_SUBMAT* submatrix = NULL;
  startClock = clock();
  CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, printTree ? &decomposition : NULL,
    (printSubmatrix || printSubmatrixElements) ? &submatrix : NULL) );

  fprintf(stderr, "Determined in %f seconds that it is %stotally unimodular.\n",
    (clock() - startClock) * 1.0 / CLOCKS_PER_SEC, isTU ? "" : "NOT ");

  if (decomposition)
  {
    CMR_CALL( CMRtudecFree(cmr, &decomposition) );
  }

  if (submatrix)
  {
    if (printSubmatrixElements)
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

    if (printSubmatrix)
    {
      startClock = clock();
      CMR_CHRMAT* violatorMatrix = NULL;
      CMR_CALL( CMRchrmatFilterSubmat(cmr, matrix, submatrix, &violatorMatrix) );
      endTime = clock();
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

int main(int argc, char** argv)
{
  FileFormat inputFormat = FILEFORMAT_UNDEFINED;
  FileFormat outputFormat = FILEFORMAT_MATRIX_DENSE;
  bool printTree = false;
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
    else if (!strcmp(argv[a], "-d"))
      printTree = true;
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
  error = testTotalUnimodularity(instanceFileName, inputFormat, outputFormat, printTree, printSubmatrixElements, printSubmatrix);

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
