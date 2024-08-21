#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <sys/time.h>
#include <stdint.h>
#include <math.h>

#include <cmr/graphic.h>

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,       /**< Whether the file format of output was defined by the user. */
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

static inline
size_t randRange(size_t first, size_t beyond)
{
  size_t N = beyond - first;
  size_t representatives = (RAND_MAX + 1u) / N;
  size_t firstInvalid = N * representatives;
  size_t x;
  do
  {
    x = rand();
  }
  while (x >= firstInvalid);
  return first + x / representatives;
}

int printUsage(const char* program)
{
  fprintf(stderr, "Usage: %s [OPTIONS] ORDER\n\n", program);
  fputs("Creates an ORDER-by-ORDER cycle matrix.\n", stderr);
  fputs("Options:\n", stderr);
  fputs("  -01        In each column with two 0s in rows 1 and 2, replace them by 1s; default: off.\n", stderr);
  fputs("  -o FORMAT  Format of output FILE; default: `dense'.\n", stderr);
  fputs("Formats for matrices: dense, sparse\n", stderr);
  return EXIT_FAILURE;
}

int compare(const void* pa, const void* pb)
{
  size_t a = *((size_t*)(pa));
  size_t b = *((size_t*)(pb));
  return a < b ? -1 : (a > b);
}

CMR_ERROR genMatrixCycle(
  size_t numRowsColumns,  /**< Number of rows and columns of matrix. */
  bool change01,          /**< In each column with two 0s in rows 1 and 2, replace them by 1s. */
  FileFormat outputFormat /**< Output file format. */
)
{
  clock_t startTime = clock();
  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  CMRchrmatCreate(cmr, &matrix, numRowsColumns, numRowsColumns, 4 * numRowsColumns);

  /* Create the nonzeros. */
  matrix->numNonzeros = 0;
  for (size_t row = 0; row < numRowsColumns; ++row)
  {
    matrix->rowSlice[row] = matrix->numNonzeros;
    if (row + 1 == numRowsColumns)
    {
      matrix->entryValues[matrix->numNonzeros] = 1;
      matrix->entryColumns[matrix->numNonzeros++] = 0;
    }
    matrix->entryValues[matrix->numNonzeros] = 1;
    matrix->entryColumns[matrix->numNonzeros++] = row;
    if (row + 1 < numRowsColumns)
    {
      matrix->entryValues[matrix->numNonzeros] = 1;
      matrix->entryColumns[matrix->numNonzeros++] = row + 1;
    }
    if (change01 && row < 2)
    {
      for (size_t c = 3; c < numRowsColumns; ++c)
      {
        matrix->entryValues[matrix->numNonzeros] = 1;
        matrix->entryColumns[matrix->numNonzeros++] = c;
      }
    }
  }
  matrix->rowSlice[numRowsColumns] = matrix->numNonzeros;

  double generationTime = (clock() - startTime) * 1.0 / CLOCKS_PER_SEC;
  fprintf(stderr, "Generated a %zux%zu matrix with %zu nonzeros in %f seconds.\n", numRowsColumns, numRowsColumns,
    matrix->numNonzeros, generationTime);

  /* Print matrix. */
  if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );
  else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatPrintSparse(cmr, matrix, stdout) );

  /* Cleanup. */
  CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

int main(int argc, char** argv)
{
  struct timeval curTime;
  gettimeofday(&curTime, NULL);
  srand(curTime.tv_usec);

  FileFormat outputFormat = FILEFORMAT_UNDEFINED;
  size_t numRowsColumns = SIZE_MAX;
  bool change01 = false;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-01"))
      change01 = true;
    else if (!strcmp(argv[a], "-o") && (a+1 < argc))
    {
      if (!strcmp(argv[a+1], "dense"))
        outputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        outputFormat = FILEFORMAT_MATRIX_SPARSE;
      else
      {
        fprintf(stderr, "Error: unknown output format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (numRowsColumns == SIZE_MAX)
    {
      char* p = NULL;
      numRowsColumns = strtoull(argv[a], &p, 10);
      if (*p != '\0')
      {
        fprintf(stderr, "Error: invalid number of rows/columns <%s>.\n\n", argv[a]);
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
    }
    else
    {
      fprintf(stderr, "Error: more than one size indicators specified: %zu %s\n\n", numRowsColumns, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (numRowsColumns == SIZE_MAX)
  {
    fputs("Error: no size indicator specified.\n", stderr);
    return printUsage(argv[0]);
  }
  else if (numRowsColumns <= 0)
  {
    fputs("Error: matrix must have at least 1 row/column.\n", stderr);
    return printUsage(argv[0]);
  }
  if (outputFormat == FILEFORMAT_UNDEFINED)
    outputFormat = FILEFORMAT_MATRIX_DENSE;

  CMR_ERROR error = genMatrixCycle(numRowsColumns, change01, outputFormat);
  switch (error)
  {
  case CMR_ERROR_INPUT:
    fputs("Input error.\n", stderr);
    return EXIT_FAILURE;
  case CMR_ERROR_MEMORY:
    fputs("Memory error.\n", stderr);
    return EXIT_FAILURE;
  default:
    return EXIT_SUCCESS;
  }
}

