#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
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
  printf("Usage: %s [OPTIONS] ROWS COLS p\n\n", program);
  puts("Creates a random ROWS-by-COLS 0/1 matrix in which each entry is 1 with probability p.\n");
  puts("Options:\n");
  puts("  -o FORMAT  Format of output FILE; default: `dense'.");
  puts("Formats for matrices: dense, sparse");
  return EXIT_FAILURE;
}


CMR_ERROR genMatrixRandom(
  size_t numRows,         /**< Number of rows of base matrix. */
  size_t numColumns,      /**< Number of columns of base matrix. */
  double probability1,    /**< Probability for a 1-entry. */
  FileFormat outputFormat /**< Output file format. */
)
{
  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );
  
  size_t estimatedNumNonzeros = 1.1 * numRows * numColumns * probability1 + 1024;

  CMR_CHRMAT* matrix = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &matrix, numRows, numColumns, estimatedNumNonzeros) );
  size_t entry = 0;
  for (size_t row = 0; row < numRows; ++row)
  {
    matrix->rowSlice[row] = entry;
    for (size_t column = 0; column < numColumns; ++column)
    {
      bool isNonzero = (rand() * 1.0 / RAND_MAX) < probability1;
      if (isNonzero)
      {
        if (entry == matrix->numNonzeros)
        {
          CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryColumns, 2*matrix->numNonzeros) );
          CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryValues, 2*matrix->numNonzeros) );
          matrix->numNonzeros *= 2;
        }
        matrix->entryColumns[entry] = column;
        matrix->entryValues[entry] = 1;
        ++entry;
      }
    }
  }
  matrix->rowSlice[numRows] = entry;
  matrix->numNonzeros = entry;

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
  size_t numRows = SIZE_MAX;
  size_t numColumns = SIZE_MAX;
  double probability1 = 0.5;
  bool readProbability1 = false;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-o") && (a+1 < argc))
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
    else if (numRows == SIZE_MAX)
    {
      char* p = NULL;
      numRows = strtoull(argv[a], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
    }
    else if (numColumns == SIZE_MAX)
    {
      char* p = NULL;
      numColumns = strtoull(argv[a], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
    }
    else if (!readProbability1)
    {
      char* p = NULL;
      probability1 = strtod(argv[a], &p);
      readProbability1 = true;
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
    }
    else
    {
      printf("Error: more than two size indicators specified: %ld %ld %s\n\n", numRows, numColumns, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (numRows == SIZE_MAX)
  {
    puts("Error: no size indicator specified.\n");
    return printUsage(argv[0]);
  }
  else if (numColumns == SIZE_MAX)
  {
    puts("Error: only one size indicator specified.\n");
    return printUsage(argv[0]);
  }
  else if (numRows <= 0 || numColumns <= 0)
  {
    puts("Error: matrix must have at least 1 row and 1 column.\n");
    return printUsage(argv[0]);
  }
  else if (!readProbability1)
  {
    puts("Error: no probability specified.\n");
    return printUsage(argv[0]);
  }
  else if (probability1 < 0.0 || readProbability1 > 1.0)
  {
    puts("Error: probability must be in [0,1].\n");
    return printUsage(argv[0]);
  }
  if (outputFormat == FILEFORMAT_UNDEFINED)
    outputFormat = FILEFORMAT_MATRIX_DENSE;

  CMR_ERROR error = genMatrixRandom(numRows, numColumns, probability1, outputFormat);
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
