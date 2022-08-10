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

typedef struct
{
  size_t row;
  size_t column;
  double value;
  char origin;
} DblNonzero;

int compareDblNonzeros(const void* pa, const void* pb)
{
  size_t aRow = ((const DblNonzero*)pa)->row;
  size_t bRow = ((const DblNonzero*)pb)->row;
  if (aRow < bRow)
    return -1;
  else if (aRow > bRow)
    return +1;

  int aColumn = ((const DblNonzero*)pa)->column;
  int bColumn = ((const DblNonzero*)pb)->column;
  if (aColumn < bColumn)
    return -1;
  else if (aColumn > bColumn)
    return +1;

  int aOrigin = ((const DblNonzero*)pa)->origin;
  int bOrigin = ((const DblNonzero*)pb)->origin;
  
  return bOrigin - aOrigin;
}

CMR_ERROR perturbMatrix(
  const char* inputMatrixFileName,  /**< Input matrix file name. */
  FileFormat inputFormat,           /**< Input file format. */
  const char* outputMatrixFileName, /**< Output matrix file name. */
  FileFormat outputFormat,          /**< Output file format. */
  size_t makeZero,                  /**< Number of nonzeros to be turned into zeros. */
  size_t makeOne,                   /**< Number of zeros to be turned into 1s. */
  size_t makeMinusOne,              /**< Number of zeros to be turned into -1s. */
  size_t flipBinary,                /**< Number of entries to be flipped over {0,1}. */
  size_t flipTernary                /**< Number of entries to be flipped over {-1,0,1}. */
)
{
  /* Read the matrix. */
  FILE* inputMatrixFile = strcmp(inputMatrixFileName, "-") ? fopen(inputMatrixFileName, "r") : stdin;
  if (!inputMatrixFile)
    return CMR_ERROR_INPUT;
  
  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_DBLMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRdblmatCreateFromSparseStream(cmr, inputMatrixFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRdblmatCreateFromDenseStream(cmr, inputMatrixFile, &matrix) );
  else
    return CMR_ERROR_INPUT;
  if (inputMatrixFile != stdin)
    fclose(inputMatrixFile);

  /* Make zeros. */
  size_t numNonzeros = matrix->numNonzeros;
  if (makeZero >= numNonzeros)
  {
    for (size_t i = 0; i < numNonzeros; ++i)
      matrix->entryValues[i] = 0.0;
  }
  else
  {
    for (size_t i = 0; i < makeZero; ++i)
    {
      size_t k = randRange(0, matrix->numNonzeros);
      if (matrix->entryValues[k])
      {
      matrix->entryValues[k] = 0.0;
      --numNonzeros;
      }
    }
  }

  numNonzeros += flipBinary + flipTernary + makeOne + makeMinusOne;
  DblNonzero* nonzeros = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &nonzeros, numNonzeros) );
  size_t entry = 0;

  /* We now flip because this may also remove entries. */
  for (size_t i = 0; i < flipBinary; ++i)
  {
    size_t row = randRange(0, matrix->numRows);
    size_t column = randRange(0, matrix->numColumns);
    size_t j = SIZE_MAX;
    CMR_CALL( CMRdblmatFindEntry(matrix, row, column, &j) );
    
    if (j == SIZE_MAX)
    {
      nonzeros[entry].row = row;
      nonzeros[entry].column = column;
      nonzeros[entry].value = 1;
      nonzeros[entry].origin = 1;
      ++entry;
    }
    else
      matrix->entryValues[j] = 0.0;
  }
  for (size_t i = 0; i < flipTernary; ++i)
  {
    size_t row = randRange(0, matrix->numRows);
    size_t column = randRange(0, matrix->numColumns);
    size_t j = SIZE_MAX;
    CMR_CALL( CMRdblmatFindEntry(matrix, row, column, &j) );
    
    if (j == SIZE_MAX)
    {
      nonzeros[entry].row = row;
      nonzeros[entry].column = column;
      nonzeros[entry].value = 2.0*randRange(0, 2) - 1.0;
      nonzeros[entry].origin = 2;
      ++entry;
    }
    else
      matrix->entryValues[j] = 0.0;
  }

  /* Now we copy the remaining entries. */
  for (size_t row = 0; row < matrix->numRows; ++row)
  {
    for (size_t i = matrix->rowSlice[row]; i < matrix->rowSlice[row+1]; ++i)
    {
      if (matrix->entryValues[i] != 0.0)
      {
        nonzeros[entry].row = row;
        nonzeros[entry].column = matrix->entryColumns[i];
        nonzeros[entry].value = matrix->entryValues[i];
        nonzeros[entry].origin = 0;
        ++entry;
      }
    }
  }

  /* We now create more. */
  for (size_t i = 0; i < makeOne; ++i)
  {
    nonzeros[entry].row = randRange(0, matrix->numRows);
    nonzeros[entry].column = randRange(0, matrix->numColumns);
    nonzeros[entry].value = 1;
    nonzeros[entry].origin = 3;
    ++entry;
  }
  for (size_t i = 0; i < makeMinusOne; ++i)
  {
    nonzeros[entry].row = randRange(0, matrix->numRows);
    nonzeros[entry].column = randRange(0, matrix->numColumns);
    nonzeros[entry].value = -1;
    nonzeros[entry].origin = 4;
    ++entry;
  }

  /* We sort all nonzeros by row and then by column. */
  qsort(nonzeros, entry, sizeof(DblNonzero), compareDblNonzeros);

  CMR_DBLMAT* result = NULL;
  CMR_CALL( CMRdblmatCreate(cmr, &result, matrix->numRows, matrix->numColumns, entry) );
  size_t previousRow = SIZE_MAX;
  size_t previousColumn = SIZE_MAX;
  numNonzeros = entry;
  entry = 0;
  for (size_t i = 0; i < numNonzeros; ++i)
  {
    size_t row = nonzeros[i].row;
    size_t column = nonzeros[i].column;
    if (row == previousRow && column == previousColumn)
      continue;

    while (previousRow < row || previousRow == SIZE_MAX)
    {
      ++previousRow;
      result->rowSlice[previousRow] = entry;
    }

    result->entryColumns[entry] = column;
    result->entryValues[entry] = nonzeros[i].value;
    ++entry;
    previousColumn = column;
  }
  while (previousRow < result->numRows || previousRow == SIZE_MAX)
  {
    ++previousRow;
    result->rowSlice[previousRow] = entry;
  }
  result->numNonzeros = entry;

  CMR_CALL( CMRfreeBlockArray(cmr, &nonzeros) );

  /* Output the matrix. */
  FILE* outputMatrixFile = strcmp(outputMatrixFileName, "-") ? fopen(outputMatrixFileName, "w") : stdout;

  if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRdblmatPrintDense(cmr, result, outputMatrixFile, '0', false) );
  else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRdblmatPrintSparse(cmr, result, outputMatrixFile) );
  if (outputMatrixFile != stdout)
    fclose(outputMatrixFile);

  /* Cleanup. */
  CMR_CALL( CMRdblmatFree(cmr, &result) );
  CMR_CALL( CMRdblmatFree(cmr, &matrix) );

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

int printUsage(const char* program)
{
  fprintf(stderr, "Usage: %s IN-MAT OUT-MAT [OPTION]...\n\n", program);
  fputs("Copies the matrix from file IN-MAT to the file OUT-MAT after applying perturbations.\n\n", stderr);
  fputs("Options:\n", stderr);
  fputs("  -i FORMAT Format of file IN-MAT, among `dense' and `sparse'; default: dense.\n", stderr);
  fputs("  -o FORMAT Format of file IN-MAT, among `dense' and `sparse'; default: same as format of IN-MAT.\n", stderr);
  fputs("  -0 NUM    Turn NUM randomly chosen nonzero entries to 0s.\n", stderr);
  fputs("  -1 NUM    Turn NUM randomly chosen zero entries into 1s.\n", stderr);
  fputs("  --1 NUM   Turn NUM randomly chosen zero entries into -1s.\n", stderr);
  fputs("  -b NUM    Flip NUM randomly chosen entries over the binary field.\n", stderr);
  fputs("  -t NUM    Flip NUM randomly chosen entries over the ternary field.\n\n", stderr);
  fputs("If IN-MAT is `-' then the input matrix is read from stdin.\n", stderr);
  fputs("If OUT-MAT is `-' then the output matrix is written to stdout.\n", stderr);
  
  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  struct timeval curTime;
  gettimeofday(&curTime, NULL);
  srand(curTime.tv_usec);

  char* inputMatrixFileName = NULL;
  char* outputMatrixFileName = NULL;
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  FileFormat outputFormat = FILEFORMAT_UNDEFINED;
  size_t makeZero = 0;
  size_t makeOne = 0;
  size_t makeMinusOne = 0;
  size_t flipBinary = 0;
  size_t flipTernary = 0;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-i") && (a+1 < argc))
    {
      if (!strcmp(argv[a+1], "dense"))
        inputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        inputFormat = FILEFORMAT_MATRIX_SPARSE;
      else
      {
        printf("Error: unknown input format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
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
    else if (!strcmp(argv[a], "-0") && (a+1 < argc))
    {
      char* p = NULL;
      makeZero = strtoul(argv[a+1], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-1") && (a+1 < argc))
    {
      char* p = NULL;
      makeOne = strtoul(argv[a+1], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      ++a;
    }
    else if (!strcmp(argv[a], "--1") && (a+1 < argc))
    {
      char* p = NULL;
      makeMinusOne = strtoul(argv[a+1], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-b") && (a+1 < argc))
    {
      char* p = NULL;
      flipBinary = strtoul(argv[a+1], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-t") && (a+1 < argc))
    {
      char* p = NULL;
      flipTernary = strtoul(argv[a+1], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      ++a;
    }
    else if (!inputMatrixFileName)
      inputMatrixFileName = argv[a];
    else if (!outputMatrixFileName)
      outputMatrixFileName = argv[a];
    else
    {
      fprintf(stderr, "Error: Three matrix files specified: %s, %s and %s\n\n", inputMatrixFileName,
        outputMatrixFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (!inputMatrixFileName)
  {
    fputs("Error: No input matrix specified.\n\n", stderr);
    return printUsage(argv[0]);
  }
  if (!outputMatrixFileName)
  {
    fputs("Error: No output matrix specified.\n\n", stderr);
    return printUsage(argv[0]);
  }
  if (outputFormat == FILEFORMAT_UNDEFINED)
    outputFormat = inputFormat;

  CMR_ERROR error = perturbMatrix(inputMatrixFileName, inputFormat, outputMatrixFileName, outputFormat, makeZero,
    makeOne, makeMinusOne, flipBinary, flipTernary);

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
