#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/time.h>

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/matroid.h>

// #define RANDOM_SEED 0 /* Define 0 for logging the seed to stderr, and something else for the actual seed. */
// #define RANDOM_SEED 133974

typedef enum
{
  TASK_COPY = 0,
  TASK_SUPPORT = 1,
  TASK_SIGNED_SUPPORT = 2,
} Task;

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,     /**< Whether the file format of input/output was defined by the user. */
  FILEFORMAT_MATRIX_DENSE = 1,  /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2  /**< Sparse matrix format. */
} FileFormat;

static
CMR_ERROR writeToFileDbl(CMR* cmr, CMR_DBLMAT* matrix, FileFormat outputFormat, const char* outputMatrixFileName,
  bool transpose)
{
  assert(matrix);

  CMR_DBLMAT* output = NULL;
  if (transpose)
    CMR_CALL( CMRdblmatTranspose(cmr, matrix, &output) );
  else
    output = matrix;

  bool outputMatrixToFile = strcmp(outputMatrixFileName, "-");
  FILE* outputMatrixFile = outputMatrixToFile ? fopen(outputMatrixFileName, "w") : stdout;

  CMR_ERROR error = CMR_OKAY;
  if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRdblmatPrintSparse(cmr, output, outputMatrixFile) );
  else if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRdblmatPrintDense(cmr, output, outputMatrixFile, '0', false) );
  else
    error = CMR_ERROR_INPUT;

  if (outputMatrixToFile)
    fclose(outputMatrixFile);

  if (transpose)
    CMR_CALL( CMRdblmatFree(cmr, &output) );

  return error;
}

static
CMR_ERROR writeToFileInt(CMR* cmr, CMR_INTMAT* matrix, FileFormat outputFormat, const char* outputMatrixFileName,
  bool transpose)
{
  assert(matrix);

  CMR_INTMAT* output = NULL;
  if (transpose)
    CMR_CALL( CMRintmatTranspose(cmr, matrix, &output) );
  else
    output = matrix;

  bool outputMatrixToFile = strcmp(outputMatrixFileName, "-");
  FILE* outputMatrixFile = outputMatrixToFile ? fopen(outputMatrixFileName, "w") : stdout;

  CMR_ERROR error = CMR_OKAY;
  if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRintmatPrintSparse(cmr, output, outputMatrixFile) );
  else if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRintmatPrintDense(cmr, output, outputMatrixFile, '0', false) );
  else
    error = CMR_ERROR_INPUT;

  if (outputMatrixToFile)
    fclose(outputMatrixFile);

  if (transpose)
    CMR_CALL( CMRintmatFree(cmr, &output) );

  return error;
}

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

static
void initRandomPermutation(size_t* permutation, size_t length)
{
  for (size_t i = 0; i < length; ++i)
    permutation[i] = i;
  for (size_t i = 0; i < length; ++i)
  {
    size_t j = randRange(i, length);
    size_t temp = permutation[i];
    permutation[i] = permutation[j];
    permutation[j] = temp;
  }
}

static
CMR_ERROR writeToFileChr(CMR* cmr, CMR_CHRMAT* matrix, FileFormat outputFormat, const char* outputMatrixFileName,
  bool transpose)
{
  assert(matrix);

  CMR_CHRMAT* output = NULL;
  if (transpose)
    CMR_CALL( CMRchrmatTranspose(cmr, matrix, &output) );
  else
    output = matrix;

  bool outputMatrixToFile = strcmp(outputMatrixFileName, "-");
  FILE* outputMatrixFile = outputMatrixToFile ? fopen(outputMatrixFileName, "w") : stdout;

  CMR_ERROR error = CMR_OKAY;
  if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatPrintSparse(cmr, output, outputMatrixFile) );
  else if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatPrintDense(cmr, output, outputMatrixFile, '0', false) );
  else
    error = CMR_ERROR_INPUT;

  if (outputMatrixToFile)
    fclose(outputMatrixFile);

  if (transpose)
    CMR_CALL( CMRchrmatFree(cmr, &output) );

  return error;
}

static
CMR_ERROR runDbl(
  const char* inputMatrixFileName,
  FileFormat inputFormat,
  const char* inputSubmatrixFileName,
  FileFormat outputFormat,
  const char* outputMatrixFileName,
  Task task,
  bool randomPermute,
  int randomPivotsType,
  size_t randomPivotsCount,
  bool transpose
)
{
  FILE* inputMatrixFile = strcmp(inputMatrixFileName, "-") ? fopen(inputMatrixFileName, "r") : stdin;
  if (!inputMatrixFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_DBLMAT* dblmatrix = NULL;
  CMR_CHRMAT* chrmatrix = NULL;
  CMR_ERROR error = CMR_OKAY;
  if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    error = CMRdblmatCreateFromSparseStream(cmr, inputMatrixFile, &dblmatrix);
  else if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    error = CMRdblmatCreateFromDenseStream(cmr, inputMatrixFile, &dblmatrix);
  else
    assert(false);

  if (inputMatrixFile != stdin)
    fclose(inputMatrixFile);

  if (error != CMR_OKAY)
  {
    CMR_CALL( CMRfreeEnvironment(&cmr) );
    return CMR_OKAY;
  }

  /* Slice a submatrix if requested. */
  if (inputSubmatrixFileName)
  {
    CMR_DBLMAT* explicitSubmatrix = NULL;
    bool inputSubmatrixFromFile = strcmp(inputSubmatrixFileName, "-");
    FILE* inputSubmatrixFile = inputSubmatrixFromFile ? fopen(inputSubmatrixFileName, "r") : stdin;
    CMR_SUBMAT* submatrix = NULL;
    size_t numRows, numColumns;
    CMR_CALL( CMRsubmatReadFromStream(cmr, &submatrix, &numRows, &numColumns, inputSubmatrixFile) );

    if (numRows != dblmatrix->numRows)
    {
      fprintf(stderr, "Error: Number of rows from IN-MAT (%zu) and IN-SUB (%zu) do not match.\n", dblmatrix->numRows,
        numRows);
      CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
      if (inputSubmatrixFileName)
        CMR_CALL( CMRdblmatFree(cmr, &explicitSubmatrix) );
      CMR_CALL( CMRdblmatFree(cmr, &dblmatrix) );
      CMR_CALL( CMRfreeEnvironment(&cmr) );
      return CMR_ERROR_INPUT;
    }

    if (numColumns != dblmatrix->numColumns)
    {
      fprintf(stderr, "Error: Number of columns from IN-MAT (%zu) and IN-SUB (%zu) do not match.\n",
        dblmatrix->numColumns, numColumns);
      CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
      if (inputSubmatrixFileName)
        CMR_CALL( CMRdblmatFree(cmr, &explicitSubmatrix) );
      CMR_CALL( CMRdblmatFree(cmr, &dblmatrix) );
      CMR_CALL( CMRfreeEnvironment(&cmr) );
      return CMR_ERROR_INPUT;
    }

    explicitSubmatrix = NULL;
    CMR_CALL( CMRdblmatSlice(cmr, dblmatrix, submatrix, &explicitSubmatrix) );
    CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
    CMR_CALL( CMRdblmatFree(cmr, &dblmatrix) );
    dblmatrix = explicitSubmatrix;
  }

  /* Permute rows if requested. */
  if (randomPermute)
  {
    size_t* permutedRowsToRows = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &permutedRowsToRows, dblmatrix->numRows) );
    initRandomPermutation(permutedRowsToRows, dblmatrix->numRows);

    size_t* columnsToPermutedColumns = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &columnsToPermutedColumns, dblmatrix->numColumns) );
    initRandomPermutation(columnsToPermutedColumns, dblmatrix->numColumns);

    CMR_DBLMAT* permuted = NULL;
    CMR_CALL( CMRdblmatCreate(cmr, &permuted, dblmatrix->numRows, dblmatrix->numColumns, dblmatrix->numNonzeros) );
    permuted->numNonzeros = 0;
    for (size_t row = 0; row < permuted->numRows; ++row)
    {
      permuted->rowSlice[row] = permuted->numNonzeros;
      size_t originalRow = permutedRowsToRows[row];
      size_t first = dblmatrix->rowSlice[originalRow];
      size_t beyond = dblmatrix->rowSlice[originalRow + 1];
      for (size_t i = first; i < beyond; ++i)
      {
        size_t originalColumn = dblmatrix->entryColumns[i];
        permuted->entryValues[permuted->numNonzeros] = dblmatrix->entryValues[i];
        permuted->entryColumns[permuted->numNonzeros++] = columnsToPermutedColumns[originalColumn];
      }
    }
    permuted->rowSlice[permuted->numRows] = permuted->numNonzeros;
    CMR_CALL( CMRdblmatSortNonzeros(cmr, permuted) );

    CMR_CALL( CMRfreeBlockArray(cmr, &columnsToPermutedColumns) );
    CMR_CALL( CMRfreeBlockArray(cmr, &permutedRowsToRows) );
    CMR_CALL( CMRdblmatFree(cmr, &dblmatrix) );
    dblmatrix = permuted;
  }

  /* Transpose the matrix if requested. */
  if (transpose)
  {
    CMR_DBLMAT* transpose = NULL;
    CMR_CALL( CMRdblmatTranspose(cmr, dblmatrix, &transpose) );
    CMR_CALL( CMRdblmatFree(cmr, &dblmatrix) );
    dblmatrix = transpose;
  }

  if (task == TASK_SUPPORT)
  {
    CMR_CALL( CMRdblmatSupport(cmr, dblmatrix, 1.0e-9, &chrmatrix) );
    CMR_CALL( CMRdblmatFree(cmr, &dblmatrix) );
  }
  else if (task == TASK_SIGNED_SUPPORT)
  {
    CMR_CALL( CMRdblmatSignedSupport(cmr, dblmatrix, 1.0e-9, &chrmatrix) );
    CMR_CALL( CMRdblmatFree(cmr, &dblmatrix) );
  }
  else
  {
    assert( task == TASK_COPY );
  }

  /* Pivot randomly if requested. */
  if (randomPivotsType >= 0 && randomPivotsType > 0)
  {
    if (!chrmatrix)
    {
      CMR_CALL( CMRdblmatToChr(cmr, dblmatrix, 1.0e-9, &chrmatrix) );
      CMR_CALL( CMRdblmatFree(cmr, &dblmatrix) );
      assert(dblmatrix == NULL);
    }

    if (chrmatrix->numNonzeros == 0)
    {
      fprintf(stderr, "WARNING: skipping pivots for zero-matrix.\n");
      randomPivotsCount = 0;
    }

    for (size_t p = 0; p < randomPivotsCount; ++p)
    {
      CMR_CHRMAT* result = NULL;
      size_t entry = randRange(0, chrmatrix->numNonzeros);
      size_t column = chrmatrix->entryColumns[entry];

      /* Search for the row of the entry. */
      size_t rowLower = 0;
      size_t rowUpper = chrmatrix->numRows;
      while (rowLower + 1 < rowUpper)
      {
        size_t rowMiddle = (rowLower + rowUpper) / 2;
        size_t middleEntry = chrmatrix->rowSlice[rowMiddle];
        if (entry >= middleEntry)
          rowLower = rowMiddle;
        else
          rowUpper = rowMiddle;
      }
      assert(entry >= chrmatrix->rowSlice[rowLower]);
      assert(entry < chrmatrix->rowSlice[rowLower + 1]);

      /* Carry out a single pivot. */
      if (randomPivotsType == 2)
        CMR_CALL( CMRchrmatBinaryPivot(cmr, chrmatrix, rowLower, column, &result) );
      else if (randomPivotsType == 3)
        CMR_CALL( CMRchrmatTernaryPivot(cmr, chrmatrix, rowLower, column, &result) );
      else
      {
        assert(false);
      }
      CMR_CALL( CMRchrmatFree(cmr, &chrmatrix) );
      chrmatrix = result;
    }
  }

  /* Finally, write and free. */
  if (dblmatrix)
  {
    CMR_CALL( writeToFileDbl(cmr, dblmatrix, outputFormat, outputMatrixFileName, false) );
    CMR_CALL( CMRdblmatFree(cmr, &dblmatrix) );
  }
  if (chrmatrix)
  {
    CMR_CALL( writeToFileChr(cmr, chrmatrix, outputFormat, outputMatrixFileName, false) );
    CMR_CALL( CMRchrmatFree(cmr, &chrmatrix) );
  }

  return CMR_OKAY;
}

static
CMR_ERROR runInt(
  const char* inputMatrixFileName,
  FileFormat inputFormat,
  const char* inputSubmatrixFileName,
  FileFormat outputFormat,
  const char* outputMatrixFileName,
  Task task,
  bool randomPermute,
  int randomPivotsType,
  size_t randomPivotsCount,
  bool transpose
)
{
  FILE* inputMatrixFile = strcmp(inputMatrixFileName, "-") ? fopen(inputMatrixFileName, "r") : stdin;
  if (!inputMatrixFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_INTMAT* intmatrix = NULL;
  CMR_CHRMAT* chrmatrix = NULL;
  CMR_ERROR error = CMR_OKAY;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
  {
    error = CMRintmatCreateFromDenseStream(cmr, inputMatrixFile, &intmatrix);
    if (error == CMR_ERROR_INPUT)
      fprintf(stderr, "Error when reading dense matrix from <%s>: %s\n", inputMatrixFileName, CMRgetErrorMessage(cmr));
  }
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
  {
    error = CMRintmatCreateFromSparseStream(cmr, inputMatrixFile, &intmatrix);
    if (error == CMR_ERROR_INPUT)
      fprintf(stderr, "Error when reading sparse matrix from <%s>: %s\n", inputMatrixFileName, CMRgetErrorMessage(cmr));
  }
  else
    assert(false);

  if (inputMatrixFile != stdin)
    fclose(inputMatrixFile);

  if (error != CMR_OKAY)
  {
    CMR_CALL( CMRfreeEnvironment(&cmr) );
    return CMR_OKAY;
  }

  /* Slice a submatrix if requested. */
  if (inputSubmatrixFileName)
  {
    CMR_INTMAT* explicitSubmatrix;
    bool inputSubmatrixFromFile = strcmp(inputSubmatrixFileName, "-");
    FILE* inputSubmatrixFile = inputSubmatrixFromFile ? fopen(inputSubmatrixFileName, "r") : stdin;
    CMR_SUBMAT* submatrix = NULL;
    size_t numRows, numColumns;
    CMR_CALL( CMRsubmatReadFromStream(cmr, &submatrix, &numRows, &numColumns, inputSubmatrixFile) );

    if (numRows != intmatrix->numRows)
    {
      fprintf(stderr, "Error: Number of rows from IN-MAT (%zu) and IN-SUB (%zu) do not match.\n", intmatrix->numRows,
        numRows);
      CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
      if (inputSubmatrixFileName)
        CMR_CALL( CMRintmatFree(cmr, &explicitSubmatrix) );
      CMR_CALL( CMRintmatFree(cmr, &intmatrix) );
      CMR_CALL( CMRfreeEnvironment(&cmr) );
      return CMR_ERROR_INPUT;
    }

    if (numColumns != intmatrix->numColumns)
    {
      fprintf(stderr, "Error: Number of columns from IN-MAT (%zu) and IN-SUB (%zu) do not match.\n", intmatrix->numColumns,
        numColumns);
      CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
      if (inputSubmatrixFileName)
        CMR_CALL( CMRintmatFree(cmr, &explicitSubmatrix) );
      CMR_CALL( CMRintmatFree(cmr, &intmatrix) );
      CMR_CALL( CMRfreeEnvironment(&cmr) );
      return CMR_ERROR_INPUT;
    }

    explicitSubmatrix = NULL;
    CMR_CALL( CMRintmatSlice(cmr, intmatrix, submatrix, &explicitSubmatrix) );
    CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
    CMR_CALL( CMRintmatFree(cmr, &intmatrix) );
    intmatrix = explicitSubmatrix;
  }

  /* Permute rows if requested. */
  if (randomPermute)
  {
    size_t* permutedRowsToRows = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &permutedRowsToRows, intmatrix->numRows) );
    initRandomPermutation(permutedRowsToRows, intmatrix->numRows);

    size_t* columnsToPermutedColumns = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &columnsToPermutedColumns, intmatrix->numColumns) );
    initRandomPermutation(columnsToPermutedColumns, intmatrix->numColumns);

    CMR_INTMAT* permuted = NULL;
    CMR_CALL( CMRintmatCreate(cmr, &permuted, intmatrix->numRows, intmatrix->numColumns, intmatrix->numNonzeros) );
    permuted->numNonzeros = 0;
    for (size_t row = 0; row < permuted->numRows; ++row)
    {
      permuted->rowSlice[row] = permuted->numNonzeros;
      size_t originalRow = permutedRowsToRows[row];
      size_t first = intmatrix->rowSlice[originalRow];
      size_t beyond = intmatrix->rowSlice[originalRow + 1];
      for (size_t i = first; i < beyond; ++i)
      {
        size_t originalColumn = intmatrix->entryColumns[i];
        permuted->entryValues[permuted->numNonzeros] = intmatrix->entryValues[i];
        permuted->entryColumns[permuted->numNonzeros++] = columnsToPermutedColumns[originalColumn];
      }
    }
    permuted->rowSlice[permuted->numRows] = permuted->numNonzeros;
    CMR_CALL( CMRintmatSortNonzeros(cmr, permuted) );

    CMR_CALL( CMRfreeBlockArray(cmr, &columnsToPermutedColumns) );
    CMR_CALL( CMRfreeBlockArray(cmr, &permutedRowsToRows) );
    CMR_CALL( CMRintmatFree(cmr, &intmatrix) );
    intmatrix = permuted;
  }

  /* Transpose the matrix if requested. */
  if (transpose)
  {
    CMR_INTMAT* transpose = NULL;
    CMR_CALL( CMRintmatTranspose(cmr, intmatrix, &transpose) );
    CMR_CALL( CMRintmatFree(cmr, &intmatrix) );
    intmatrix = transpose;
  }

  if (task == TASK_SUPPORT)
  {
    CMR_CALL( CMRintmatSupport(cmr, intmatrix, &chrmatrix) );
    CMR_CALL( CMRintmatFree(cmr, &intmatrix) );
    assert(intmatrix == NULL);
  }
  else if (task == TASK_SIGNED_SUPPORT)
  {
    CMR_CALL( CMRintmatSignedSupport(cmr, intmatrix, &chrmatrix) );
    CMR_CALL( CMRintmatFree(cmr, &intmatrix) );
    assert(intmatrix == NULL);
  }
  else
  {
    assert( task == TASK_COPY );
  }

  /* Pivot randomly if requested. */
  if (randomPivotsType >= 0 && randomPivotsType > 0)
  {
    if (!chrmatrix)
    {
      CMR_CALL( CMRintmatToChr(cmr, intmatrix, &chrmatrix) );
      CMR_CALL( CMRintmatFree(cmr, &intmatrix) );
      assert(intmatrix == NULL);
    }

    if (chrmatrix->numNonzeros == 0)
    {
      fprintf(stderr, "WARNING: skipping pivots for zero-matrix.\n");
      randomPivotsCount = 0;
    }

    for (size_t p = 0; p < randomPivotsCount; ++p)
    {
      CMR_CHRMAT* result = NULL;
      size_t entry = randRange(0, chrmatrix->numNonzeros);
      size_t column = chrmatrix->entryColumns[entry];

      /* Search for the row of the entry. */
      size_t rowLower = 0;
      size_t rowUpper = chrmatrix->numRows;
      while (rowLower + 1 < rowUpper)
      {
        size_t rowMiddle = (rowLower + rowUpper) / 2;
        size_t middleEntry = chrmatrix->rowSlice[rowMiddle];
        if (entry >= middleEntry)
          rowLower = rowMiddle;
        else
          rowUpper = rowMiddle;
      }
      assert(entry >= chrmatrix->rowSlice[rowLower]);
      assert(entry < chrmatrix->rowSlice[rowLower + 1]);

      /* Carry out a single pivot. */
      if (randomPivotsType == 2)
        CMR_CALL( CMRchrmatBinaryPivot(cmr, chrmatrix, rowLower, column, &result) );
      else if (randomPivotsType == 3)
        CMR_CALL( CMRchrmatTernaryPivot(cmr, chrmatrix, rowLower, column, &result) );
      else
      {
        assert(false);
      }
      CMR_CALL( CMRchrmatFree(cmr, &chrmatrix) );
      chrmatrix = result;
    }
  }

  /* Finally, write and free. */
  if (intmatrix)
  {
    CMR_CALL( writeToFileInt(cmr, intmatrix, outputFormat, outputMatrixFileName, false) );
    CMR_CALL( CMRintmatFree(cmr, &intmatrix) );
  }
  if (chrmatrix)
  {
    CMR_CALL( writeToFileChr(cmr, chrmatrix, outputFormat, outputMatrixFileName, false) );
    CMR_CALL( CMRchrmatFree(cmr, &chrmatrix) );
  }

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

int printUsage(const char* program)
{
  fputs("Usage:\n", stderr);
  fprintf(stderr, "%s IN-MAT OUT-MAT [OPTION]...\n\n", program);
  fputs("  copies the matrix from file IN-MAT to file OUT-MAT, potentially applying certain operations.\n\n", stderr);
  fputs("Options:\n", stderr);
  fputs("  -i FORMAT Format of file IN-MAT, among `dense' and `sparse'; default: dense.\n", stderr);
  fputs("  -o FORMAT Format of file OUT-MAT, among `dense' and `sparse'; default: same format as of IN-MAT.\n", stderr);
  fputs("  -S IN-SUB Consider the submatrix of IN-MAT specified in file IN-SUB instead of IN-MAT itself;"
    " can be combined with other operations.\n", stderr);
  fputs("  -t        Transpose the matrix; can be combined with other operations.\n", stderr);
  fputs("  -c        Compute the support matrix instead of copying.\n", stderr);
  fputs("  -C        Compute the signed support matrix instead of copying.\n", stderr);
  fputs("  -r        Randomize the output matrix by randomly permuting rows/columns.\n", stderr);
  fputs("  -R2 NUM   Randomize the output matrix by performing NUM random binary pivots.\n", stderr);
  fputs("  -R3 NUM   Randomize the output matrix by performing NUM random ternary pivots.\n", stderr);
  fputs("  -d        Use double arithmetic instead of integers.\n\n", stderr);
  fputs("If IN-MAT is `-' then the input matrix is read from stdin.\n", stderr);
  fputs("If OUT-MAT is `-' then the output matrix is written to stdout.\n", stderr);

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  Task task = TASK_COPY;
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  FileFormat outputFormat = FILEFORMAT_UNDEFINED;
  bool randomPermute = false;
  int randomPivotsType = -1;
  size_t randomPivotsCount = 0;
  bool transpose = false;
  bool doubleArithmetic = false;
  char* inputMatrixFileName = NULL;
  char* inputSubmatrixFileName = NULL;
  char* outputMatrixFileName = NULL;
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
        fprintf(stderr, "Error: Unknown input format <%s>.\n\n", argv[a+1]);
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
        fprintf(stderr, "Error: Unknown output format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-S") && a+1 < argc)
      inputSubmatrixFileName = argv[++a];
    else if (!strcmp(argv[a], "-c"))
      task = TASK_SUPPORT;
    else if (!strcmp(argv[a], "-C"))
      task = TASK_SIGNED_SUPPORT;
    else if (!strcmp(argv[a], "-r"))
      randomPermute = true;
    else if (!strcmp(argv[a], "-R2") && a + 1 < argc)
    {
      randomPivotsType = 2;
      char* p;
      randomPivotsCount = strtoull(argv[a+1], &p, 10);
      if (*p != '\0' || randomPivotsCount == 0)
      {
        fprintf(stderr, "Error: invalid number of binary pivots <%s>.\n", argv[a+1]);
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a++;
    }
    else if (!strcmp(argv[a], "-R3") && a + 1 < argc)
    {
      randomPivotsType = 3;
      char* p;
      randomPivotsCount = strtoull(argv[a+1], &p, 10);
      if (*p != '\0' || randomPivotsCount == 0)
      {
        fprintf(stderr, "Error: invalid number of ternary pivots <%s>.\n", argv[a+1]);
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a++;
    }
    else if (!strcmp(argv[a], "-t"))
      transpose = true;
    else if (!strcmp(argv[a], "-d"))
      doubleArithmetic = true;
    else if (!inputMatrixFileName)
      inputMatrixFileName = argv[a];
    else if (!outputMatrixFileName)
      outputMatrixFileName = argv[a];
    else
    {
      fprintf(stderr, "Error: Three matrix files <%s>, <%s> and <%s> specified.\n\n", inputMatrixFileName,
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

  if (randomPermute || (randomPivotsType >= 0 && randomPivotsCount > 0))
  {
    struct timeval curTime;
    gettimeofday(&curTime, NULL);
    unsigned int seed = curTime.tv_usec;
#ifdef RANDOM_SEED
#if RANDOM_SEED > 0
    seed = RANDOM_SEED;
#endif
    fprintf(stderr, "Random seed: %u\n", seed);
#endif /* RANDOM_SEED */
    srand(seed);
  }

  CMR_ERROR error;
  if (doubleArithmetic)
  {
    error = runDbl(inputMatrixFileName, inputFormat, inputSubmatrixFileName, outputFormat, outputMatrixFileName, task,
      randomPermute, randomPivotsType, randomPivotsCount, transpose);
  }
  else
  {
    error = runInt(inputMatrixFileName, inputFormat, inputSubmatrixFileName, outputFormat, outputMatrixFileName, task,
      randomPermute, randomPivotsType, randomPivotsCount, transpose);
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
