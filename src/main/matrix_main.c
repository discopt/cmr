#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <cmr/matrix.h>

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

CMR_ERROR runDbl(
  const char* inputMatrixFileName,
  FileFormat inputFormat,
  const char* inputSubmatrixFileName,
  FileFormat outputFormat,
  const char* outputMatrixFileName,
  Task task,
  bool transpose
)
{
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

  CMR_DBLMAT* explicitSubmatrix;
  if (inputSubmatrixFileName)
  {
    bool inputSubmatrixFromFile = strcmp(inputSubmatrixFileName, "-");
    FILE* inputSubmatrixFile = inputSubmatrixFromFile ? fopen(inputSubmatrixFileName, "r") : stdin;
    CMR_SUBMAT* submatrix = NULL;
    size_t numRows, numColumns;
    CMR_CALL( CMRsubmatReadFromStream(cmr, &submatrix, &numRows, &numColumns, inputSubmatrixFile) );

    if (numRows != matrix->numRows)
    {
      fprintf(stderr, "Error: Number of rows from IN-MAT (%zu) and IN-SUB (%zu) do not match.\n", matrix->numRows,
        numRows);
      CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
      if (inputSubmatrixFileName)
        CMR_CALL( CMRdblmatFree(cmr, &explicitSubmatrix) );
      CMR_CALL( CMRdblmatFree(cmr, &matrix) );
      CMR_CALL( CMRfreeEnvironment(&cmr) );
      return CMR_ERROR_INPUT;
    }

    if (numColumns != matrix->numColumns)
    {
      fprintf(stderr, "Error: Number of columns from IN-MAT (%zu) and IN-SUB (%zu) do not match.\n", matrix->numColumns,
        numColumns);
      CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
      if (inputSubmatrixFileName)
        CMR_CALL( CMRdblmatFree(cmr, &explicitSubmatrix) );
      CMR_CALL( CMRdblmatFree(cmr, &matrix) );
      CMR_CALL( CMRfreeEnvironment(&cmr) );
      return CMR_ERROR_INPUT;
    }

    explicitSubmatrix = NULL;
    CMR_CALL( CMRdblmatZoomSubmat(cmr, matrix, submatrix, &explicitSubmatrix) );
    CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
  }
  else
    explicitSubmatrix = matrix;

  if (task == TASK_SUPPORT)
  {
    CMR_CHRMAT* result = NULL;
    CMR_CALL( CMRdblmatSupport(cmr, explicitSubmatrix, 1.0e-9, &result) );
    CMR_CALL( writeToFileChr(cmr, result, outputFormat, outputMatrixFileName, transpose) );
    CMR_CALL( CMRchrmatFree(cmr, &result) );
  }
  else if (task == TASK_SIGNED_SUPPORT)
  {
    CMR_CHRMAT* result = NULL;
    CMR_CALL( CMRdblmatSignedSupport(cmr, explicitSubmatrix, 1.0e-9, &result) );
    CMR_CALL( writeToFileChr(cmr, result, outputFormat, outputMatrixFileName, transpose) );
    CMR_CALL( CMRchrmatFree(cmr, &result) );
  }
  else
    CMR_CALL( writeToFileDbl(cmr, explicitSubmatrix, outputFormat, outputMatrixFileName, transpose) );

  /* Cleanup. */

  if (inputSubmatrixFileName)
    CMR_CALL( CMRdblmatFree(cmr, &explicitSubmatrix) );
  CMR_CALL( CMRdblmatFree(cmr, &matrix) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

CMR_ERROR runInt(
  const char* inputMatrixFileName,
  FileFormat inputFormat,
  const char* inputSubmatrixFileName,
  FileFormat outputFormat,
  const char* outputMatrixFileName,
  Task task,
  bool transpose
)
{
  FILE* inputMatrixFile = strcmp(inputMatrixFileName, "-") ? fopen(inputMatrixFileName, "r") : stdin;
  if (!inputMatrixFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_INTMAT* matrix = NULL;
  CMR_ERROR error = CMR_OKAY;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
  {
    error = CMRintmatCreateFromDenseStream(cmr, inputMatrixFile, &matrix);
    if (error == CMR_ERROR_INPUT)
      fprintf(stderr, "Error when reading dense matrix from <%s>: %s\n", inputMatrixFileName, CMRgetErrorMessage(cmr));
  }
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
  {
    error = CMRintmatCreateFromSparseStream(cmr, inputMatrixFile, &matrix);
    if (error == CMR_ERROR_INPUT)
      fprintf(stderr, "Error when reading sparse matrix from <%s>: %s\n", inputMatrixFileName, CMRgetErrorMessage(cmr));
  }
  else
    assert(false);
  if (inputMatrixFile != stdin)
    fclose(inputMatrixFile);

  if (error == CMR_OKAY)
  {
    CMR_INTMAT* explicitSubmatrix;
    if (inputSubmatrixFileName)
    {
      bool inputSubmatrixFromFile = strcmp(inputSubmatrixFileName, "-");
      FILE* inputSubmatrixFile = inputSubmatrixFromFile ? fopen(inputSubmatrixFileName, "r") : stdin;
      CMR_SUBMAT* submatrix = NULL;
      size_t numRows, numColumns;
      CMR_CALL( CMRsubmatReadFromStream(cmr, &submatrix, &numRows, &numColumns, inputSubmatrixFile) );

      if (numRows != matrix->numRows)
      {
        fprintf(stderr, "Error: Number of rows from IN-MAT (%zu) and IN-SUB (%zu) do not match.\n", matrix->numRows,
          numRows);
        CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
        if (inputSubmatrixFileName)
          CMR_CALL( CMRintmatFree(cmr, &explicitSubmatrix) );
        CMR_CALL( CMRintmatFree(cmr, &matrix) );
        CMR_CALL( CMRfreeEnvironment(&cmr) );
        return CMR_ERROR_INPUT;
      }

      if (numColumns != matrix->numColumns)
      {
        fprintf(stderr, "Error: Number of columns from IN-MAT (%zu) and IN-SUB (%zu) do not match.\n", matrix->numColumns,
          numColumns);
        CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
        if (inputSubmatrixFileName)
          CMR_CALL( CMRintmatFree(cmr, &explicitSubmatrix) );
        CMR_CALL( CMRintmatFree(cmr, &matrix) );
        CMR_CALL( CMRfreeEnvironment(&cmr) );
        return CMR_ERROR_INPUT;
      }

      explicitSubmatrix = NULL;
      CMR_CALL( CMRintmatZoomSubmat(cmr, matrix, submatrix, &explicitSubmatrix) );
      CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
    }
    else
      explicitSubmatrix = matrix;

    if (task == TASK_SUPPORT)
    {
      CMR_CHRMAT* result = NULL;
      CMR_CALL( CMRintmatSupport(cmr, explicitSubmatrix, &result) );
      CMR_CALL( writeToFileChr(cmr, result, outputFormat, outputMatrixFileName, transpose) );
      CMR_CALL( CMRchrmatFree(cmr, &result) );
    }
    else if (task == TASK_SIGNED_SUPPORT)
    {
      CMR_CHRMAT* result = NULL;
      CMR_CALL( CMRintmatSignedSupport(cmr, explicitSubmatrix, &result) );
      CMR_CALL( writeToFileChr(cmr, result, outputFormat, outputMatrixFileName, transpose) );
      CMR_CALL( CMRchrmatFree(cmr, &result) );
    }
    else
      CMR_CALL( writeToFileInt(cmr, explicitSubmatrix, outputFormat, outputMatrixFileName, transpose) );

    if (inputSubmatrixFileName)
      CMR_CALL( CMRintmatFree(cmr, &explicitSubmatrix) );
  }

  CMR_CALL( CMRintmatFree(cmr, &matrix) );
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
  fputs("  -S IN-SUB Consider the submatrix of IN-MAT specified in file IN-SUB instead of IN-MAT itself; can be combined with other operations.\n",
    stderr);
  fputs("  -t        Transpose the matrix; can be combined with other operations.\n", stderr);
  fputs("  -c        Compute the support matrix instead of copying.\n", stderr);
  fputs("  -C        Compute the signed support matrix instead of copying.\n", stderr);
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

  CMR_ERROR error;
  if (doubleArithmetic)
    error = runDbl(inputMatrixFileName, inputFormat, inputSubmatrixFileName, outputFormat, outputMatrixFileName, task, transpose);
  else
    error = runInt(inputMatrixFileName, inputFormat, inputSubmatrixFileName, outputFormat, outputMatrixFileName, task, transpose);
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
