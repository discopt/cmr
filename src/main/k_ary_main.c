#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include <cmr/matrix.h>

typedef enum
{
  TASK_RECOGNIZE = 1, /**< Check whether the matrix is binary, ternary or integer. */
  TASK_SUBMATRIX = 2  /**< Find a large binary or ternarry submatrix. */
} Task;

typedef enum
{
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

static
CMR_ERROR recognize(
  const char* inputMatrixFileName,
  FileFormat inputFormat,
  bool testBinary,
  bool testTernary,
  bool testInteger,
  double epsilon)
{
  clock_t readClock = clock();
  FILE* inputMatrixFile = strcmp(inputMatrixFileName, "-") ? fopen(inputMatrixFileName, "r") : stdin;
  if (!inputMatrixFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_DBLMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRdblmatCreateFromDenseStream(cmr, inputMatrixFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRdblmatCreateFromSparseStream(cmr, inputMatrixFile, &matrix) );
  if (inputMatrixFile != stdin)
    fclose(inputMatrixFile);
  fprintf(stderr, "Read %zux%zu matrix with %zu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  bool isInteger = testInteger;
  bool isBinary = testBinary;
  bool isTernary = testTernary;

  for (size_t e = 0; e < matrix->numNonzeros && (isInteger || isBinary || isTernary); ++e)
  {
    double value = matrix->entryValues[e];
    double rounded = round(value);
    if (fabs(value - rounded) > epsilon)
    {
      isInteger = isBinary = isTernary = false;
    }
    else if (rounded >= 1.5 || rounded <= -1.5)
    {
      isBinary = isTernary = false;
    }
    else if (rounded <= -0.5)
    {
      isBinary = false;
    }
  }

  if (testInteger)
  {
    printf("Matrix IS%s integer.\n", isInteger ? "" : " NOT");
  }
  if (testTernary)
  {
    printf("Matrix IS%s ternary.\n", isTernary ? "" : " NOT");
  }
  if (testBinary)
  {
    printf("Matrix IS%s binary.\n", isBinary ? "" : " NOT");
  }

  CMR_CALL( CMRdblmatFree(cmr, &matrix) );

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

static
CMR_ERROR findLargeSubmatrix(
  const char* inputMatrixFileName,
  FileFormat inputFormat,
  bool ternary,
  double epsilon,
  const char* outputSubmatrixFileName
)
{
  clock_t readClock = clock();
  FILE* inputMatrixFile = strcmp(inputMatrixFileName, "-") ? fopen(inputMatrixFileName, "r") : stdin;
  if (!inputMatrixFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_DBLMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRdblmatCreateFromDenseStream(cmr, inputMatrixFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRdblmatCreateFromSparseStream(cmr, inputMatrixFile, &matrix) );
  if (inputMatrixFile != stdin)
    fclose(inputMatrixFile);
  fprintf(stderr, "Read %zux%zu matrix with %zu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  CMR_SUBMAT* submatrix = NULL;
  if (ternary)
  {
    CMR_CALL( CMRdblmatFindTernarySubmatrix(cmr, matrix, epsilon, &submatrix) );
    bool outputSubmatrixToFile = strcmp(outputSubmatrixFileName, "-");
    fprintf(stderr, "Writing large ternary submatrix to %s%s%s.\n", outputSubmatrixToFile ? "file <" : "",
      outputSubmatrixToFile ? outputSubmatrixFileName : "stdout", outputSubmatrixToFile ? ">" : "");    
    CMR_CALL( CMRsubmatWriteToFile(cmr, submatrix, matrix->numRows, matrix->numColumns, outputSubmatrixFileName) );
  }
  else
  {
    CMR_CALL( CMRdblmatFindBinarySubmatrix(cmr, matrix, epsilon, &submatrix) );
    bool outputSubmatrixToFile = strcmp(outputSubmatrixFileName, "-");
    fprintf(stderr, "Writing large binary submatrix to %s%s%s.\n", outputSubmatrixToFile ? "file <" : "",
      outputSubmatrixToFile ? outputSubmatrixFileName : "stdout", outputSubmatrixToFile ? ">" : "");    
    CMR_CALL( CMRsubmatWriteToFile(cmr, submatrix, matrix->numRows, matrix->numColumns, outputSubmatrixFileName) );
  }

  CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
  CMR_CALL( CMRdblmatFree(cmr, &matrix) );
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
  fputs("  (1) determines whether the matrix given in file IN-MAT is integer (resp. binary or ternary).\n", stderr);
  fputs("\n", stderr);

  fprintf(stderr, "%s IN-MAT -R OUT-SUB [OPTION]...\n", program);
  fputs("  (2) finds a large binary (resp. ternary) submatrix of the matrix given in file IN-MAT.\n",
    stderr);
  fputs("\n", stderr);

  fputs("Options specific to (1):\n", stderr);
  fputs("  -b         Test whether the matrix is binary, i.e., has entries in {0,+1}.\n", stderr);
  fputs("  -t         Test whether the matrix is ternary, i.e., has entries in {-1,0,+1}.\n", stderr);
  fputs("  -I         Test whether the matrix is integer.\n", stderr);
  fputs("\n", stderr);

  fputs("Options specific to (2):\n", stderr);
  fputs("  -b         Find a large binary submatrix, i.e., one with only entries in {0,+1}.\n", stderr);
  fputs("  -t         Find a large ternary submatrix, i.e., one with only entries in {-1,0,+1}.\n", stderr);
  fputs("\n", stderr);

  fputs("Common options:\n", stderr);
  fputs("  -i FORMAT    Format of file IN-MAT; default: dense.\n", stderr);
  fputs("  -e EPSILON   Allows rounding of numbers up to tolerance EPSILON; default: 1.0e-9.\n", stderr);
  fputs("\n", stderr);

  fputs("Advanced options:\n", stderr);
  fputs("  --stats              Print statistics about the computation to stderr.\n", stderr);
  fputs("  --time-limit LIMIT   Allow at most LIMIT seconds for the computation.\n", stderr);
  fputs("\n", stderr);

  fputs("Formats for matrices: dense, sparse\n", stderr);
  fputs("If IN-MAT is `-' then the matrix is read from stdin.\n", stderr);
  fputs("If OUT-SUB is `-' then the submatrix is written to stdout.\n", stderr);

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  Task task = TASK_RECOGNIZE;
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  bool integer = false;
  bool binary = false;
  bool ternary = false;
  char* inputMatrixFileName = NULL;
  char* outputSubmatrixFileName = NULL;
  double epsilon = 1.0e-9;
  double timeLimit = DBL_MAX;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-b"))
      binary = true;
    else if (!strcmp(argv[a], "-t"))
      ternary = true;
    else if (!strcmp(argv[a], "-I"))
      integer = true;
    else if (!strcmp(argv[a], "-i") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        inputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        inputFormat = FILEFORMAT_MATRIX_SPARSE;
      else
      {
        fprintf(stderr, "Error: Unknown input file format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-e") && a+1 < argc)
    {
      char* p;
      epsilon = strtod(argv[a+1], &p);
      if (*p != '\0' || epsilon < 0.0 || epsilon > 0.5)
      {
        fprintf(stderr, "Error: Invalid tolerance <%s>", argv[a+1]);
        return printUsage(argv[0]);
      }
      a++;
    }
    else if (!strcmp(argv[a], "-R") && a+1 < argc)
    {
      outputSubmatrixFileName = argv[++a];
      task = TASK_SUBMATRIX;
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
      printf("Error: Two input matrix files <%s> and <%s> specified.\n\n", inputMatrixFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (!inputMatrixFileName)
  {
    fputs("Error: No input matrix specified.\n\n", stderr);
    return printUsage(argv[0]);
  }

  CMR_ERROR error;
  if (task == TASK_RECOGNIZE)
  {
    error = recognize(inputMatrixFileName, inputFormat, binary, ternary, integer, epsilon);
  }
  else if (task == TASK_SUBMATRIX)
  {
    if (integer)
    {
      fputs("Error: Option -I is invalid for large submatrix search.\n\n", stderr);
      return printUsage(argv[0]);
    }
    if ((binary ? 1 : 0) + (ternary ? 1 : 0) != 1)
    {
      fputs("Error: Either -b or -t must be specified.\n\n", stderr);
      return printUsage(argv[0]);
    }

    error = findLargeSubmatrix(inputMatrixFileName, inputFormat, ternary, epsilon, outputSubmatrixFileName);
  }
  else
  {
    error = CMR_ERROR_INPUT;
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
