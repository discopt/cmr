#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <float.h>

#include <cmr/matrix.h>
#include <cmr/ctu.h>

typedef enum
{
  TASK_RECOGNIZE = 1, /**< Determine if given matrix is ctu. */
  TASK_APPLY = 2      /**< Apply complement operations. */
} Task;

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

/**
 * \brief Tests matrix from a file for total unimodularity.
 */

static
CMR_ERROR testComplementTotalUnimodularity(
  const char* inputMatrixFileName,  /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,           /**< Format of the input matrix. */
  FileFormat outputFormat,          /**< Format of the output submatrix. */
  char* outputOperationsFileName,   /**< File name for the operations; may be `-' for stdout. */
  char* outputMatrixFileName,       /**< File name for the matrix; may be `-' for stdout. */
  bool printStats,                  /**< Whether to print statistics to stderr. */
  double timeLimit                  /**< Time limit to impose. */
)
{
  clock_t readClock = clock();
  FILE* inputMatrixFile = strcmp(inputMatrixFileName, "-") ? fopen(inputMatrixFileName, "r") : stdin;
  if (!inputMatrixFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_CHRMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, inputMatrixFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, inputMatrixFile, &matrix) );
  if (inputMatrixFile != stdin)
    fclose(inputMatrixFile);
  fprintf(stderr, "Read %zux%zu matrix with %zu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Actual test. */

  bool isCTU;
  size_t complementRow = SIZE_MAX;
  size_t complementColumn = SIZE_MAX;
  CMR_CTU_STATISTICS stats;
  CMR_CALL( CMRstatsComplementTotalUnimodularityInit(&stats) );
  CMR_CALL( CMRtestComplementTotalUnimodularity(cmr, matrix, &isCTU, &complementRow, &complementColumn, &stats,
    timeLimit) );

  fprintf(stderr, "Matrix %scomplement totally unimodular.\n", isCTU ? "IS " : "IS NOT ");
  if (printStats)
    CMR_CALL( CMRstatsComplementTotalUnimodularityPrint(stderr, &stats, NULL) );

  if (complementRow < SIZE_MAX || complementColumn < SIZE_MAX)
  {
    if (outputOperationsFileName)
    {
      bool outputOperationsToFile = strcmp(outputOperationsFileName, "-");
      FILE* outputOperationsFile = outputOperationsToFile ? fopen(outputOperationsFileName, "w") : stdout;
      fprintf(stderr, "Writing complement operations to %s%s%s.\n", outputOperationsToFile ? "file <" : "",
        outputOperationsToFile ? outputOperationsFileName : "stdout", outputOperationsToFile ? ">" : "");

      if (complementRow < SIZE_MAX)
        fprintf(outputOperationsFile, "Complement row %zu\n", complementRow + 1);
      if (complementColumn < SIZE_MAX)
        fprintf(outputOperationsFile, "Complement column %zu\n", complementColumn + 1);

      if (outputOperationsFile)
        fclose(outputOperationsFile);
    }

    if (outputMatrixFileName)
    {
      bool outputMatrixToFile = strcmp(outputMatrixFileName, "-");
      FILE* outputMatrixFile = outputMatrixToFile ? fopen(outputMatrixFileName, "w") : stdout;
      fprintf(stderr, "Writing complemented non-totally unimodular matrix to %s%s%s in %s format.\n",
        outputMatrixToFile ? "file <" : "", outputMatrixToFile ? outputMatrixFileName : "stdout",
        outputMatrixToFile ? ">" : "", outputFormat == FILEFORMAT_MATRIX_DENSE ? "dense" : "sparse");

      CMR_CHRMAT* complemented = NULL;
      CMR_CALL( CMRcomplementRowColumn(cmr, matrix, complementRow, complementColumn, &complemented) );

      if (outputFormat == FILEFORMAT_MATRIX_DENSE)
        CMR_CALL( CMRchrmatPrintDense(cmr, complemented, stdout, '0', false) );
      else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
        CMR_CALL( CMRchrmatPrintSparse(cmr, complemented, stdout) );
      else
        assert(false);
  
      CMR_CALL( CMRchrmatFree(cmr, &complemented) );

      if (outputMatrixToFile)
        fclose(outputMatrixFile);
    }
  }

  /* Cleanup. */

  CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}


static
CMR_ERROR complementMatrix(
  const char* inputMatrixFileName,  /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,           /**< Format of the input matrix. */
  FileFormat outputFormat,          /**< Format of the output submatrix. */
  size_t complementRow,             /**< Complement row. */
  size_t complementColumn,          /**< Complement column. */
  char* outputMatrixFileName        /**< File name for the matrix; may be `-' for stdout. */
)
{
  FILE* inputMatrixFile = strcmp(inputMatrixFileName, "-") ? fopen(inputMatrixFileName, "r") : stdin;
  if (!inputMatrixFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_CHRMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, inputMatrixFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, inputMatrixFile, &matrix) );
  if (inputMatrixFile != stdin)
    fclose(inputMatrixFile);
  fprintf(stderr, "Read %zux%zu matrix with %zu nonzeros.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros);

  /* Carry out complementing. */
  CMR_CHRMAT* complemented = NULL;
  CMR_CALL( CMRcomplementRowColumn(cmr, matrix, complementRow, complementColumn, &complemented) );

  bool outputMatrixToFile = strcmp(outputMatrixFileName, "-");
  FILE* outputMatrixFile = outputMatrixToFile ? fopen(outputMatrixFileName, "w") : stdout;
  fprintf(stderr, "Writing complemented matrix to %s%s%s in %s format.\n",
    outputMatrixToFile ? "file <" : "", outputMatrixToFile ? outputMatrixFileName : "stdout",
    outputMatrixToFile ? ">" : "", outputFormat == FILEFORMAT_MATRIX_DENSE ? "dense" : "sparse");

  if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatPrintDense(cmr, complemented, outputMatrixFile, '0', false) );
  else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatPrintSparse(cmr, complemented, outputMatrixFile) );

  if (outputMatrixToFile)
    fclose(outputMatrixFile);

  CMR_CALL( CMRchrmatFree(cmr, &complemented) );
  
  /* Cleanup. */

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
  fprintf(stderr, "%s IN-MAT [OPTION]...\n\n", program);
  fputs("  (1) determines whether the matrix given in file IN-MAT is complement totally unimodular.\n\n", stderr);
  fprintf(stderr, "%s -c IN-MAT OUT-MAT [OPTION]...\n\n", program);
  fputs("  (2) applies a sequence of row or column complement operations the matrix given in file IN-MAT and writes the result to OUT-MAT.\n\n\n",
    stderr);
  fputs("Options specific to (1):\n", stderr);
  fputs("  -n OUT-OPS  Write complement operations that leads to a non-totally-unimodular matrix to file OUT-OPS; default: skip computation.\n", stderr);
  fputs("  -N OUT-MAT  Write a complemented matrix that is non-totally-unimodular to file OUT-MAT; default: skip computation.\n", stderr);
  fputs("Options specific to (2):\n", stderr);
  fputs("  -r ROW    Apply row complement operation to row ROW.\n", stderr);
  fputs("  -c COLUMN Apply column complement operation to column COLUMN.\n", stderr);
  fputs("Common options:\n", stderr);
  fputs("  -i FORMAT   Format of file IN-MAT, among `dense' and `sparse'; default: dense.\n", stderr);
  fputs("  -o FORMAT   Format of file OUT-MAT, among `dense' and `sparse'; default: same as for IN-MAT.\n", stderr);
  fputs("  -s          Print statistics about the computation to stderr.\n\n", stderr);
  fputs("Advanced options:\n", stderr);
  fputs("  --time-limit LIMIT   Allow at most LIMIT seconds for the computation.\n\n", stderr);
  fputs("If IN-MAT is `-' then the matrix is read from stdin.\n", stderr);
  fputs("If OUT-OPS or OUT-MAT is `-` then the list of operations (resp. the matrix) is written to stdout.\n", stderr);

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  Task task = TASK_RECOGNIZE;
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  FileFormat outputFormat = FILEFORMAT_UNDEFINED;
  size_t complementRow = SIZE_MAX;
  size_t complementColumn = SIZE_MAX;
  char* inputMatrixFileName = NULL;
  char* outputMatrixFileName = NULL;
  char* outputOperationsFileName = NULL;
  bool printStats = false;
  double timeLimit = DBL_MAX;
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
        fprintf(stderr, "Error: Invalid complement row <%s>", argv[a+1]);
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
        fprintf(stderr, "Error: Invalid complement column <%s>", argv[a+1]);
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a++;
    }
    else if (!strcmp(argv[a], "-n") && a+1 < argc)
      outputOperationsFileName = argv[++a];
    else if (!strcmp(argv[a], "-N") && a+1 < argc)
    {
      if (outputMatrixFileName && task == TASK_APPLY)
      {
        fprintf(stderr, "Error: Option -N is invalid for operation application.\n\n");
        return printUsage(argv[0]);
      }
      outputMatrixFileName = argv[++a];
    }
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
        fprintf(stderr, "Error: Unknown input file format <%s>.\n\n", argv[a+1]);
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
    else if (!outputMatrixFileName)
      outputMatrixFileName = argv[a];
    else if (task == TASK_APPLY)
    {
      printf("Error: Two output files <%s> and <%s> specified.\n\n", outputMatrixFileName, argv[a]);
      return printUsage(argv[0]);
    }
    else
    {
      assert(task == TASK_RECOGNIZE);
      fprintf(stderr, "Error: Cannot choose to apply operations after supplying option -N.\n\n");
      return printUsage(argv[0]);
    }
  }

  if (!inputMatrixFileName)
  {
    fputs("Error: No input file specified.\n\n", stderr);
    return printUsage(argv[0]);
  }

  if (outputFormat == FILEFORMAT_UNDEFINED)
    outputFormat = inputFormat;

  CMR_ERROR error;
  if (task == TASK_RECOGNIZE)
  {
    error = testComplementTotalUnimodularity(inputMatrixFileName, inputFormat, outputFormat, outputOperationsFileName,
      outputMatrixFileName, printStats, timeLimit);
  }
  else
  {
    assert(task == TASK_APPLY);
    error = complementMatrix(inputMatrixFileName, inputFormat, outputFormat, complementRow, complementColumn,
      outputMatrixFileName);
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
