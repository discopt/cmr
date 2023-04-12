#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <float.h>

#include <cmr/matrix.h>
#include <cmr/camion.h>

typedef enum
{
  TASK_CHECK = 1, /**< Check for being Camion-signed. */
  TASK_SIGN = 2   /**< Camion-sign the matrix. */
} Task;

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,       /**< Whether the file format of input/output was defined by the user. */
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

/**
 * \brief Tests matrix from a file for being Camion-signed.
 */

static
CMR_ERROR checkCamionSigned(
  const char* inputMatrixFileName,      /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,               /**< Format of the input matrix. */
  const char* outputSubmatrixFileName,  /**< File name of output file for non-Camion submatrix. */
  bool printStats,                      /**< Whether to print statistics to stderr. */
  double timeLimit                      /**< Time limit to impose. */
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
  fprintf(stderr, "Read %lux%lu matrix with %lu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Actual test. */

  bool isCamion;
  CMR_SUBMAT* submatrix = NULL;
  CMR_CAMION_STATISTICS stats;
  CMR_CALL( CMRstatsCamionInit(&stats) );
  CMR_CALL( CMRtestCamionSigned(cmr, matrix, &isCamion,
    outputSubmatrixFileName ? &submatrix : NULL, printStats ? &stats : NULL, timeLimit) );

  fprintf(stderr, "Matrix %sCamion-signed.\n", isCamion ? "IS " : "IS NOT ");
  if (printStats)
    CMR_CALL( CMRstatsCamionPrint(stderr, &stats, NULL) );

  if (submatrix)
  {
    if (outputSubmatrixFileName)
    {
      bool outputSubmatrixToFile = strcmp(outputSubmatrixFileName, "-");
      fprintf(stderr, "Writing minimal non-Camion submatrix to %s%s%s.\n", outputSubmatrixToFile ? "file <" : "",
        outputSubmatrixToFile ? outputSubmatrixFileName : "stdout", outputSubmatrixToFile ? ">" : "");

      assert(submatrix);
      CMR_CALL( CMRsubmatWriteToFile(cmr, submatrix, matrix->numRows, matrix->numColumns, outputSubmatrixFileName) );
    }
  }

  /* Cleanup. */

  CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
  CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}


/**
 * \brief Camion-signs a matrix from a file.
 */

static
CMR_ERROR computeCamionSigned(
  const char* inputMatrixFileName,  /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,           /**< Format of the input matrix. */
  const char* outputMatrixFileName, /**< File name of output file for Camion-signed matrix. */
  FileFormat outputFormat,          /**< Format of the output matrix. */
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
  fprintf(stderr, "Read %lux%lu matrix with %lu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Actual signing. */
  
  CMR_CAMION_STATISTICS stats;
  CMR_CALL( CMRstatsCamionInit(&stats) );
  CMR_CALL( CMRcomputeCamionSigned(cmr, matrix, NULL, NULL, &stats, timeLimit) );
  if (printStats)
    CMR_CALL( CMRstatsCamionPrint(stderr, &stats, NULL) );

  /* Write to file. */

  bool outputMatrixToFile = strcmp(outputMatrixFileName, "-");
  FILE* outputMatrixFile = outputMatrixToFile ? fopen(outputMatrixFileName, "w") : stdout;
  fprintf(stderr, "Writing Camion-signed matrix to %s%s%s in %s format.\n", outputMatrixToFile ? "file <" : "",
    outputMatrixToFile ? outputMatrixFileName : "stdout", outputMatrixToFile ? ">" : "",
    outputFormat == FILEFORMAT_MATRIX_DENSE ? "dense" : "sparse");
  if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatPrintDense(cmr, matrix, outputMatrixFile, '0', false) );
  else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatPrintSparse(cmr, matrix, outputMatrixFile) );
  else
    assert(false);
  if (outputMatrixToFile)
    fclose(outputMatrixFile);

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
  fputs("  (1) determines whether the matrix given in file IN-MAT is Camion-signed.\n\n", stderr);
  fprintf(stderr, "%s IN-MAT -S OUT-MAT [OPTION]...\n\n", program);
  fputs("  (2) modifies the signs of the matrix given in file IN-MAT such that it is Camion-signed and writes the resulting new matrix to file OUT-MAT.\n\n\n",
    stderr);
  fputs("Options specific to (1):\n", stderr);
  fputs("  -N NON-SUB   Write a minimal non-Camion submatrix to file NON-SUB; default: skip computation.\n\n", stderr);
  fputs("Options specific to (2):\n", stderr);
  fputs("  -o FORMAT    Format of file OUT-MAT, among `dense' and `sparse'; default: same as format of IN-MAT.\n\n",
    stderr);
  fputs("Common options:\n", stderr);
  fputs("  -i FORMAT    Format of file IN-MAT, among `dense' and `sparse'; default: dense.\n", stderr);
  fputs("  -s           Print statistics about the computation to stderr.\n\n", stderr);
  fputs("Advanced options:\n", stderr);
  fputs("  --time-limit LIMIT   Allow at most LIMIT seconds for the computation.\n\n", stderr);
  fputs("If IN-MAT is `-' then the matrix is read from stdin.\n", stderr);
  fputs("If NON-SUB or OUT-MAT is `-' then the submatrix (resp. the Camion-signed matrix) is written to stdout.\n",
    stderr);

  return CMR_OKAY;
}

int main(int argc, char** argv)
{
  Task task = TASK_CHECK;
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  FileFormat outputFormat = FILEFORMAT_UNDEFINED;
  char* inputMatrixFileName = NULL;
  char* outputSubmatrixFileName = NULL;
  char* outputMatrixFileName = NULL;
  bool printStats = false;
  double timeLimit = DBL_MAX;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-N") && a+1 < argc)
      outputSubmatrixFileName = argv[++a];
    else if (!strcmp(argv[a], "-S") && a+1 < argc)
    {
      outputMatrixFileName = argv[++a];
      task = TASK_SIGN;
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
    else
    {
      printf("Error: Two input files <%s> and <%s> specified.\n\n", inputMatrixFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (!inputMatrixFileName)
  {
    fputs("Error: No input file specified.\n\n", stderr);
    return printUsage(argv[0]);
  }
  if (outputSubmatrixFileName && outputMatrixFileName)
  {
    fputs("Error: Options -S and -N cannot be specified at the same time.\n\n", stderr);
    return printUsage(argv[0]);
  }

  if (outputFormat == FILEFORMAT_UNDEFINED)
    outputFormat = inputFormat;
  
  CMR_ERROR error;
  if (task == TASK_CHECK)
  {
    error = checkCamionSigned(inputMatrixFileName, inputFormat, outputSubmatrixFileName, printStats, timeLimit);
  }
  else if (task == TASK_SIGN)
  {
    error = computeCamionSigned(inputMatrixFileName, inputFormat, outputMatrixFileName, outputFormat, printStats,
      timeLimit);
  }
  else
    assert(false);

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
