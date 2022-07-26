#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <cmr/matrix.h>

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,       /**< Whether the file format of input/output was defined by the user. */
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

static
CMR_ERROR run(
  const char* instanceFileName,
  FileFormat inputFormat,
  bool testBinary,
  bool testTernary,
  bool testInteger,
  double epsilon)
{
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_DBLMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRdblmatCreateFromSparseStream(cmr, instanceFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRdblmatCreateFromDenseStream(cmr, instanceFile, &matrix) );
  else
    return CMR_ERROR_INPUT;
  if (instanceFile != stdin)
    fclose(instanceFile);

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

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... MATRIX\n\n", program);
  puts("Inspects entries of MATRIX.");
  puts("\nOptions:");
  puts("  -i FORMAT  Format of MATRIX file, among {dense, sparse}; default: dense.");
  puts("  -b         Tests whether the matrix is binary, i.e., has entries in {0,+1}.");
  puts("  -t         Tests whether the matrix is ternary, i.e., has entries in {-1,0,+1}.");
  puts("  -I         Tests whether the matrix is integer.");
  puts("  -t EPSILON Allows rounding of numbers up to tolerance EPSILON; default: 1.0e-9.");
  puts("If MATRIX is `-', then the matrix will be read from stdin.");
  
  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  /* Parse command line options. */
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  bool testInteger = false;
  bool testBinary = false;
  bool testTernary = false;
  char* instanceFileName = NULL;
  double epsilon = 1.0e-9;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-b"))
      testBinary = true;
    else if (!strcmp(argv[a], "-t"))
      testTernary = true;
    else if (!strcmp(argv[a], "-I"))
      testInteger = true;
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
    else if (!strcmp(argv[a], "-t") && a+1 < argc)
    {
      char* p;
      epsilon = strtod(argv[a+1], &p);
      if (*p != '\0' || epsilon < 0.0 || epsilon > 0.5)
      {
        fprintf(stderr, "Error: invalid tolerance <%s>", argv[a+1]);
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a++;
    }
    else if (!instanceFileName)
      instanceFileName = argv[a];
    else
    {
      printf("Error: Two input matrix files <%s> and <%s> specified.\n\n", instanceFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (!instanceFileName)
  {
    puts("No input matrix specified.\n");
    return printUsage(argv[0]);
  }

  CMR_ERROR error;
  error = run(instanceFileName, inputFormat, testBinary, testTernary, testInteger, epsilon);
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
