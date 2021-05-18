#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <tu/matrix.h>
#include <tu/sign.h>

typedef enum
{
  UNDEFINED = 0,
  DENSE = 1,
  SPARSE = 2
} Format;

typedef enum
{
  COPY = 0,
  SUPPORT = 1,
  SIGNED_SUPPORT = 2
} Task;

static
TU_ERROR printDbl(TU* tu, TU_DBLMAT* matrix, Format outputFormat, bool transpose)
{
  assert(matrix);

  TU_DBLMAT* output = NULL;
  if (transpose)
  {
    TU_CALL( TUdblmatTranspose(tu, matrix, &output) );
  }
  else
    output = matrix;

  TU_ERROR error = TU_OKAY;
  if (outputFormat == SPARSE)
    TU_CALL( TUdblmatPrintSparse(stdout, output) );
  else if (outputFormat == DENSE)
    TU_CALL( TUdblmatPrintDense(stdout, output, '0', false) );
  else
    error = TU_ERROR_INPUT;

  if (transpose)
    TU_CALL( TUdblmatFree(tu, &output) );

  return error;
}

static
TU_ERROR printInt(TU* tu, TU_INTMAT* matrix, Format outputFormat, bool transpose)
{
  assert(matrix);

  TU_INTMAT* output = NULL;
  if (transpose)
  {
    TU_CALL( TUintmatTranspose(tu, matrix, &output) );
  }
  else
    output = matrix;

  TU_ERROR error = TU_OKAY;
  if (outputFormat == SPARSE)
    TU_CALL( TUintmatPrintSparse(stdout, output) );
  else if (outputFormat == DENSE)
    TU_CALL( TUintmatPrintDense(stdout, output, '0', false) );
  else
    error = TU_ERROR_INPUT;

  if (transpose)
    TU_CALL( TUintmatFree(tu, &output) );

  return error;
}

static
TU_ERROR printChr(TU* tu, TU_CHRMAT* matrix, Format outputFormat, bool transpose)
{
  assert(matrix);

  TU_CHRMAT* output = NULL;
  if (transpose)
  {
    TU_CALL( TUchrmatTranspose(tu, matrix, &output) );
  }
  else
    output = matrix;

  TU_ERROR error = TU_OKAY;
  if (outputFormat == SPARSE)
    TU_CALL( TUchrmatPrintSparse(stdout, output) );
  else if (outputFormat == DENSE)
    TU_CALL( TUchrmatPrintDense(stdout, output, '0', false) );
  else
    error = TU_ERROR_INPUT;

  if (transpose)
    TU_CALL( TUchrmatFree(tu, &output) );

  return error;
}

TU_ERROR runDbl(const char* instanceFileName, Format inputFormat, Format outputFormat, Task task, bool transpose)
{
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return TU_ERROR_INPUT;

  TU* tu = NULL;
  TU_CALL( TUcreateEnvironment(&tu) );

  TU_DBLMAT* matrix = NULL;
  if (inputFormat == SPARSE)
    TU_CALL( TUdblmatCreateFromSparseStream(tu, &matrix, instanceFile) );
  else if (inputFormat == DENSE)
    TU_CALL( TUdblmatCreateFromDenseStream(tu, &matrix, instanceFile) );
  else
    return TU_ERROR_INPUT;
  if (instanceFile != stdin)
    fclose(instanceFile);

  if (task == SUPPORT)
  {
    TU_CHRMAT* result = NULL;
    TU_CALL( TUsupportDbl(tu, matrix, 1.0e-9, &result) );
    TU_CALL( printChr(tu, result, outputFormat, transpose) );
    TU_CALL( TUchrmatFree(tu, &result) );
  }
  else if (task == SIGNED_SUPPORT)
  {
    TU_CHRMAT* result = NULL;
    TU_CALL( TUsignedSupportDbl(tu, matrix, 1.0e-9, &result) );
    TU_CALL( printChr(tu, result, outputFormat, transpose) );
    TU_CALL( TUchrmatFree(tu, &result) );
  }
  else
  {
    TU_CALL( printDbl(tu, matrix, outputFormat, transpose) );
  }

  TU_CALL( TUdblmatFree(tu, &matrix) );

  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
}

TU_ERROR runInt(const char* instanceFileName, Format inputFormat, Format outputFormat, Task task, bool transpose)
{
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return TU_ERROR_INPUT;

  TU* tu = NULL;
  TU_CALL( TUcreateEnvironment(&tu) );

  TU_INTMAT* matrix = NULL;
  if (inputFormat == SPARSE)
    TU_CALL( TUintmatCreateFromSparseStream(tu, &matrix, instanceFile) );
  else if (inputFormat == DENSE)
    TU_CALL( TUintmatCreateFromDenseStream(tu, &matrix, instanceFile) );
  else
    return TU_ERROR_INPUT;
  if (instanceFile != stdin)
    fclose(instanceFile);

  if (task == SUPPORT)
  {
    TU_CHRMAT* result = NULL;
    TU_CALL( TUsupportInt(tu, matrix, &result) );
    TU_CALL( printChr(tu, result, outputFormat, transpose) );
    TU_CALL( TUchrmatFree(tu, &result) );
  }
  else if (task == SIGNED_SUPPORT)
  {
    TU_CHRMAT* result = NULL;
    TU_CALL( TUsignedSupportInt(tu, matrix, &result) );
    TU_CALL( printChr(tu, result, outputFormat, transpose) );
    TU_CALL( TUchrmatFree(tu, &result) );
  }
  else
  {
    TU_CALL( printInt(tu, matrix, outputFormat, transpose) );
  }

  TU_CALL( TUintmatFree(tu, &matrix) );

  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
}

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... MATRIX\n\n", program);
  puts("Copies MATRIX, potentially applying an operation.");
  puts("\nOptions:");
  puts("  -i, --input FORMAT  Format of MATRIX file, among {dense, sparse}; default: dense.");
  puts("  -o, --output FORMAT Format of output, among {dense, sparse}; default: same as input.");
  puts("  -s, --support       Create support matrix instead of copying.");
  puts("  -t, --transpose     Output transposed matrix (can be combined with other operations).");
  puts("  -S, --sign          Create signed support matrix instead of copying.");
  puts("  -d, --double        Use double arithmetic.");
  puts("If MATRIX is `-', then the matrix will be read from stdin.");
  
  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  /* Parse command line options. */
  Format inputFormat = DENSE;
  Format outputFormat = UNDEFINED;
  Task task = COPY;
  bool transpose = false;
  bool doubleArithmetic = false;
  char* instanceFileName = NULL;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h") || !strcmp(argv[a], "--help"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if ((!strcmp(argv[a], "-i") || !strcmp(argv[a], "--input")) && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        inputFormat = DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        inputFormat = SPARSE;
      else
      {
        printf("Error: unknown input format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if ((!strcmp(argv[a], "-o") || !strcmp(argv[a], "--output")) && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        outputFormat = DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        outputFormat = SPARSE;
      else
      {
        printf("Error: unknown output format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-s") || !strcmp(argv[a], "--support"))
      task = SUPPORT;
    else if (!strcmp(argv[a], "-S") || !strcmp(argv[a], "--sign"))
      task = SIGNED_SUPPORT;
    else if (!strcmp(argv[a], "-t") || !strcmp(argv[a], "--transpose"))
      transpose = true;
    else if (!strcmp(argv[a], "-d") || !strcmp(argv[a], "--double"))
      doubleArithmetic = true;
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
  if (outputFormat == UNDEFINED)
    outputFormat = inputFormat;

  TU_ERROR error;
  if (doubleArithmetic)
    error = runDbl(instanceFileName, inputFormat, outputFormat, task, transpose);
  else
    error = runInt(instanceFileName, inputFormat, outputFormat, task, transpose);
  switch (error)
  {
  case TU_ERROR_INPUT:
    puts("Input error.");
    return EXIT_FAILURE;
  case TU_ERROR_MEMORY:
    puts("Memory error.");
    return EXIT_FAILURE;
  default:
    return EXIT_SUCCESS;
  }
}
