#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

TU_ERROR runDbl(const char* instanceFileName, Format inputFormat, Format outputFormat, Task task)
{
  FILE* instanceFile = fopen(instanceFileName, "r");
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
  fclose(instanceFile);

  if (task == SUPPORT)
  {
    TU_CHRMAT* result = NULL;
    TU_CALL( TUsupportDbl(tu, matrix, 1.0e-9, &result) );
    if (outputFormat == SPARSE)
      TU_CALL( TUchrmatPrintSparse(stdout, result) );
    else if (outputFormat == DENSE)
      TU_CALL( TUchrmatPrintDense(stdout, result, '0', false) );
    else
      return TU_ERROR_INPUT;
    TU_CALL( TUchrmatFree(tu, &result) );
  }
  else if (task == SIGNED_SUPPORT)
  {
    TU_CHRMAT* result = NULL;
    TU_CALL( TUsignedSupportDbl(tu, matrix, 1.0e-9, &result) );
    if (outputFormat == SPARSE)
      TU_CALL( TUchrmatPrintSparse(stdout, result) );
    else if (outputFormat == DENSE)
      TU_CALL( TUchrmatPrintDense(stdout, result, '0', false) );
    else
      return TU_ERROR_INPUT;
    TU_CALL( TUchrmatFree(tu, &result) );
  }
  else
  {
    if (outputFormat == SPARSE)
      TU_CALL( TUdblmatPrintSparse(stdout, matrix) );
    else if (outputFormat == DENSE)
      TU_CALL( TUdblmatPrintDense(stdout, matrix, '0', false) );
    else
      return TU_ERROR_INPUT;
  }

  TU_CALL( TUdblmatFree(tu, &matrix) );

  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
}

TU_ERROR runInt(const char* instanceFileName, Format inputFormat, Format outputFormat, Task task)
{
  FILE* instanceFile = fopen(instanceFileName, "r");
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
  fclose(instanceFile);

  if (task == SUPPORT)
  {
    TU_CHRMAT* result = NULL;
    TU_CALL( TUsupportInt(tu, matrix, &result) );
    if (outputFormat == SPARSE)
      TU_CALL( TUchrmatPrintSparse(stdout, result) );
    else if (outputFormat == DENSE)
      TU_CALL( TUchrmatPrintDense(stdout, result, '0', false) );
    else
      return TU_ERROR_INPUT;
    TU_CALL( TUchrmatFree(tu, &result) );
  }
  else if (task == SIGNED_SUPPORT)
  {
    TU_CHRMAT* result = NULL;
    TU_CALL( TUsignedSupportInt(tu, matrix, &result) );
    if (outputFormat == SPARSE)
      TU_CALL( TUchrmatPrintSparse(stdout, result) );
    else if (outputFormat == DENSE)
      TU_CALL( TUchrmatPrintDense(stdout, result, '0', false) );
    else
      return TU_ERROR_INPUT;
    TU_CALL( TUchrmatFree(tu, &result) );
  }
  else
  {
    if (outputFormat == SPARSE)
      TU_CALL( TUintmatPrintSparse(stdout, matrix) );
    else if (outputFormat == DENSE)
      TU_CALL( TUintmatPrintDense(stdout, matrix, '0', false) );
    else
      return TU_ERROR_INPUT;
  }

  TU_CALL( TUintmatFree(tu, &matrix) );

  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
}

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... MATRIX\n\n", program);
  puts("Copies MATRIX into another format.");
  puts("\nOptions:");
  puts("  -i FORMAT  Format of MATRIX file, among {dense, sparse}; default: dense.");
  puts("  -o FORMAT  Format of output, among {dense, sparse}; default: same as input.");
  puts("  -s         Create support matrix instead of copying.");
  puts("  -t         Create signed support matrix instead of copying.");
  puts("  -d         Use double arithmetic.");
  
  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  /* Parse command line options. */
  Format inputFormat = DENSE;
  Format outputFormat = UNDEFINED;
  Task task = COPY;
  bool doubleArithmetic = false;
  char* instanceFileName = NULL;
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
    else if (!strcmp(argv[a], "-o") && a+1 < argc)
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
    else if (!strcmp(argv[a], "-s"))
      task = SUPPORT;
    else if (!strcmp(argv[a], "-t"))
      task = SIGNED_SUPPORT;
    else if (!strcmp(argv[a], "-d"))
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
    error = runDbl(instanceFileName, inputFormat, outputFormat, task);
  else
    error = runInt(instanceFileName, inputFormat, outputFormat, task);
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
