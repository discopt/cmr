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

TU_ERROR run(const char* instanceFileName, Format inputFormat, Format outputFormat)
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

  TU_CHRMAT* signedSupport = NULL;
  TU_CALL( TUsignedSupportInt(tu, matrix, &signedSupport) );
  TU_CALL( TUintmatFree(tu, &matrix) );

  bool alreadySigned;
  TU_CALL( TUcorrectSignChr(tu, signedSupport, &alreadySigned, NULL) );

  if (outputFormat == SPARSE)
    TU_CALL( TUchrmatPrintSparse(stdout, signedSupport) );
  else if (outputFormat == DENSE)
    TU_CALL( TUchrmatPrintDense(stdout, signedSupport, '0', false) );
  else
    return TU_ERROR_INPUT;

  TU_CALL( TUchrmatFree(tu, &signedSupport) );

  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
}

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... MATRIX\n\n", program);
  puts("Outputs the support matrix of the integer MATRIX, signed according to Camion's criterion.");
  puts("If MATRIX is already properly signed, then it won't be modified.");
  puts("\nOptions:");
  puts("  -i FORMAT  Format of MATRIX file, among {dense, sparse}; default: dense.");
  puts("  -o FORMAT  Format of output, among {dense, sparse}; default: same as input.");

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  /* Parse command line options. */
  Format inputFormat = DENSE;
  Format outputFormat = UNDEFINED;
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

  TU_ERROR error = run(instanceFileName, inputFormat, outputFormat);
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
