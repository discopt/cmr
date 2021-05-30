#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <tu/matrix.h>
#include <tu/regular.h>

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... MATRIX\n\n", program);
  puts("Tests matrix for total unimodularity or related properties.");
  puts("\nOptions:");
  puts("  -i FORMAT  Format of MATRIX file, among {dense, sparse}; default: dense.");

  return EXIT_FAILURE;
}

TU_ERROR run(const char* instanceFileName, bool sparse)
{
  FILE* instanceFile = fopen(instanceFileName, "r");
  if (!instanceFile)
    return TU_ERROR_INPUT;

  TU* tu = NULL;
  TU_CALL( TUcreateEnvironment(&tu) );

  TU_CHRMAT* matrix = NULL;
  if (sparse)
    TU_CALL( TUchrmatCreateFromSparseStream(tu, &matrix, instanceFile) );
  else
    TU_CALL( TUchrmatCreateFromDenseStream(tu, &matrix, instanceFile) );
  fclose(instanceFile);

  TU_DEC* dec = NULL;
  TU_CALL( TUtestBinaryRegular(tu, matrix, NULL, NULL, true, true, NULL, &dec) );

  TU_CALL( TUchrmatFree(tu, &matrix) );

  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
}

int main(int argc, char** argv)
{
  bool sparse = false;
  char* instanceFileName = NULL;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-s"))
      sparse = true;
    else if (!instanceFileName)
      instanceFileName = argv[a];
    else
    {
      printf("Multiple input files <%s> and <%s> specified.\n", instanceFileName, argv[a]);
      return EXIT_FAILURE;
    }
  }

  if (!instanceFileName)
  {
    printf("No input file specified.\n");
    return EXIT_FAILURE;
  }

  TU_ERROR error = run(instanceFileName, sparse);
  switch (error)
  {
  case TU_ERROR_INPUT:
    printf("Input error in file <%s>\n", instanceFileName);
    return EXIT_FAILURE;
  case TU_ERROR_MEMORY:
    printf("Memory error.\n");
    return EXIT_FAILURE;
  case TU_ERROR_INVALID:
    printf("Invalid call.\n");
    return EXIT_FAILURE;
  default:
    return EXIT_SUCCESS;
  }
}
