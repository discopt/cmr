#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <tu/matrix.h>
#include <tu/graphic.h>

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

  if (!TUisTernaryChr(tu, matrix, NULL))
  {
    printf("Input matrix is not ternary.\n");
    TU_CALL( TUchrmatFree(tu, &matrix) );
    TU_CALL( TUfreeEnvironment(&tu) );
    return TU_OKAY;
  }

  TU_CHRMAT* transpose = NULL;
  TU_CALL( TUchrmatTranspose(tu, matrix, &transpose) );

  bool isGraphic;
  TU_GRAPH* graph = NULL;
  TU_GRAPH_EDGE* basis = NULL;
  TU_GRAPH_EDGE* cobasis = NULL;
  TU_CALL( TUtestTernaryGraphic(tu, transpose, &isGraphic, &graph, &basis, &cobasis, NULL, NULL) );

  if (isGraphic)
  {
    printf("Input matrix is graphic.\n");

    TU_CHRMAT* checkMatrix = NULL;
    TU_CALL( TUcomputeGraphBinaryRepresentationMatrix(tu, graph, &checkMatrix, NULL, matrix->numRows, basis,
      matrix->numColumns, cobasis, NULL) );

    if (!TUchrmatCheckEqual(matrix, checkMatrix))
      printf("ERROR: computed representation matrix does NOT agree with input matrix!\n");

    TU_CALL( TUchrmatFree(tu, &checkMatrix) );

    TU_CALL( TUfreeBlockArray(tu, &cobasis) );
    TU_CALL( TUfreeBlockArray(tu, &basis) );
    TU_CALL( TUgraphFree(tu, &graph) );
  }
  else
    printf("Input matrix is NOT graphic.\n");

  TU_CALL( TUchrmatFree(tu, &transpose) );
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
  default:
    return EXIT_SUCCESS;
  }
}
