#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <cmr/matrix.h>
#include <cmr/regular.h>
#include <cmr/graphic.h>
#include <cmr/network.h>

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... MATRIX\n\n", program);
  puts("Tests matrix for total unimodularity or related properties.");
  puts("\nOptions:");
  puts("  -i FORMAT  Format of MATRIX file, among {dense, sparse}; default: dense.");

  return EXIT_FAILURE;
}

CMR_ERROR run(const char* instanceFileName, bool sparse)
{
  FILE* instanceFile = fopen(instanceFileName, "r");
  if (!instanceFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  if (sparse)
    CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, &matrix, instanceFile) );
  else
    CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, &matrix, instanceFile) );
  fclose(instanceFile);

  CMR_TU_DEC* dec = NULL;
  CMR_CALL( CMRtestBinaryRegular(cmr, matrix, NULL, NULL, true, true, NULL, &dec) );

  CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  CMR_CHRMAT* transpose = NULL;
  CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

  bool isGraphic;
  CMR_GRAPH* graph = NULL;
  CMR_GRAPH_EDGE* basis = NULL;
  CMR_GRAPH_EDGE* cobasis = NULL;
  CMR_CALL( CMRtestNetworkMatrix(cmr, transpose, &isGraphic, &graph, &basis, &cobasis, NULL, NULL) );

  if (isGraphic)
  {
    printf("Input matrix is graphic.\n");

    CMR_CHRMAT* checkMatrix = NULL;
    CMR_CALL( CMRcomputeGraphicMatrix(cmr, graph, &checkMatrix, NULL, matrix->numRows, basis,
      matrix->numColumns, cobasis, NULL) );

    if (!CMRchrmatCheckEqual(matrix, checkMatrix))
      printf("ERROR: computed representation matrix does NOT agree with input matrix!\n");

    CMR_CALL( CMRchrmatFree(cmr, &checkMatrix) );

    CMR_CALL( CMRfreeBlockArray(cmr, &cobasis) );
    CMR_CALL( CMRfreeBlockArray(cmr, &basis) );
    CMR_CALL( CMRgraphFree(cmr, &graph) );
  }
  else
    printf("Input matrix is NOT graphic.\n");

  CMR_CALL( CMRchrmatFree(cmr, &transpose) );
  CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
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

  CMR_ERROR error = run(instanceFileName, sparse);
  switch (error)
  {
  case CMR_ERROR_INPUT:
    printf("Input error in file <%s>\n", instanceFileName);
    return EXIT_FAILURE;
  case CMR_ERROR_MEMORY:
    printf("Memory error.\n");
    return EXIT_FAILURE;
  case CMR_ERROR_INVALID:
    printf("Invalid call.\n");
    return EXIT_FAILURE;
  default:
    return EXIT_SUCCESS;
  }
}
