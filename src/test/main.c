#include <stdlib.h>
#include <stdio.h>

#include <tu/matrix.h>
#include <tu/graphic.h>

TU_ERROR run(const char* instanceFileName)
{
  FILE* instanceFile = fopen(instanceFileName, "r");

  TU* tu = NULL;
  TU_CALL( TUcreateEnvironment(&tu) );

  TU_CHRMAT* matrix = NULL;
  TU_CALL( TUchrmatCreateFromDenseStream(tu, &matrix, instanceFile) );
  fclose(instanceFile);

  TU_CHRMAT* transpose = NULL;
  TU_CALL( TUchrmatTranspose(tu, matrix, &transpose) );

  bool isGraphic;
  TU_GRAPH* graph = NULL;
  TU_GRAPH_EDGE* basis = NULL;
  TU_GRAPH_EDGE* cobasis = NULL;
  TU_CALL( TUtestGraphicness(tu, transpose, &isGraphic, &graph, &basis, &cobasis, NULL) );

  if (isGraphic)
  {
    printf("Input matrix is graphic.\n");

    TU_CHRMAT* checkMatrix = NULL;
    TU_CALL( TUconvertGraphToBinaryMatrix(tu, graph, &checkMatrix, matrix->numRows, basis, matrix->numColumns,
      cobasis) );

    if (TUchrmatCheckEqual(matrix, checkMatrix))
      printf("Check: computed representation matrix agrees with input matrix.\n");
    else
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
  TU_ERROR error = run(argv[1]);
  switch (error)
  {
  case TU_ERROR_INPUT:
    printf("Input error in file <%s>\n", argv[1]);
    return EXIT_FAILURE;
  case TU_ERROR_MEMORY:
    printf("Memory error.\n");
    return EXIT_FAILURE;
  default:
    return EXIT_SUCCESS;
  }
}
