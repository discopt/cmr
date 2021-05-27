#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <tu/matrix.h>
#include <tu/graphic.h>

TU_ERROR findLargeTernarySubmatrix(TU* tu, TU_DBLMAT* matrix, TU_CHRMAT** result)
{
  assert(tu);
  assert(matrix);
  assert(result);

  return TU_OKAY;
}

TU_ERROR run(const char* instanceFileName)
{
  FILE* instanceFile = fopen(instanceFileName, "r");
  if (!instanceFile)
    return TU_ERROR_INPUT;

  TU* tu = NULL;
  TU_CALL( TUcreateEnvironment(&tu) );

  TU_DBLMAT* inputMatrix = NULL;
  TU_CALL( TUdblmatCreateFromDenseStream(tu, &inputMatrix, instanceFile) );
  fclose(instanceFile);
  TU_CALL( TUdblmatPrintDense(stdout, inputMatrix, '0', true) );

  TU_CHRMAT* ternaryMatrix = NULL;
  TU_CALL( findLargeTernarySubmatrix(tu, inputMatrix, &ternaryMatrix) );

  TU_CHRMAT* transpose = NULL;
  TU_CALL( TUchrmatTranspose(tu, ternaryMatrix, &transpose) );

  TU_SUBMAT* submatrix = NULL;
  TU_CALL( TUtestBinaryGraphicColumnSubmatrixGreedy(tu, transpose, NULL, &submatrix) );

  TU_CHRMAT* graphicMatrix = NULL;
  TU_CALL( TUchrmatFilterSubmat(tu, ternaryMatrix, submatrix, &graphicMatrix) );
  TU_CALL( TUsubmatFree(tu, &submatrix) );

//   TU_CALL( TUchrmatPrintDense(stdout, graphicMatrix, '0', true) );

  printf("Found a %dx%d graphic submatrix.\n", graphicMatrix->numRows, graphicMatrix->numColumns);

  TU_CALL( TUchrmatFree(tu, &transpose) );
  TU_CALL( TUchrmatFree(tu, &ternaryMatrix) );
  TU_CALL( TUchrmatFree(tu, &graphicMatrix) );
  TU_CALL( TUdblmatFree(tu, &inputMatrix) );

  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
}

int main(int argc, char** argv)
{
  if (argc == 1)
  {
    printf("No input file specified.\n");
    return EXIT_FAILURE;
  }

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
