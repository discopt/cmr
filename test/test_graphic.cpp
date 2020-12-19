#include <gtest/gtest.h>

#include "common.h"
#include <tu/graphic.h>

void testGraphicMatrix(
  TU* tu,             /**< \ref TU environment. */
  TU_CHRMAT* matrix   /**< Matrix to be used for testing. */
)
{
  TU_GRAPH* graph = NULL;
  TU_GRAPH_EDGE* basis = NULL;
  TU_GRAPH_EDGE* cobasis = NULL;
  bool isGraphic;
  ASSERT_FALSE( TUtestGraphicnessChr(tu, matrix, &isGraphic, &graph, &basis, &cobasis, NULL) );

  ASSERT_TRUE( isGraphic );
  ASSERT_TRUE( basis );
  ASSERT_TRUE( cobasis );

  TU_CHRMAT* result = NULL;
  TU_CALL( TUconvertGraphToBinaryMatrix(tu, graph, &result, matrix->numRows, basis,
    matrix->numColumns, cobasis) );

  ASSERT_TRUE( result );

  if (TUchrmatCheckEqual(matrix, result))
  {
    printf("The representation matrix of represented graph is equal to input matrix.\n");
  }
  else
  {
    printf("Input matrix:\n");
    ASSERT_FALSE( TUchrmatPrintDense(stdout, matrix, ' ', true) );

    printf("Graph:\n");
    ASSERT_FALSE( TUgraphPrint(stdout, graph) );

    printf("Representation matrix:\n");
    ASSERT_FALSE( TUchrmatPrintDense(stdout, result, ' ', true) );
  }
  

  ASSERT_FALSE( TUgraphFree(tu, &graph) );
  ASSERT_FALSE( TUfreeBlockArray(tu, &basis) );
  ASSERT_FALSE( TUfreeBlockArray(tu, &cobasis) );
  ASSERT_FALSE( TUchrmatFree(tu, &result) );
}

TEST(Graphic, Polygon)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* A = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "4 1 "
    "1 "
    "1 "
    "0 "
    "1 "
  ) );
  testGraphicMatrix(tu, A);
  ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, PolygonPlusEdge)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* A = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "4 2 "
    "1 0 "
    "1 1 "
    "1 0 "
    "1 1 "
  ) );
  testGraphicMatrix(tu, A);
  ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}
