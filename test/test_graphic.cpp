#include <gtest/gtest.h>

#include <stdlib.h>

#include "common.h"
#include <tu/graphic.h>
#include <tu/tdec.h>

void testGraphicMatrix(
  TU* tu,             /**< \ref TU environment. */
  TU_CHRMAT* matrix,  /**< Matrix to be used for testing. */
  int mergeLeafBonds  /**< Leaf bonds of the t-decomposition are merged (1: at the end; 2: after each column). */
)
{
  TU_GRAPH* graph = NULL;
  ASSERT_TU_CALL( TUgraphCreateEmpty(tu, &graph, 0, 0) );
  TU_GRAPH_EDGE* basis = NULL;
  ASSERT_TU_CALL( TUallocBlockArray(tu, &basis, matrix->numRows) );
  TU_GRAPH_EDGE* cobasis = NULL;
  ASSERT_TU_CALL( TUallocBlockArray(tu, &cobasis, matrix->numColumns) );
  bool isGraphic;
  TU_CHRMAT* transpose = NULL;
  ASSERT_TU_CALL( TUchrmatTranspose(tu, matrix, &transpose) );

  ASSERT_TU_CALL( testGraphicnessTDecomposition(tu, matrix, transpose, &isGraphic, graph, basis,
    cobasis, NULL, mergeLeafBonds) );

  ASSERT_TRUE( isGraphic );
  ASSERT_TRUE( basis );
  ASSERT_TRUE( cobasis );

  ASSERT_TU_CALL( TUchrmatFree(tu, &transpose) );

  TU_CHRMAT* result = NULL;
  ASSERT_TU_CALL( TUconvertGraphToBinaryMatrix(tu, graph, &result, matrix->numRows, basis,
    matrix->numColumns, cobasis) );
  ASSERT_TRUE( result );

  if (TUchrmatCheckEqual(matrix, result))
  {
//     printf("The representation matrix of represented graph is equal to input matrix.\n");
  }
  else
  {
    printf("Input matrix:\n");
    ASSERT_TU_CALL( TUchrmatPrintDense(stdout, matrix, ' ', true) );
  
    printf("Graph:\n");
    ASSERT_TU_CALL( TUgraphPrint(stdout, graph) );

    printf("Representation matrix:\n");
    ASSERT_TU_CALL( TUchrmatPrintDense(stdout, result, ' ', true) );

    printf("Basis:");
    for (int r = 0; r < matrix->numRows; ++r)
      printf(" %d", basis[r]);
    printf("\n");

    printf("Cobasis:");
    for (int c = 0; c < matrix->numColumns; ++c)
      printf(" %d", cobasis[c]);
    printf("\n");
  }

  ASSERT_TRUE( TUchrmatCheckEqual(matrix, result) );
  ASSERT_TRUE( isGraphic );

  ASSERT_TU_CALL( TUgraphFree(tu, &graph) );
  ASSERT_TU_CALL( TUfreeBlockArray(tu, &basis) );
  ASSERT_TU_CALL( TUfreeBlockArray(tu, &cobasis) );
  ASSERT_TU_CALL( TUchrmatFree(tu, &result) );
}

void testNongraphicMatrix(
  TU* tu,             /**< \ref TU environment. */
  TU_CHRMAT* matrix,  /**< Matrix to be used for testing. */
  int mergeLeafBonds  /**< Leaf bonds of the t-decomposition are merged (1: at the end; 2: after each column). */
)
{
  TU_GRAPH* graph = NULL;
  ASSERT_TU_CALL( TUgraphCreateEmpty(tu, &graph, 0, 0) );
  bool isGraphic;
  TU_CHRMAT* transpose = NULL;
  ASSERT_TU_CALL( TUchrmatTranspose(tu, matrix, &transpose) );
  ASSERT_TU_CALL( testGraphicnessTDecomposition(tu, matrix, transpose, &isGraphic, graph, NULL, NULL, NULL,
    mergeLeafBonds) );
  ASSERT_FALSE( isGraphic );
  ASSERT_TU_CALL( TUchrmatFree(tu, &transpose) );
  ASSERT_TU_CALL( TUgraphFree(tu, &graph) );
}

void testMatrix(
  TU* tu,             /**< \ref TU environment. */
  TU_CHRMAT* matrix,  /**< Matrix to be used for testing. */
  int mergeLeafBonds  /**< Leaf bonds of the t-decomposition are merged (1: at the end; 2: after each column). */
)
{
  TU_GRAPH* graph = NULL;
  ASSERT_TU_CALL( TUgraphCreateEmpty(tu, &graph, 0, 0) );
  bool isGraphic;
  TU_CHRMAT* transpose = NULL;
  ASSERT_TU_CALL( TUchrmatTranspose(tu, matrix, &transpose) );
  ASSERT_TU_CALL( testGraphicnessTDecomposition(tu, matrix, transpose, &isGraphic, graph, NULL, NULL, NULL,
    mergeLeafBonds) );
  ASSERT_TU_CALL( TUchrmatFree(tu, &transpose) );
  ASSERT_TU_CALL( TUgraphFree(tu, &graph) );
}

TEST(Graphic, TypingManyChildrenTerminals)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 10 "
    "1 0 1 1 1 0 1 0 0 0 "
    "0 0 1 0 0 1 1 0 0 0 "
    "0 1 0 1 0 0 1 1 0 0 "
    "0 0 0 0 1 0 1 1 0 0 "
    "0 0 1 0 0 0 1 0 1 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInnerBondWithEdge)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 6 "
    "1 0 0 0 0 1 "
    "0 0 0 1 0 1 "
    "0 0 1 0 1 1 "
    "1 0 1 1 0 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInnerPolygonDoubleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 3 "
    "1 0 1 "
    "0 1 1 "
    "1 1 1 "
    "1 1 0 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingRootPolygonDoubleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 3 "
    "1 1 0 "
    "1 1 1 "
    "1 0 1 "
    "0 1 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingAnyPrimeDegree3)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "1 1 0 1 "
    "0 0 1 0 "
    "1 0 1 1 "
    "0 1 1 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}


TEST(Graphic, TypingAnyPrimeManyPaths)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 7 "
    "0 0 0 1 1 1 1 "
    "0 1 1 1 1 0 0 "
    "0 0 0 1 1 0 1 "
    "0 1 1 1 0 0 1 "
    "0 1 0 0 0 1 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingRootPrimeNoPathsDisjointSingleChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 4 "
    "0 1 1 0 "
    "0 0 1 1 "
    "1 0 1 0 "
    "1 1 0 0 "
    "1 1 0 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}


TEST(Graphic, TypingRootPrimeOnePathSingleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "1 1 0 0 "
    "0 1 1 0 "
    "0 1 0 1 "
    "1 0 1 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingRootPrimeOnePathTwoSingleChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 4 "
    "0 1 1 1 "
    "0 0 1 1 "
    "0 1 1 0 "
    "1 0 1 0 "
    "1 1 0 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}


TEST(Graphic, TypingRootPrimeOnePathDoubleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 6 "
    "0 1 1 1 0 1 "
    "1 1 1 1 1 1 "
    "0 1 0 0 1 0 "
    "0 0 1 0 1 1 "
    "1 1 1 0 1 0 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingRootPrimeTwoPaths)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "3 4 "
    "1 1 1 0 "
    "1 1 0 1 "
    "0 1 1 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInternalPrimeNoPathNoSingleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 7 "
    "1 1 0 0 0 1 0 "
    "1 1 0 1 1 0 0 "
    "1 0 0 1 0 0 0 "
    "0 0 0 0 1 1 0 "
    "1 1 0 1 1 1 0 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInternalPrimeNoPathOneSingleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "1 0 1 1 "
    "0 1 0 1 "
    "0 1 1 0 "
    "1 1 0 0 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInternalPrimeNoPathTwoSingleChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 4 "
    "0 1 1 0 "
    "1 0 1 1 "
    "1 0 0 1 "
    "1 1 0 0 "
    "0 1 0 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}


TEST(Graphic, TypingInternalPrimeOnePathOneSingleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "1 0 1 1 "
    "0 1 0 1 "
    "0 1 1 1 "
    "1 1 0 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );

  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 5 "
    "0 0 0 1 1 "
    "1 0 0 0 1 "
    "1 1 0 1 1 "
    "1 1 0 0 1 "
    "0 1 0 1 0 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );

  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "0 1 1 1 "
    "1 1 0 1 "
    "1 1 0 0 "
    "1 0 1 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInternalPrimeOnePathTwoSingleChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 4 "
    "0 1 1 0 "
    "1 0 1 1 "
    "1 0 0 1 "
    "1 1 0 0 "
    "0 1 0 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );

  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "6 6 "
    "0 1 1 0 1 0 "
    "0 1 1 1 1 1 "
    "1 0 0 1 0 0 "
    "0 0 1 1 0 0 "
    "1 0 1 1 0 0 "
    "1 1 0 1 0 0 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );

  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 4 "
    "0 1 0 1 "
    "1 0 1 1 "
    "1 0 0 1 "
    "0 1 1 0 "
    "1 1 1 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInternalPrimeOnePathNoChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "3 4 "
    "1 1 0 1 "
    "1 1 1 0 "
    "1 0 1 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInternalPrimeOnePathDoubleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 8 "
    "1 1 1 0 0 0 1 1 "
    "0 0 1 1 1 1 0 0 "
    "0 0 0 0 0 1 1 0 "
    "0 1 0 1 0 0 1 0 "
    "0 0 0 0 1 0 1 0 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInternalPrimeTwoPathsNonadjacentParent)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 6 "
    "0 0 0 0 0 1 "
    "1 1 0 1 1 0 "
    "1 1 1 0 0 0 "
    "0 1 0 1 1 0 "
    "0 1 1 0 1 0 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );

  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 6 "
    "0 0 1 0 1 1 "
    "0 1 1 0 0 1 "
    "0 0 0 1 0 0 "
    "0 1 1 1 0 0 "
    "1 0 0 1 1 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInternalPrimeTwoPathOneSingleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 6 "
    "0 0 0 1 0 1 "
    "1 1 1 0 1 1 "
    "1 0 0 1 1 1 "
    "1 1 0 0 0 1 "
    "0 0 0 1 1 0 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInternalPrimeTwoPathTwoSingleChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 6 "
    "0 1 1 1 1 1 "
    "0 0 0 0 1 1 "
    "0 1 1 0 1 1 "
    "0 1 0 0 0 1 "
    "1 0 1 1 0 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, TypingInternalPrimeTwoPathsDoubleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "7 6 "
    "1 1 1 1 1 1 "
    "0 0 1 0 0 0 "
    "1 1 1 1 1 0 "
    "0 1 0 1 1 1 "
    "0 0 0 1 0 1 "
    "1 0 0 0 1 1 "
    "0 1 0 0 0 1 "
  ) );
  testNongraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, RandomMatrix)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  
  srand(0);
  const int numMatrices = 1000;
  const int numRows = 10;
  const int numColumns = 20;
  const double probability = 0.3;

  for (int i = 0; i < numMatrices; ++i)
  {
    TU_CHRMAT* A = NULL;
    TUchrmatCreate(tu, &A, numRows, numColumns, numRows * numColumns);

    A->numNonzeros = 0;
    for (int row = 0; row < numRows; ++row)
    {
      A->rowStarts[row] = A->numNonzeros;
      for (int column = 0; column < numColumns; ++column)
      {
        if ((rand() * 1.0 / RAND_MAX) < probability)
        {
          A->entryColumns[A->numNonzeros] = column;
          A->entryValues[A->numNonzeros] = 1;
          A->numNonzeros++;
        }
      }
    }
    A->rowStarts[numRows] = A->numNonzeros;
    
    TUchrmatPrintDense(stdout, A, '0', false);

    testMatrix(tu, A, 0);

    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootBondNoChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "1 1 "
    "1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootBondTwoSingleChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "3 3 "
    "1 1 0 "
    "0 1 1 "
    "1 0 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootBondTwoSingleChildrenSplit)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "3 4 "
    "1 1 1 0 "
    "1 0 0 1 "
    "0 0 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateInternalBondOneSingleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "3 3 "
    "1 0 1 "
    "1 1 0 "
    "0 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPolygonNoChildrenParent)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "3 3 "
    "1 1 1 "
    "0 1 1 "
    "0 1 0 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPolygonNoChildrenHamiltonianPath)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "3 3 "
    "1 0 0 "
    "1 1 1 "
    "1 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPolygonNoChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "3 2 "
    "1 0 "
    "1 1 "
    "1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}


TEST(Graphic, UpdateRootPolygonOneSingleChildParent)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "3 3 "
    "1 0 1 "
    "1 1 0 "
    "0 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPolygonOneSingleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 3 "
    "1 1 0 "
    "1 1 1 "
    "1 0 1 "
    "0 0 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPolygonTwoSingleChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "1 1 0 0 "
    "1 1 1 1 "
    "0 1 0 1 "
    "0 0 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPolygonTwoSingleChildrenParent)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "1 0 1 1 "
    "0 1 0 1 "
    "1 1 1 0 "
    "0 0 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateLeafPolygon)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "3 3 "
    "1 0 1 "
    "1 1 0 "
    "0 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateInternalPolygonOneSingleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "1 0 0 1 "
    "1 1 0 1 "
    "1 1 1 0 "
    "0 0 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPrimeParentOnly)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "1 0 1 1 "
    "0 1 1 0 "
    "1 1 1 0 "
    "0 1 0 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPrimeParentJoinsPaths)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "1 1 1 1 "
    "1 1 0 1 "
    "0 0 1 0 "
    "0 1 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPrimeNoChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "3 4 "
    "1 0 1 1 "
    "1 1 0 0 "
    "0 1 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPrimeOneSingleChild)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "1 1 0 0 "
    "0 1 0 1 "
    "1 1 1 0 "
    "1 0 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPrimeTwoSingleChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 5 "
    "1 0 0 1 1 "
    "1 1 0 0 1 "
    "0 1 1 0 1 "
    "0 1 0 1 0 "
    "0 0 0 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateRootPrimeTwoParallelSingleChildren)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 5 "
    "1 1 0 1 0 "
    "1 0 0 0 1 "
    "1 0 1 1 0 "
    "0 0 0 1 1 "
    "0 1 1 0 0 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}


TEST(Graphic, UpdateLeafPrime)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "4 4 "
    "1 1 0 1 "
    "1 1 0 0 "
    "1 0 1 0 "
    "0 1 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateInnerPrimeOnePath)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 4 "
    "0 1 1 1 "
    "1 0 1 0 "
    "1 1 0 0 "
    "0 0 1 1 "
    "1 0 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, UpdateInnerPrimeNoPath)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "5 5 "
    "1 1 0 0 1 "
    "0 0 1 1 0 "
    "1 0 1 0 0 "
    "0 1 1 0 1 "
    "0 1 0 1 1 "
  ) );
  testGraphicMatrix(tu, matrix, 2);
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}


TEST(Graphic, UpdateRandomGraph)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  srand(1);
  const int numGraphs = 1000;
  const int numNodes = 6;
  const int numEdges = 10;

  TU_GRAPH_NODE* nodes = NULL;
  ASSERT_TU_CALL( TUallocBlockArray(tu, &nodes, numNodes) );
  for (int i = 0; i < numGraphs; ++i)
  {
    TU_GRAPH* graph = NULL;
    ASSERT_TU_CALL( TUgraphCreateEmpty(tu, &graph, numNodes, numEdges) );
    for (int v = 0; v < numNodes; ++v)
      ASSERT_TU_CALL( TUgraphAddNode(tu, graph, &nodes[v]) );

    for (int e = 0; e < numEdges; ++e)
    {
      int u = (rand() * 1.0 / RAND_MAX) * numNodes;
      int v = (rand() * 1.0 / RAND_MAX) * numNodes;
      ASSERT_TU_CALL( TUgraphAddEdge(tu, graph, nodes[u], nodes[v], NULL) );
    }

//     printf("\n\n\nGraph:\n");
//     ASSERT_TU_CALL( TUgraphPrint(stdout, graph) );

    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( TUconvertGraphToBinaryMatrix(tu, graph, &A, 0, NULL, 0, NULL) );
    
//     FILE* f = fopen("instance.mat", "w");
    ASSERT_TU_CALL( TUchrmatPrintDense(stdout, A, '0', false) );
//     fclose(f);

    ASSERT_TU_CALL( TUgraphFree(tu, &graph) );

    testGraphicMatrix(tu, A, 0);

    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }
  ASSERT_TU_CALL( TUfreeBlockArray(tu, &nodes) );

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}
