#include <gtest/gtest.h>

#include <stdlib.h>

#include "common.h"

#include <cmr/graphic.h>
#include <cmr/network.h>

void testBinaryGraphicMatrix(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix /**< Matrix to be used for testing. */
)
{
  bool isGraphic;
  CMR_GRAPH* graph = NULL;
  CMR_GRAPH_EDGE* basis = NULL;
  CMR_GRAPH_EDGE* cobasis = NULL;
  CMR_CHRMAT* transpose = NULL;
  ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

  ASSERT_CMR_CALL( CMRtestBinaryGraphic(cmr, transpose, &isGraphic, &graph, &basis, &cobasis, NULL) );

  ASSERT_TRUE( isGraphic );
  ASSERT_TRUE( basis );
  ASSERT_TRUE( cobasis );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );

  CMR_CHRMAT* result = NULL;
  bool isCorrectBasis;
  ASSERT_CMR_CALL( CMRcomputeGraphBinaryRepresentationMatrix(cmr, graph, &result, NULL, matrix->numRows, basis,
    matrix->numColumns, cobasis, &isCorrectBasis) );
  ASSERT_TRUE( isCorrectBasis );
  ASSERT_TRUE( result );

  if (!CMRchrmatCheckEqual(matrix, result))
  {
    printf("Input matrix:\n");
    ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, stdout, matrix, ' ', true) );
  
    printf("Graph:\n");
    ASSERT_CMR_CALL( CMRgraphPrint(stdout, graph) );

    printf("Representation matrix:\n");
    ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, stdout, result, ' ', true) );

    printf("Basis:");
    for (int r = 0; r < matrix->numRows; ++r)
      printf(" %d", basis[r]);
    printf("\n");

    printf("Cobasis:");
    for (int c = 0; c < matrix->numColumns; ++c)
      printf(" %d", cobasis[c]);
    printf("\n");
  }

  ASSERT_TRUE( CMRchrmatCheckEqual(matrix, result) );
  ASSERT_TRUE( isGraphic );

  ASSERT_CMR_CALL( CMRgraphFree(cmr, &graph) );
  ASSERT_CMR_CALL( CMRfreeBlockArray(cmr, &basis) );
  ASSERT_CMR_CALL( CMRfreeBlockArray(cmr, &cobasis) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &result) );
}

void testBinaryNongraphicMatrix(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix /**< Matrix to be used for testing. */
)
{
  bool isGraphic;
  CMR_CHRMAT* transpose = NULL;
  ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );
  ASSERT_CMR_CALL( CMRtestBinaryGraphic(cmr, transpose, &isGraphic, NULL, NULL, NULL, NULL) );
  ASSERT_FALSE( isGraphic );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
}

void testBinaryMatrix(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix /**< Matrix to be used for testing. */
)
{
  CMR_GRAPH* graph = NULL;
  bool isGraphic;
  CMR_CHRMAT* transpose = NULL;
  ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );
  ASSERT_CMR_CALL( CMRtestBinaryGraphic(cmr, transpose, &isGraphic, &graph, NULL, NULL, NULL) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
  if (graph)
    ASSERT_CMR_CALL( CMRgraphFree(cmr, &graph) );
}

TEST(Graphic, TypingManyChildrenTerminals)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 10 "
    "1 0 1 1 1 0 1 0 0 0 "
    "0 0 1 0 0 1 1 0 0 0 "
    "0 1 0 1 0 0 1 1 0 0 "
    "0 0 0 0 1 0 1 1 0 0 "
    "0 0 1 0 0 0 1 0 1 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerParallelWithEdge)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 6 "
    "1 0 0 0 0 1 "
    "0 0 0 1 0 1 "
    "0 0 1 0 1 1 "
    "1 0 1 1 0 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerSeriesDoubleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 3 "
    "1 0 1 "
    "0 1 1 "
    "1 1 1 "
    "1 1 0 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingRootSeriesDoubleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 3 "
    "1 1 0 "
    "1 1 1 "
    "1 0 1 "
    "0 1 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingAnyRigidDegree3)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 1 0 1 "
    "0 0 1 0 "
    "1 0 1 1 "
    "0 1 1 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingAnyRigidManyPaths)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 7 "
    "0 0 0 1 1 1 1 "
    "0 1 1 1 1 0 0 "
    "0 0 0 1 1 0 1 "
    "0 1 1 1 0 0 1 "
    "0 1 0 0 0 1 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingRootRigidNoPathsDisjointSingleChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 4 "
    "0 1 1 0 "
    "0 0 1 1 "
    "1 0 1 0 "
    "1 1 0 0 "
    "1 1 0 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


TEST(Graphic, TypingRootRigidOnePathSingleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 1 0 0 "
    "0 1 1 0 "
    "0 1 0 1 "
    "1 0 1 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingRootRigidOnePathTwoSingleChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 4 "
    "0 1 1 1 "
    "0 0 1 1 "
    "0 1 1 0 "
    "1 0 1 0 "
    "1 1 0 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


TEST(Graphic, TypingRootRigidOnePathDoubleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 6 "
    "0 1 1 1 0 1 "
    "1 1 1 1 1 1 "
    "0 1 0 0 1 0 "
    "0 0 1 0 1 1 "
    "1 1 1 0 1 0 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingRootRigidTwoPaths)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 4 "
    "1 1 1 0 "
    "1 1 0 1 "
    "0 1 1 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerRigidNoPathNoSingleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 7 "
    "1 1 0 0 0 1 0 "
    "1 1 0 1 1 0 0 "
    "1 0 0 1 0 0 0 "
    "0 0 0 0 1 1 0 "
    "1 1 0 1 1 1 0 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerRigidNoPathOneSingleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 0 1 1 "
    "0 1 0 1 "
    "0 1 1 0 "
    "1 1 0 0 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerRigidNoPathTwoSingleChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 4 "
    "0 1 1 0 "
    "1 0 1 1 "
    "1 0 0 1 "
    "1 1 0 0 "
    "0 1 0 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


TEST(Graphic, TypingInnerRigidOnePathOneSingleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 0 1 1 "
    "0 1 0 1 "
    "0 1 1 1 "
    "1 1 0 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
    "0 0 0 1 1 "
    "1 0 0 0 1 "
    "1 1 0 1 1 "
    "1 1 0 0 1 "
    "0 1 0 1 0 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "0 1 1 1 "
    "1 1 0 1 "
    "1 1 0 0 "
    "1 0 1 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerRigidOnePathTwoSingleChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 4 "
    "0 1 1 0 "
    "1 0 1 1 "
    "1 0 0 1 "
    "1 1 0 0 "
    "0 1 0 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 6 "
    "0 1 1 0 1 0 "
    "0 1 1 1 1 1 "
    "1 0 0 1 0 0 "
    "0 0 1 1 0 0 "
    "1 0 1 1 0 0 "
    "1 1 0 1 0 0 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 4 "
    "0 1 0 1 "
    "1 0 1 1 "
    "1 0 0 1 "
    "0 1 1 0 "
    "1 1 1 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerRigidOnePathNoChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 4 "
    "1 1 0 1 "
    "1 1 1 0 "
    "1 0 1 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerRigidOnePathDoubleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 8 "
    "1 1 1 0 0 0 1 1 "
    "0 0 1 1 1 1 0 0 "
    "0 0 0 0 0 1 1 0 "
    "0 1 0 1 0 0 1 0 "
    "0 0 0 0 1 0 1 0 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerRigidTwoPathsNonadjacentParent)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 6 "
    "0 0 0 0 0 1 "
    "1 1 0 1 1 0 "
    "1 1 1 0 0 0 "
    "0 1 0 1 1 0 "
    "0 1 1 0 1 0 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 6 "
    "0 0 1 0 1 1 "
    "0 1 1 0 0 1 "
    "0 0 0 1 0 0 "
    "0 1 1 1 0 0 "
    "1 0 0 1 1 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerRigidTwoPathOneSingleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 6 "
    "0 0 0 1 0 1 "
    "1 1 1 0 1 1 "
    "1 0 0 1 1 1 "
    "1 1 0 0 0 1 "
    "0 0 0 1 1 0 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerRigidTwoPathTwoSingleChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 6 "
    "0 1 1 1 1 1 "
    "0 0 0 0 1 1 "
    "0 1 1 0 1 1 "
    "0 1 0 0 0 1 "
    "1 0 1 1 0 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, TypingInnerRigidTwoPathsDoubleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "7 6 "
    "1 1 1 1 1 1 "
    "0 0 1 0 0 0 "
    "1 1 1 1 1 0 "
    "0 1 0 1 1 1 "
    "0 0 0 1 0 1 "
    "1 0 0 0 1 1 "
    "0 1 0 0 0 1 "
  ) );
  testBinaryNongraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, RandomMatrix)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  
  srand(0);
  const int numMatrices = 1000;
  const int numRows = 100;
  const int numColumns = 200;
  const double probability = 0.3;

  for (int i = 0; i < numMatrices; ++i)
  {
    CMR_CHRMAT* A = NULL;
    CMRchrmatCreate(cmr, &A, numRows, numColumns, numRows * numColumns);

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
    
    /* CMRchrmatPrintDense(stdout, A, '0', false); */

    testBinaryMatrix(cmr, A);

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &A) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootParallelNoChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "1 1 "
    "1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootParallelTwoSingleChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 3 "
    "1 1 0 "
    "0 1 1 "
    "1 0 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootParallelTwoSingleChildrenSplit)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 4 "
    "1 1 1 0 "
    "1 0 0 1 "
    "0 0 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateInnerParallelOneSingleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 3 "
    "1 0 1 "
    "1 1 0 "
    "0 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootSeriesNoChildrenParent)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 3 "
    "1 1 1 "
    "0 1 1 "
    "0 1 0 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootSeriesNoChildrenHamiltonianPath)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 3 "
    "1 0 0 "
    "1 1 1 "
    "1 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootSeriesNoChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 2 "
    "1 0 "
    "1 1 "
    "1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


TEST(Graphic, UpdateRootSeriesOneSingleChildParent)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 3 "
    "1 0 1 "
    "1 1 0 "
    "0 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootSeriesOneSingleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 3 "
    "1 1 0 "
    "1 1 1 "
    "1 0 1 "
    "0 0 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootSeriesTwoSingleChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 1 0 0 "
    "1 1 1 1 "
    "0 1 0 1 "
    "0 0 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootSeriesTwoSingleChildrenParent)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 0 1 1 "
    "0 1 0 1 "
    "1 1 1 0 "
    "0 0 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateLeafSeries)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 3 "
    "1 0 1 "
    "1 1 0 "
    "0 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateInnerSeriesOneSingleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 0 0 1 "
    "1 1 0 1 "
    "1 1 1 0 "
    "0 0 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootRigidParentOnly)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 0 1 1 "
    "0 1 1 0 "
    "1 1 1 0 "
    "0 1 0 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootRigidParentJoinsPaths)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 1 1 1 "
    "1 1 0 1 "
    "0 0 1 0 "
    "0 1 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootRigidNoChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 4 "
    "1 0 1 1 "
    "1 1 0 0 "
    "0 1 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootRigidOneSingleChild)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 1 0 0 "
    "0 1 0 1 "
    "1 1 1 0 "
    "1 0 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootRigidTwoSingleChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
    "1 0 0 1 1 "
    "1 1 0 0 1 "
    "0 1 1 0 1 "
    "0 1 0 1 0 "
    "0 0 0 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRootRigidTwoParallelSingleChildren)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
    "1 1 0 1 0 "
    "1 0 0 0 1 "
    "1 0 1 1 0 "
    "0 0 0 1 1 "
    "0 1 1 0 0 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


TEST(Graphic, UpdateLeafRigid)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 1 0 1 "
    "1 1 0 0 "
    "1 0 1 0 "
    "0 1 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateInnerRigidOnePath)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 4 "
    "0 1 1 1 "
    "1 0 1 0 "
    "1 1 0 0 "
    "0 0 1 1 "
    "1 0 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateInnerRigidNoPath)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
    "1 1 0 0 1 "
    "0 0 1 1 0 "
    "1 0 1 0 0 "
    "0 1 1 0 1 "
    "0 1 0 1 1 "
  ) );
  testBinaryGraphicMatrix(cmr, matrix);
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, UpdateRandomGraph)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  srand(1);
  const int numGraphs = 1000;
  const int numNodes = 21;
  const int numEdges = 40;

  CMR_GRAPH_NODE* nodes = NULL;
  ASSERT_CMR_CALL( CMRallocBlockArray(cmr, &nodes, numNodes) );
  for (int i = 0; i < numGraphs; ++i)
  {
    CMR_GRAPH* graph = NULL;
    ASSERT_CMR_CALL( CMRgraphCreateEmpty(cmr, &graph, numNodes, numEdges) );
    for (int v = 0; v < numNodes; ++v)
      ASSERT_CMR_CALL( CMRgraphAddNode(cmr, graph, &nodes[v]) );

    for (int e = 0; e < numEdges; ++e)
    {
      int u = (rand() * 1.0 / RAND_MAX) * numNodes;
      int v = (rand() * 1.0 / RAND_MAX) * numNodes;
      ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, nodes[u], nodes[v], NULL) );
    }

    CMR_CHRMAT* A = NULL;
    ASSERT_CMR_CALL( CMRcomputeGraphBinaryRepresentationMatrix(cmr, graph, &A, NULL, 0, NULL, 0, NULL, NULL) );

    ASSERT_CMR_CALL( CMRgraphFree(cmr, &graph) );

    testBinaryGraphicMatrix(cmr, A);

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &A) );
  }
  ASSERT_CMR_CALL( CMRfreeBlockArray(cmr, &nodes) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Graphic, RepresentationMatrix)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_GRAPH* graph = NULL;
  ASSERT_CMR_CALL( CMRgraphCreateEmpty(cmr, &graph, 0, 0) );

  /**
   * Spanning forest: 
   *
   *        c2                h7 --  
   *        ^                 ^   \
   *        |                 \   /
   * a0 --> b1 <-- d3 --> e4   ---
   *
   *     -->
   *  f5 --> g6
   *     -->
   */
  
  CMR_GRAPH_NODE a, b, c, d, e, f, g, h;
  ASSERT_CMR_CALL( CMRgraphAddNode(cmr, graph, &a) );
  ASSERT_CMR_CALL( CMRgraphAddNode(cmr, graph, &b) );
  ASSERT_CMR_CALL( CMRgraphAddNode(cmr, graph, &c) );
  ASSERT_CMR_CALL( CMRgraphAddNode(cmr, graph, &d) );
  ASSERT_CMR_CALL( CMRgraphAddNode(cmr, graph, &e) );
  ASSERT_CMR_CALL( CMRgraphAddNode(cmr, graph, &f) );
  ASSERT_CMR_CALL( CMRgraphAddNode(cmr, graph, &g) );
  ASSERT_CMR_CALL( CMRgraphAddNode(cmr, graph, &h) );
  
  CMR_GRAPH_EDGE basis[5];
  ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, a, b, &basis[0]) );
  ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, b, c, &basis[1]) );
  ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, d, b, &basis[2]) );
  ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, d, e, &basis[3]) );
  ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, f, g, &basis[4]) );

  CMR_GRAPH_EDGE cobasis[6];
  ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, c, a, &cobasis[0]) );
  ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, a, e, &cobasis[1]) );
  ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, e, c, &cobasis[2]) );
  ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, f, g, &cobasis[3]) );
  ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, g, f, &cobasis[4]) );
  ASSERT_CMR_CALL( CMRgraphAddEdge(cmr, graph, h, h, &cobasis[5]) );

  bool isCorrectBasis = false;
  CMR_CHRMAT* matrix = NULL;
  CMR_CHRMAT* check = NULL;

  /* Check binary representation matrix. */
  ASSERT_CMR_CALL( CMRcomputeGraphBinaryRepresentationMatrix(cmr, graph, &matrix, NULL, 5, basis, 6, cobasis, &isCorrectBasis) );
  ASSERT_TRUE( isCorrectBasis );
  stringToCharMatrix(cmr, &check, "5 6 "
    "1 1 0 0 0 0 "
    "1 0 1 0 0 0 "
    "0 1 1 0 0 0 "
    "0 1 1 0 0 0 "
    "0 0 0 1 1 0 "
  );
  ASSERT_TRUE( CMRchrmatCheckEqual(matrix, check) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  /* Check ternary representation matrix. */
  bool edgesReversed[] = { false, false, false, false, false, false, false, false, false, false, false };
  ASSERT_CMR_CALL( CMRcomputeNetworkMatrix(cmr, graph, &matrix, NULL, edgesReversed, 5, basis, 6, cobasis,
    &isCorrectBasis) );
  ASSERT_TRUE( isCorrectBasis );
  stringToCharMatrix(cmr, &check, "5 6 "
    "-1  1  0  0  0  0 "
    "-1  0  1  0  0  0 "
    "0  -1  1  0  0  0 "
    "0   1 -1  0  0  0 "
    "0   0  0  1 -1  0 "
  );
  ASSERT_TRUE( CMRchrmatCheckEqual(matrix, check) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRgraphFree(cmr, &graph) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
