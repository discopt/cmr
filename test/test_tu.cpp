#include <gtest/gtest.h>

#include "common.h"

#include <cmr/tu.h>
#include <cmr/separation.h>
#include <cmr/graphic.h>

TEST(TotallyUnimodular, OneSum)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* K_3_3 = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3, "5 4 "
      " 1 1 0  0 "
      " 1 1 1  0 "
      " 1 0 0 -1 "
      " 0 1 1  1 "
      " 0 0 1  1 "
    ) );

    CMR_CHRMAT* K_3_3_dual = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3_dual, "4 5 "
      " 1 1 1 0 0 "
      " 1 1 0 -1 0 "
      " 0 1 0 -1 -1 "
      " 0 0 1 1 1 "
    ) );

    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRoneSum(cmr, K_3_3, K_3_3_dual, &matrix) );

    bool isTU;
    CMR_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, NULL, NULL) );

    ASSERT_CMR_CALL( CMRdecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_TRUE( isTU );
    ASSERT_FALSE( CMRdecHasMatrix(dec) ); /* Default settings should mean that the matrix is not copied. */
    ASSERT_FALSE( CMRdecHasTranspose(dec) ); /* Default settings should mean that the transpose is never computed. */
    ASSERT_EQ( CMRdecIsSum(dec, NULL, NULL), 1 );
    ASSERT_EQ( CMRdecNumChildren(dec), 2 );
    ASSERT_FALSE( CMRdecIsGraphic(CMRdecChild(dec, 0)) );
    ASSERT_TRUE( CMRdecIsCographic(CMRdecChild(dec, 0)) );
    ASSERT_TRUE( CMRdecIsGraphic(CMRdecChild(dec, 1)) );
    ASSERT_FALSE( CMRdecIsCographic(CMRdecChild(dec, 1)) );
    
    ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3_dual) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, SeriesParallelTwoSeparation)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* K_3_3 = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3, "5 4 "
      " 1 1 0  0 "
      " 1 1 1  0 "
      " 1 0 0 -1 "
      " 0 1 1  1 "
      " 0 0 1  1 "
    ) );

    CMR_CHRMAT* K_3_3_dual = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3_dual, "4 5 "
      " 1 1  1 0 0 "
      " 1 1  0 1 0 "
      " 0 1  0 1 1 "
      " 0 0 -1 1 1 "
    ) );

    CMR_CHRMAT* twoSum = NULL;
    ASSERT_CMR_CALL( CMRtwoSum(cmr, K_3_3, K_3_3_dual, CMRrowToElement(1), CMRcolumnToElement(1), &twoSum) );

    size_t rowPermutations[] = { 4, 6, 5, 7, 0, 1, 2, 3 };
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRchrmatPermute(cmr, twoSum, rowPermutations, NULL, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &twoSum) );

    CMRchrmatPrintDense(cmr, K_3_3, stdout, '0', true);
    CMRchrmatPrintDense(cmr, K_3_3_dual, stdout, '0', true);
    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, NULL, NULL) );

    ASSERT_CMR_CALL( CMRdecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_TRUE( isTU );
    ASSERT_FALSE( CMRdecHasMatrix(dec) ); /* Default settings should mean that the matrix is not copied. */
    ASSERT_TRUE( CMRdecHasTranspose(dec) ); /* As we test for graphicness, the transpose is constructed. */
    ASSERT_EQ( CMRdecIsSum(dec, NULL, NULL), 2 );
    ASSERT_EQ( CMRdecNumChildren(dec), 2 );
    ASSERT_FALSE( CMRdecIsGraphic(CMRdecChild(dec, 0)) );
    ASSERT_TRUE( CMRdecIsCographic(CMRdecChild(dec, 0)) );
    ASSERT_TRUE( CMRdecIsGraphic(CMRdecChild(dec, 1)) );
    ASSERT_FALSE( CMRdecIsCographic(CMRdecChild(dec, 1)) );

    ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3_dual) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, NestedMinorSearchTwoSeparation)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* K_3_3 = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3, "5 4 "
      " 1 1 0  0 "
      " 1 1 1  0 "
      " 1 0 0 -1 "
      " 0 1 1  1 "
      " 0 0 1  1 "
    ) );

    CMR_CHRMAT* K_3_3_dual = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3_dual, "4 5 "
      " 1 1  1 0 0 "
      " 1 1  0 1 0 "
      " 0 1  0 1 1 "
      " 0 0 -1 1 1 "
    ) );

    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRtwoSum(cmr, K_3_3, K_3_3_dual, CMRrowToElement(1), CMRcolumnToElement(1), &matrix) );

    bool isTU;
    CMR_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, NULL, NULL) );

    ASSERT_CMR_CALL( CMRdecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_TRUE( isTU );
    ASSERT_FALSE( CMRdecHasMatrix(dec) ); /* Default settings should mean that the matrix is not copied. */
    ASSERT_TRUE( CMRdecHasTranspose(dec) ); /* As we test for graphicness, the transpose is constructed. */
    ASSERT_EQ( CMRdecIsSum(dec, NULL, NULL), 2 );
    ASSERT_EQ( CMRdecNumChildren(dec), 2 );
    ASSERT_TRUE( CMRdecIsGraphic(CMRdecChild(dec, 0)) );
    ASSERT_FALSE( CMRdecIsCographic(CMRdecChild(dec, 0)) );
    ASSERT_FALSE( CMRdecIsGraphic(CMRdecChild(dec, 1)) );
    ASSERT_TRUE( CMRdecIsCographic(CMRdecChild(dec, 1)) );
    
    ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3_dual) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, NestedMinorPivotsOneRowOneColumn)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 6 "
    " 1 0 1 0 0 0 "
    " 1 1 0 0 0 1 "
    " 0 1 1 0 0 0 "
    " 0 0 0 0 1 1 "
    " 0 0 0 1 1 0 "
    " 0 1 1 1 0 0 "
  ) );

  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

  bool isTU;
  CMR_DEC* dec = NULL;
  CMR_TU_PARAMETERS params;
  ASSERT_CMR_CALL( CMRparamsTotalUnimodularityInit(&params) );
  params.regular.fastGraphicness = false;
  ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, &params, NULL) );

  ASSERT_CMR_CALL( CMRdecPrint(cmr, dec, stdout, 0, true, true, true) );
  
  ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, NestedMinorPivotsTwoRowsOneColumn)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "7 6 "
    " 1 0 1 0 0 0 "
    " 1 1 0 0 0 0 "
    " 0 1 1 0 0 0 "
    " 0 0 1 0 0 1 "
    " 0 0 0 0 1 1 "
    " 0 0 0 1 1 0 "
    " 0 1 1 1 0 0 "
  ) );

  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

  bool isTU;
  CMR_DEC* dec = NULL;
  CMR_TU_PARAMETERS params;
  ASSERT_CMR_CALL( CMRparamsTotalUnimodularityInit(&params) );
  params.regular.fastGraphicness = false;
  ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, &params, NULL) );

  ASSERT_CMR_CALL( CMRdecPrint(cmr, dec, stdout, 0, true, true, true) );
  
  ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, NestedMinorPivotsOneRowTwoColumns)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 7 "
    " 1 0 1 0 0 0 0 "
    " 1 1 0 0 0 0 1 "
    " 0 1 1 1 0 0 1 "
    " 0 0 0 0 0 1 1 "
    " 0 0 0 0 1 1 0 "
    " 0 0 0 1 1 0 0 "
  ) );

  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

  bool isTU;
  CMR_DEC* dec = NULL;
  CMR_TU_PARAMETERS params;
  ASSERT_CMR_CALL( CMRparamsTotalUnimodularityInit(&params) );
  params.regular.fastGraphicness = false;
  ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, &params, NULL) );

  ASSERT_CMR_CALL( CMRdecPrint(cmr, dec, stdout, 0, true, true, true) );
  
  ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, NestedMinorPivotsTwoSeparation)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "11 11 "
    " 1 0 1 0 0 0 0 0 0 0 0 "
    " 1 1 0 0 0 1 0 0 0 0 0 "
    " 0 1 1 0 0 0 0 0 0 0 0 "
    " 0 0 0 0 1 1 0 0 0 0 0 "
    " 0 0 0 1 1 0 0 0 0 0 1 "
    " 0 1 1 1 0 0 0 0 0 1 0 "
    " 0 0 0 0 0 0 0 1 1 0 0 "
    " 0 0 0 0 0 0 1 1 0 0 0 "
    " 0 1 1 0 0 0 1 0 0 0 0 "
    " 0 1 1 0 0 0 0 0 1 0 0 "
    " 0 0 0 0 0 0 0 0 0 1 1 "
  ) );

  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

  bool isTU;
  CMR_DEC* dec = NULL;
  CMR_TU_PARAMETERS params;
  ASSERT_CMR_CALL( CMRparamsTotalUnimodularityInit(&params) );
  params.regular.fastGraphicness = false;
  ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, &params, NULL) );

  ASSERT_CMR_CALL( CMRdecPrint(cmr, dec, stdout, 0, true, true, true) );
  
  ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

/**
 * \brief Tests a matrix that is graphic.
 */

static
void testSequenceGraphicness(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix, /**< Matrix to test. */
  bool knowNetwork    /**< Whether the matrix is graphic. */
)
{
  printf("Testing matrix for graphicness via nested minor sequence:\n");
  ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );

  bool isTU;
  CMR_DEC* dec = NULL;
  CMR_TU_PARAMETERS params;
  ASSERT_CMR_CALL( CMRparamsTotalUnimodularityInit(&params) );
  params.regular.fastGraphicness = false;
  ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, &params, NULL) );

  if (knowNetwork)
  {
    ASSERT_TRUE( isTU );
    ASSERT_TRUE( CMRdecIsGraphic(dec) );
    CMR_CHRMAT* networkMatrix = NULL;
    bool isForest;
    ASSERT_CMR_CALL( CMRcomputeNetworkMatrix(cmr, CMRdecGraph(dec), &networkMatrix, NULL, CMRdecGraphArcsReversed(dec),
      CMRdecGraphSizeForest(dec), CMRdecGraphForest(dec), CMRdecGraphSizeCoforest(dec), CMRdecGraphCoforest(dec), &isForest) );

// TODO: Decomposition graph is currently undirected.
//     if (!CMRchrmatCheckEqual(matrix, networkMatrix))
//     {
//       printf("Computed digraph with different representation matrix:\n");
//       CMRchrmatPrintDense(cmr, networkMatrix, stdout, '0', true);
//     }
//     ASSERT_TRUE( CMRchrmatCheckEqual(matrix, networkMatrix) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &networkMatrix) );
    ASSERT_TRUE(isForest);
  }
  else
  {
    ASSERT_FALSE( CMRdecIsGraphic(dec) );
  }

  ASSERT_CMR_CALL( CMRdecPrint(cmr, dec, stdout, 0, true, true, true) );

  ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
}

TEST(TotallyUnimodular, SequenceGraphicnessWheel)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 6 "
      "1 0 0 0 0 1 "
      "1 1 0 0 0 0 "
      "0 1 1 0 0 0 "
      "0 0 1 1 0 0 "
      "0 0 0 1 1 0 "
      "0 0 0 0 1 1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
      "1  0 1 0 "
      "1 -1 0 0 "
      "0  1 1 1 "
      "0  0 1 1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
      "1 1 1 0 "
      "1 1 0 0 "
      "0 1 1 1 "
      "0 0 1 1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }
  
  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
      "1 0 1 0 "
      "1 1 1 0 "
      "0 1 1 1 "
      "0 0 1 1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
      "1 0 1 0 "
      "1 1 0 0 "
      "1 1 1 1 "
      "0 0 1 1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, SequenceGraphicnessOneRowOneColumn)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "7 7 "
      "1 0 0 1   0   0   0 "
      "1 1 0 0   0   0   0 "
      "0 1 1 0   0   0   0 "
      "0 0 1 1   1   1   1 "
      "                    "
      "0 0 0 1   1   1   1 "
      "                    "
      "0 0 1 1   1   0   0 "
      "                    "
      "0 0 0 0   1   0   1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, SequenceGraphicnessTwoRowsOneColumn)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "8 6 "
      "1 0 0 1   0   0 "
      "1 1 0 0   0   0 "
      "0 1 1 0   0   0 "
      "0 0 1 1   0   0 "
      "                "
      "0 0 0 1   1   0 "
      "0 0 1 1   1   0 "
      "                "
      "0 0 0 1   0   1 "
      "0 0 0 0  -1   1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, SequenceGraphicnessOneRowTwoColumns)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 8 "
      "1 0 0 1   0 0   0 0 "
      "1 1 0 0   0 0   0 0 "
      "0 1 1 0  -1 0  -1 0 "
      "0 0 1 1   0 1   0 1 "
      "                    "
      "0 0 0 0   1 1   1 1 "
      "                    "
      "0 0 0 0   0 0   1 1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, SequenceGraphicnessOneColumn)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 6 "
      "1 0 0 1   0   0 "
      "1 1 0 0   0   0 "
      "0 1 1 0   0   0 "
      "0 0 1 1   0   1 "
      "                "
      "0 0 0 1   1   0 "
      "0 0 1 1   1   1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, SequenceGraphicnessOneRow)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 5 "
      "1 1  0   0   0 "
      "0 1  1   0   1 "
      "1 0 -1   1   0 "
      "              "
      "1 0  0   1   0 "
      "              "
      "0 1  0   0   1 "
      "              "
      "0 1  1   0   0 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
      "1 1  0   0  1 " 
      "0 1  1   1  0 "
      "1 0 -1   0  0 "
      "             "
      "0 1  0   1  1 "
      "             "
      "0 0  0   1  1 "
    ) );
    testSequenceGraphicness(cmr, matrix, false);
  }

  {
    /* Runs into addition of a single row, at least 1 articulation point, but none begin part of all fundamental cycles
     * induced by 1-edges. */
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "10 9 "
      "1 1  0  0 0 0 0 0 0 "
      "0 1  1  1 0 0 0 0 0 "
      "0 1  1  0 0 0 0 0 0 "
      "1 0 -1  0 0 0 0 0 0 "
      "0 0  1  1 0 0 0 0 1 "
      "1 0 -1 -1 1 0 0 0 0 "
      "0 0  0  1 0 1 0 0 0 "
      "0 0  0  0 0 1 1 1 0 "
      "0 0  0  0 0 1 0 0 0 "
      "0 0  0  0 0 0 1 0 0 "
    ) );
    testSequenceGraphicness(cmr, matrix, false);
  }
  
  {
    /* Runs into addition of a single row with a unique 1 articulation point, but for which the auxiliary graph is not
     * bipartite. */
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "8 10 "
      "1 0 1 0 1 0 0  0 0  0 "
      "1 1 1 1 1 0 0  0 0  0 "
      "0 0 0 0 1 1 0  0 0  0 "
      "0 1 0 1 1 0 0  1 1  0 "
      "0 1 1 0 1 1 1  0 0  0 "
      "0 0 0 0 0 0 1 -1 0  0 "
      "0 0 0 1 0 0 0  0 0  1 "
      "0 0 1 0 0 0 0  0 0 -1 "
    ) );
    testSequenceGraphicness(cmr, matrix, false);
  }

  {
    /* Runs into addition of a single row with a unique 1 articulation point and bipartite auxiliary graph, and a
     * 1-edge adjacent to the split node. */
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
      "1 1 0  0 1 "
      "1 1 1  0 1 "
      "0 1 1  1 0 "
      "1 0 0 -1 0 "
      "0 1 1  1 1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, R10)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
      "1 0 0 1 1 "
      "1 1 0 0 1 "
      "0 1 1 0 1 "
      "0 0 1 1 1 "
      "1 1 1 1 1 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, NULL, NULL) );
    ASSERT_TRUE( CMRdecIsRegular(dec) );
    ASSERT_EQ( CMRdecNumChildren(dec), 0 );
    ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
      "1 1  0  0 1 "
      "1 1 -1  0 0 "
      "0 1 -1 -1 0 "
      "0 0  1  1 1 "
      "1 0  0  1 1 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, NULL, NULL) );
    ASSERT_TRUE( dec );
    ASSERT_TRUE( CMRdecIsRegular(dec) );
    ASSERT_EQ( CMRdecNumChildren(dec), 0 );
    ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


/**
 * \brief Tests a 3-separable matrix for regularity.
 */

static
void testEnumerate(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix, /**< Matrix to test. */
  bool knowRegular    /**< Whether the matrix is regular. */
)
{
  printf("Testing matrix for regularity with a 3-separation:\n");
  ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );

  bool isTU;
  CMR_DEC* dec = NULL;
  CMR_TU_PARAMETERS params;
  ASSERT_CMR_CALL( CMRparamsTotalUnimodularityInit(&params) );
  params.regular.fastGraphicness = false;
  ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, &params, NULL) );

  if (knowRegular)
  {
    ASSERT_TRUE( isTU );
  }
  else
  {
    ASSERT_FALSE( isTU );
  }

  ASSERT_CMR_CALL( CMRdecPrint(cmr, dec, stdout, 0, true, true, true) );

  ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
}

TEST(TotallyUnimodular, EnumerateRanksZeroTwo)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
      "1 0 1 1 "
      "1 1 1 0 "
      "0 0 1 1 "
      "0 1 1 1 "
    ) );
    testEnumerate(cmr, matrix, false);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, EnumerateRanksOneOne)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
      "1 1 1 1 "
      "0 1 1 1 "
      "1 0 0 1 "
      "1 1 0 0 "
    ) );
    testEnumerate(cmr, matrix, false);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, EnumerateRanksTwoZero)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
      "1 1 0 1 0 "
      "1 1 1 0 1 "
      "1 1 1 0 0 "
      "1 1 1 1 1 "
      "1 0 0 0 1 "
    ) );
    testEnumerate(cmr, matrix, false);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, R12)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 6 "
      "1  0 1  1 0 0 "
      "0  1 1  1 0 0 "
      "1  0 1  0 1 1 "
      "0 -1 0 -1 1 1 "
      "1  0 1  0 1 0 "
      "0 -1 0 -1 0 1 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, NULL, NULL, NULL) );
    ASSERT_TRUE( CMRdecIsRegular(dec) );
    ASSERT_EQ( CMRdecNumChildren(dec), 2 );
    size_t graphicChildren = (CMRdecIsGraphic(CMRdecChild(dec, 0)) ? 2 : 0)
      + (CMRdecIsGraphic(CMRdecChild(dec, 1)) ? 1 : 0);
    size_t cographicChildren = (CMRdecIsCographic(CMRdecChild(dec, 0)) ? 2 : 0)
      + (CMRdecIsCographic(CMRdecChild(dec, 1)) ? 1 : 0);
    ASSERT_EQ( CMRdecNumChildren(CMRdecChild(dec, 0)), 0 );
    ASSERT_EQ( CMRdecNumChildren(CMRdecChild(dec, 1)), 0 );
    ASSERT_EQ( graphicChildren + cographicChildren, 3 );
    ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TotallyUnimodular, ForbiddenSubmatrix)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "14 14 "
      "1 1 1 0 1 0 1 0  1 1 1 1 1 1 "
      "1 0 1 0 1 0 1 0  1 1 1 1 1 0 "
      "0 1 1 0 0 0 0 0  0 0 0 0 0 0 "
      "0 1 1 1 0 0 0 0  0 0 0 0 0 0 "
      "0 1 1 1 1 0 0 0  0 0 0 0 0 0 "
      "0 1 1 1 1 1 0 0  0 0 0 0 0 0 "
      "0 1 1 1 1 1 1 0  0 0 0 0 0 0 "
      "0 1 1 1 1 1 1 1  0 0 0 0 0 0 "
      "0 1 1 1 1 1 1 1  1 0 0 0 0 0 "
      "0 0 0 0 0 0 0 0  1 1 0 0 0 0 "
      "0 0 0 0 0 0 0 0  0 1 1 0 0 0 "
      "0 0 0 0 0 0 0 0  0 0 1 1 0 0 "
      "0 0 0 0 0 0 0 0  0 0 0 1 1 0 "
      "0 0 0 0 0 0 0 0  0 0 0 0 1 1 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_DEC* dec = NULL;
    CMR_SUBMAT* forbiddenSubmatrix = NULL;
    ASSERT_CMR_CALL( CMRtestTotalUnimodularity(cmr, matrix, &isTU, &dec, &forbiddenSubmatrix, NULL, NULL) );
    ASSERT_FALSE( CMRdecIsRegular(dec) );
    ASSERT_EQ( forbiddenSubmatrix->numRows, 8 );
    ASSERT_EQ( forbiddenSubmatrix->numColumns, 8 );
    // TODO: Compute determinant once implemented.
    ASSERT_CMR_CALL( CMRsubmatFree(cmr, &forbiddenSubmatrix) );
    ASSERT_CMR_CALL( CMRdecFree(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
