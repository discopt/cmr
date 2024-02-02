#include <gtest/gtest.h>

#include "common.h"

#include <cmr/tu.h>
#include <cmr/separation.h>
#include <cmr/graphic.h>

TEST(TU, OneSum)
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
    CMR_MATROID_DEC* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.regular.planarityCheck = true;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_TRUE( isTU );
    ASSERT_FALSE( CMRmatroiddecHasTranspose(dec) ); /* Default settings should mean that the transpose is never computed. */
    ASSERT_EQ( CMRmatroiddecType(dec), CMR_MATROID_DEC_TYPE_ONE_SUM );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 2UL );
    ASSERT_LT( CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 0))
      * CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 1)), 0);
    ASSERT_LT( CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 0))
      * CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 1)), 0);

    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3_dual) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, SeriesParallelTwoSeparation)
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
    ASSERT_CMR_CALL( CMRtwoSum(cmr, K_3_3, K_3_3_dual, CMRrowToElement(1), CMRcolumnToElement(1), 3, &twoSum) );

    size_t rowPermutations[] = { 4, 6, 5, 7, 0, 1, 2, 3 };
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRchrmatPermute(cmr, twoSum, rowPermutations, NULL, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &twoSum) );

    CMRchrmatPrintDense(cmr, K_3_3, stdout, '0', true);
    CMRchrmatPrintDense(cmr, K_3_3_dual, stdout, '0', true);
    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_MATROID_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_TRUE( isTU );
    ASSERT_TRUE( CMRmatroiddecHasTranspose(dec) ); /* As we test for graphicness, the transpose is constructed. */
    ASSERT_EQ( CMRmatroiddecType(dec), CMR_MATROID_DEC_TYPE_TWO_SUM );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 2UL );
    ASSERT_GT( CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 0)), 0 );
    ASSERT_LT( CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 0)), 1 );
    ASSERT_LT( CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 1)), 0 );
    ASSERT_GT( CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 1)), 0 );

    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3_dual) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, NestedMinorSearchTwoSeparation)
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
    ASSERT_CMR_CALL( CMRtwoSum(cmr, K_3_3, K_3_3_dual, CMRrowToElement(1), CMRcolumnToElement(1), 3, &matrix) );

    bool isTU;
    CMR_MATROID_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_TRUE( isTU );
    ASSERT_TRUE( CMRmatroiddecHasTranspose(dec) ); /* As we test for graphicness, the transpose is constructed. */
    ASSERT_EQ( CMRmatroiddecType(dec), CMR_MATROID_DEC_TYPE_TWO_SUM );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 2UL );
    ASSERT_GT( CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 0)), 0 );
    ASSERT_LE( CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 0)), 0 );
    ASSERT_LT( CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 1)), 0 );
    ASSERT_GT( CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 1)), 0 );
    
    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3_dual) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, NestedMinorSearchTwoSeparationViolator)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "8 8 "
      " 1 1  0  0 0  0 0 0 "
      " 1 0  0 -1 0  0 0 0 "
      " 0 1  1  1 0  0 0 0 "
      " 0 0  1  1 0  0 0 0 "
      " 1 1  1  0 1  1 0 0 "
      " 1 1  1  0 1  0 1 0 "
      " 1 1 -1  0 0  0 1 1 "
      " 0 0  0  0 0 -1 1 1 "
    ) );

    bool isTU;
    CMR_MATROID_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_FALSE( isTU );
    ASSERT_EQ( CMRmatroiddecType(dec), CMR_MATROID_DEC_TYPE_SUBMATRIX );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 1UL );
    ASSERT_EQ( CMRmatroiddecType(CMRmatroiddecChild(dec, 0)), CMR_MATROID_DEC_TYPE_DETERMINANT );

    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, NestedMinorPivotsOneRowOneColumn)
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
  CMR_MATROID_DEC* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
  
  ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, NestedMinorPivotsTwoRowsOneColumn)
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
  CMR_MATROID_DEC* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  params.regular.directGraphicness = false;
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
  
  ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, NestedMinorPivotsOneRowTwoColumns)
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
  CMR_MATROID_DEC* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  params.regular.directGraphicness = false;
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
  
  ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, NestedMinorPivotsTwoSeparation)
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
  CMR_MATROID_DEC* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  params.regular.directGraphicness = false;
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
  
  ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

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
  CMR_MATROID_DEC* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  params.regular.directGraphicness = false;
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  if (knowNetwork)
  {
    ASSERT_TRUE( isTU );
    ASSERT_GT( CMRmatroiddecGraphicness(dec), 0 );
    CMR_CHRMAT* networkMatrix = NULL;
    bool isForest;
    ASSERT_CMR_CALL( CMRcomputeNetworkMatrix(cmr, CMRmatroiddecGraph(dec), &networkMatrix, NULL,
      CMRmatroiddecGraphArcsReversed(dec), CMRmatroiddecGraphSizeForest(dec), CMRmatroiddecGraphForest(dec),
      CMRmatroiddecGraphSizeCoforest(dec), CMRmatroiddecGraphCoforest(dec), &isForest) );

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
    ASSERT_LT( CMRmatroiddecGraphicness(dec), 0 );
  }

//   ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );

  ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
}

TEST(TU, SeqGraphicWheel)
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

TEST(TU, SeqGraphicOneRowOneColumn)
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

TEST(TU, SeqGraphicTwoRowsOneColumn)
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

TEST(TU, SeqGraphicOneRowTwoColumns)
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

TEST(TU, SeqGraphicnOneColumn)
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

TEST(TU, SeqGraphicOneRow)
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

TEST(TU, R10)
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
    CMR_MATROID_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );
    ASSERT_GT( CMRmatroiddecRegularity(dec), 0 );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 0UL );
    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  printf("\n\n");

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
    CMR_MATROID_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );
    ASSERT_TRUE( dec );
    ASSERT_GT( CMRmatroiddecRegularity(dec), 0 );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 0UL );
    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
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
  CMR_MATROID_DEC* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  params.regular.directGraphicness = false;
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  if (knowRegular)
  {
    ASSERT_TRUE( isTU );
  }
  else
  {
    ASSERT_FALSE( isTU );
  }

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );

  ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
}

TEST(TU, EnumerateRanksZeroTwo)
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

TEST(TU, EnumerateRanksOneOne)
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

TEST(TU, EnumerateRanksTwoZero)
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

TEST(TU, ThreeSumWideWideR12)
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
    CMR_MATROID_DEC* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.regular.threeSumStrategy = CMR_MATROID_DEC_THREESUM_FLAG_DISTRIBUTED_RANKS
      | CMR_MATROID_DEC_THREESUM_FLAG_FIRST_WIDE | CMR_MATROID_DEC_THREESUM_FLAG_SECOND_WIDE;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_GT( CMRmatroiddecRegularity(dec), 0 );
    ASSERT_LT( CMRmatroiddecGraphicness(dec), 0 );
    ASSERT_EQ( CMRmatroiddecType(dec), CMR_MATROID_DEC_TYPE_PIVOTS );
    ASSERT_EQ( CMRmatroiddecNumPivots(dec), 1UL );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 1UL );

    CMR_MATROID_DEC* child = CMRmatroiddecChild(dec, 0);

    ASSERT_EQ( CMRmatroiddecType(child), CMR_MATROID_DEC_TYPE_THREE_SUM );
    ASSERT_EQ( CMRmatroiddecNumChildren(child), 2UL );

    CMR_MATROID_DEC* grandChild1 = CMRmatroiddecChild(child, 0);
    CMR_MATROID_DEC* grandChild2 = CMRmatroiddecChild(child, 1);

    ASSERT_LT( CMRmatroiddecGraphicness(grandChild1), 0 );
    ASSERT_GT( CMRmatroiddecCographicness(grandChild1), 0 );

    ASSERT_GT( CMRmatroiddecGraphicness(grandChild2), 0 );

    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, ForbiddenSubmatrixWideWide)
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
    CMR_MATROID_DEC* dec = NULL;
    CMR_SUBMAT* forbiddenSubmatrix = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.regular.threeSumStrategy = CMR_MATROID_DEC_THREESUM_FLAG_DISTRIBUTED_RANKS
      | CMR_MATROID_DEC_THREESUM_FLAG_FIRST_WIDE | CMR_MATROID_DEC_THREESUM_FLAG_SECOND_WIDE;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, &forbiddenSubmatrix, &params, NULL, DBL_MAX) );

    ASSERT_LT( CMRmatroiddecRegularity(dec), 0 );
    ASSERT_EQ( forbiddenSubmatrix->numRows, 8UL );
    ASSERT_EQ( forbiddenSubmatrix->numColumns, 8UL );
    // TODO: Compute determinant once implemented.
    ASSERT_CMR_CALL( CMRsubmatFree(cmr, &forbiddenSubmatrix) );
    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, PartitionAlgorithm)
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
    ASSERT_CMR_CALL( CMRtwoSum(cmr, K_3_3, K_3_3_dual, CMRrowToElement(1), CMRcolumnToElement(1), 3, &twoSum) );

    size_t rowPermutations[] = { 4, 6, 5, 7, 0, 1, 2, 3 };
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRchrmatPermute(cmr, twoSum, rowPermutations, NULL, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &twoSum) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_MATROID_DEC* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.algorithm = CMR_TU_ALGORITHM_PARTITION;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );
    ASSERT_EQ( dec, (CMR_MATROID_DEC*) NULL );

    ASSERT_TRUE( isTU );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3_dual) );
  }

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
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.algorithm = CMR_TU_ALGORITHM_PARTITION;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, NULL, NULL, &params, NULL, DBL_MAX) );
    ASSERT_FALSE( isTU );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
