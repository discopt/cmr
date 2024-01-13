#include <gtest/gtest.h>

#include "common.h"

#include <cmr/regular.h>
#include <cmr/separation.h>
#include <cmr/graphic.h>

TEST(Regular, OneSum)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* K_3_3 = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3, "5 4 "
      " 1 1 0 0 "
      " 1 1 1 0 "
      " 1 0 0 1 "
      " 0 1 1 1 "
      " 0 0 1 1 "
    ) );

    CMR_CHRMAT* K_3_3_dual = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3_dual, "4 5 "
      " 1 1 1 0 0 "
      " 1 1 0 1 0 "
      " 0 1 0 1 1 "
      " 0 0 1 1 1 "
    ) );

    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRoneSum(cmr, K_3_3, K_3_3_dual, &matrix) );

    bool isRegular;
    CMR_MATROID_DEC* dec = NULL;
    CMR_REGULAR_PARAMS params;
    ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
    params.planarityCheck = true;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_TRUE( isRegular );
    ASSERT_FALSE( CMRmatroiddecHasTranspose(dec) ); /* Default settings should mean that the transpose is never computed. */
    ASSERT_EQ( CMRmatroiddecType(dec), CMR_MATROID_DEC_TYPE_ONE_SUM );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 2UL );
    ASSERT_LT( CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 0))
      * CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 1)), 0 );
    ASSERT_LT( CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 0))
      * CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 1)), 0 );
    
    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3_dual) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Regular, SeriesParallelTwoSeparation)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* K_3_3 = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3, "5 4 "
      " 1 1 0 0 "
      " 1 1 1 0 "
      " 1 0 0 1 "
      " 0 1 1 1 "
      " 0 0 1 1 "
    ) );

    CMR_CHRMAT* K_3_3_dual = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3_dual, "4 5 "
      " 1 1 1 0 0 "
      " 1 1 0 1 0 "
      " 0 1 0 1 1 "
      " 0 0 1 1 1 "
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

    bool isRegular;
    CMR_MATROID_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_TRUE( isRegular );
    ASSERT_TRUE( CMRmatroiddecHasTranspose(dec) ); /* As we test for graphicness, the transpose is constructed. */
    ASSERT_EQ( CMRmatroiddecType(dec), CMR_MATROID_DEC_TYPE_TWO_SUM );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 2UL );
    ASSERT_EQ( CMRmatroiddecType(CMRmatroiddecChild(dec, 0)), CMR_MATROID_DEC_TYPE_GRAPH );
    ASSERT_EQ( CMRmatroiddecType(CMRmatroiddecChild(dec, 1)), CMR_MATROID_DEC_TYPE_COGRAPH );
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

TEST(Regular, NestedMinorSearchTwoSeparation)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* K_3_3 = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3, "5 4 "
      " 1 1 0 0 "
      " 1 1 1 0 "
      " 1 0 0 1 "
      " 0 1 1 1 "
      " 0 0 1 1 "
    ) );

    CMR_CHRMAT* K_3_3_dual = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &K_3_3_dual, "4 5 "
      " 1 1 1 0 0 "
      " 1 1 0 1 0 "
      " 0 1 0 1 1 "
      " 0 0 1 1 1 "
    ) );

    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRtwoSum(cmr, K_3_3, K_3_3_dual, CMRrowToElement(1), CMRcolumnToElement(1), 3, &matrix) );

//     CMRchrmatPrintDense(cmr, K_3_3, stdout, '0', true);
//     CMRchrmatPrintDense(cmr, K_3_3_dual, stdout, '0', true);
//     CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isRegular;
    CMR_MATROID_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_TRUE( isRegular );
    ASSERT_FALSE( CMRmatroiddecHasMatrix(dec) ); /* Default settings should mean that the matrix is not copied. */
    ASSERT_TRUE( CMRmatroiddecHasTranspose(dec) ); /* As we test for graphicness, the transpose is constructed. */
    ASSERT_EQ( CMRmatroiddecType(dec), CMR_MATROID_DEC_TYPE_TWO_SUM );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 2UL );
    int graphic = (CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 0)) ? 2 : 0) + (CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 1)) ? 1 : 0);
    int cographic = (CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 0)) ? 1 : 0) + (CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 1)) ? 2 : 0);
    ASSERT_EQ( graphic, cographic );
    ASSERT_NE( graphic, 0);
    ASSERT_NE( graphic, 3);

    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3_dual) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Regular, NestedMinorPivotsOneRowOneColumn)
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

  bool isRegular;
  CMR_MATROID_DEC* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.directGraphicness = false;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
  
  ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Regular, NestedMinorPivotsTwoRowsOneColumn)
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

  bool isRegular;
  CMR_MATROID_DEC* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.directGraphicness = false;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
  
  ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Regular, NestedMinorPivotsOneRowTwoColumns)
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

  bool isRegular;
  CMR_MATROID_DEC* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.directGraphicness = false;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );
  
  ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Regular, NestedMinorPivotsTwoSeparation)
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

  bool isRegular;
  CMR_MATROID_DEC* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.directGraphicness = false;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

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
  bool knowGraphic    /**< Whether the matrix is graphic. */
)
{
  printf("Testing matrix for graphicness via nested minor sequence:\n");
  ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );

  bool isRegular;
  CMR_MATROID_DEC* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.directGraphicness = false;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

  if (knowGraphic)
  {
    ASSERT_TRUE( isRegular );
    ASSERT_TRUE( CMRmatroiddecGraphicness(dec) );
    CMR_CHRMAT* graphicMatrix = NULL;
    bool isForest;
    ASSERT_CMR_CALL( CMRcomputeGraphicMatrix(cmr, CMRmatroiddecGraph(dec), &graphicMatrix, NULL,
      CMRmatroiddecGraphSizeForest(dec), CMRmatroiddecGraphForest(dec), CMRmatroiddecGraphSizeCoforest(dec),
      CMRmatroiddecGraphCoforest(dec), &isForest) );
    if (!CMRchrmatCheckEqual(matrix, graphicMatrix))
    {
      printf("Computed graph with different representation matrix:\n");
      CMRchrmatPrintDense(cmr, graphicMatrix, stdout, '0', true);
    }
    ASSERT_TRUE( CMRchrmatCheckEqual(matrix, graphicMatrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &graphicMatrix) );
    ASSERT_TRUE(isForest);
  }
  else
  {
    ASSERT_FALSE( CMRmatroiddecGraphicness(dec) );
  }

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );

  ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
}

TEST(Regular, SequenceGraphicnessWheel)
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
      "1 0 1 0 "
      "1 1 0 0 "
      "0 1 1 1 "
      "0 0 1 1 "
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

TEST(Regular, SequenceGraphicnessOneRowOneColumn)
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

TEST(Regular, SequenceGraphicnessTwoRowsOneColumn)
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
      "0 0 0 0   1   1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Regular, SequenceGraphicnessOneRowTwoColumns)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 8 "
      "1 0 0 1   0 0   0 0 "
      "1 1 0 0   0 0   0 0 "
      "0 1 1 0   1 0   1 0 "
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

TEST(Regular, SequenceGraphicnessOneColumn)
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

TEST(Regular, SequenceGraphicnessOneRow)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 5 "
      "1 1 0   0   0 "
      "0 1 1   0   1 "
      "1 0 1   1   0 "
      "              "
      "1 0 0   1   0 "
      "              "
      "0 1 0   0   1 "
      "              "
      "0 1 1   0   0 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
      "1 1 0   0  1 " 
      "0 1 1   1  0 "
      "1 0 1   0  0 "
      "             "
      "0 1 0   1  1 "
      "             "
      "0 0 0   1  1 "
    ) );
    testSequenceGraphicness(cmr, matrix, false);
  }

  {
    /* Runs into addition of a single row, at least 1 articulation point, but none begin part of all fundamental cycles
     * induced by 1-edges. */
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "10 9 "
      "1 1 0 0 0 0 0 0 0 "
      "0 1 1 1 0 0 0 0 0 "
      "0 1 1 0 0 0 0 0 0 "
      "1 0 1 0 0 0 0 0 0 "
      "0 0 1 1 0 0 0 0 1 "
      "1 0 1 1 1 0 0 0 0 "
      "0 0 0 1 0 1 0 0 0 "
      "0 0 0 0 0 1 1 1 0 "
      "0 0 0 0 0 1 0 0 0 "
      "0 0 0 0 0 0 1 0 0 "
    ) );
    testSequenceGraphicness(cmr, matrix, false);
  }
  
  {
    /* Runs into addition of a single row with a unique 1 articulation point, but for which the auxiliary graph is not
     * bipartite. */
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "8 10 "
      "1 0 1 0 1 0 0 0 0 0 "
      "1 1 1 1 1 0 0 0 0 0 "
      "0 0 0 0 1 1 0 0 0 0 "
      "0 1 0 1 1 0 0 1 1 0 "
      "0 1 1 0 1 1 1 0 0 0 "
      "0 0 0 0 0 0 1 1 0 0 "
      "0 0 0 1 0 0 0 0 0 1 "
      "0 0 1 0 0 0 0 0 0 1 "
    ) );
    testSequenceGraphicness(cmr, matrix, false);
  }

  {
    /* Runs into addition of a single row with a unique 1 articulation point and bipartite auxiliary graph, and a
     * 1-edge adjacent to the split node. */
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
      "1 1 0 0 1 "
      "1 1 1 0 1 "
      "0 1 1 1 0 "
      "1 0 0 1 0 "
      "0 1 1 1 1 "
    ) );
    testSequenceGraphicness(cmr, matrix, true);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Regular, R10)
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

    bool isRegular;
    CMR_MATROID_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, NULL, NULL, DBL_MAX) );
    ASSERT_GT( CMRmatroiddecRegularity(dec), 0 );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 0UL );
    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  printf("\n\n");

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
      "1 1 0 0 1 "
      "1 1 1 0 0 "
      "0 1 1 1 0 "
      "0 0 1 1 1 "
      "1 0 0 1 1 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isRegular;
    CMR_MATROID_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, NULL, NULL, DBL_MAX) );
    ASSERT_GT( CMRmatroiddecRegularity(dec), 0 );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 0 );
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

  bool isRegular;
  CMR_MATROID_DEC* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.directGraphicness = false;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

  if (knowRegular)
  {
    ASSERT_TRUE( isRegular );
  }
  else
  {
    ASSERT_FALSE( isRegular );
  }

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true) );

  ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
}

TEST(Regular, EnumerateRanksZeroTwo)
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

TEST(Regular, EnumerateRanksOneOne)
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

TEST(Regular, EnumerateRanksTwoZero)
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

TEST(Regular, R12)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 6 "
      "1 0 1 1 0 0 "
      "0 1 1 1 0 0 "
      "1 0 1 0 1 1 "
      "0 1 0 1 1 1 "
      "1 0 1 0 1 0 "
      "0 1 0 1 0 1 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isRegular;
    CMR_MATROID_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, NULL, NULL, DBL_MAX) );
    ASSERT_GT( CMRmatroiddecRegularity(dec), 0 );
    ASSERT_EQ( CMRmatroiddecNumChildren(dec), 2UL );
    size_t graphicChildren = (CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 0)) ? 2 : 0)
      + (CMRmatroiddecGraphicness(CMRmatroiddecChild(dec, 1)) ? 1 : 0);
    size_t cographicChildren = (CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 0)) ? 2 : 0)
      + (CMRmatroiddecCographicness(CMRmatroiddecChild(dec, 1)) ? 1 : 0);
    ASSERT_EQ( CMRmatroiddecNumChildren(CMRmatroiddecChild(dec, 0)), 0UL );
    ASSERT_EQ( CMRmatroiddecNumChildren(CMRmatroiddecChild(dec, 1)), 0UL );
    ASSERT_EQ( graphicChildren + cographicChildren, 3UL );
    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
