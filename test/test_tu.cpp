// #define MASSIVE_RANDOM /* Uncomment to test a large number of random matrices. */

#include <gtest/gtest.h>

#include "common.h"

#include <cmr/tu.h>
#include <cmr/separation.h>
#include <cmr/graphic.h>
#include <cmr/linear_algebra.h>

static
CMR_ERROR checkDecompositionTreePartition(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec, /**< Root of decomposition tree. */
  bool *psuccess        /**< Pointer for storing success. */
)
{
  assert(cmr);
  assert(dec);
  assert(psuccess);

  if (!CMRmatroiddecIsTernary(dec))
    return CMR_OKAY;

  if (CMRmatroiddecRegularity(dec) != 0)
  {
    CMR_TU_PARAMS params;
    CMR_CALL( CMRtuParamsInit(&params) );
    params.algorithm = CMR_TU_ALGORITHM_PARTITION;
    bool isTU;
    CMR_CALL( CMRtuTest(cmr, CMRmatroiddecGetMatrix(dec), &isTU, NULL, NULL, &params, NULL, DBL_MAX) );

    if (isTU != (CMRmatroiddecRegularity(dec) > 0))
    {
      printf("========== The following node has wrong regularity status! ==========\n");
      CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 2, false, true, true, true, true, true) );
      *psuccess = false;
    }
  }

  for (size_t c = 0; c < CMRmatroiddecNumChildren(dec); ++c)
  {
    CMR_CALL( checkDecompositionTreePartition(cmr, CMRmatroiddecChild(dec, c), psuccess) );
  }

  return CMR_OKAY;
}


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

    ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true, true, true, true) );
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

    ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true, true, true, true) );
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

    ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true, true, true, true) );
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

    ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true, true, true, true) );
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

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true, true, true, true) );
  
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

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true, true, true, true) );
  
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

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true, true, true, true) );
  
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

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true, true, true, true) );
  
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
    ASSERT_CMR_CALL( CMRnetworkComputeMatrix(cmr, CMRmatroiddecGraph(dec), &networkMatrix, NULL,
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

  ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 0, true, true, true, true, true, true) );

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

TEST(TU, Fano)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "3 4 "
      "1 1 0 1 "
      "0 1 1 1 "
      "1 0 1 1 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_MATROID_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_EQ( CMRmatroiddecType(dec), CMR_MATROID_DEC_TYPE_IRREGULAR );

    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
  
TEST(TU, FanoDual)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 3 "
      "1 1 0 "
      "0 1 1 "
      "1 0 1 "
      "1 1 1 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_MATROID_DEC* dec = NULL;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_EQ( CMRmatroiddecType(dec), CMR_MATROID_DEC_TYPE_IRREGULAR );

    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
  
TEST(TU, SubmatrixAlgorithm)
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
    params.algorithm = CMR_TU_ALGORITHM_SUBMATRIX;
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
    params.algorithm = CMR_TU_ALGORITHM_SUBMATRIX;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, NULL, NULL, &params, NULL, DBL_MAX) );
    ASSERT_FALSE( isTU );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

#if defined(MASSIVE_RANDOM)

TEST(TU, Random)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix,
"8 8 "
"0 0 1 1 0 0 0 0 "
"0 0 1 1 1 1 0 0 "
"0 0 0 0 0 0 1 0 "
"0 0 1 1 0 0 1 0 "
"0 0 0 1 0 0 1 0 "
"0 1 1 1 0 0 0 0 "
"0 1 0 0 0 0 -1 0 "
"0 1 0 0 -1 -1 1 0 "
  ) );

  size_t repetitions = 10000;
  size_t numRows = 20;
  size_t numColumns = 20;
  double probability1 = 0.2;
  size_t maxSizePartition = 100;

  for (size_t r = 0; r < repetitions; ++r)
  {
    if (!matrix)
    {
      size_t estimatedNumNonzeros = 1.1 * numRows * numColumns * probability1 + 1024;
      ASSERT_CMR_CALL( CMRchrmatCreate(cmr, &matrix, numRows, numColumns, estimatedNumNonzeros) );
      size_t entry = 0;
      for (size_t row = 0; row < numRows; ++row)
      {
        matrix->rowSlice[row] = entry;
        for (size_t column = 0; column < numColumns; ++column)
        {
          bool isNonzero = (rand() * 1.0 / RAND_MAX) < probability1;
          if (isNonzero)
          {
            if (entry == matrix->numNonzeros)
            {
              ASSERT_CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryColumns, 2*matrix->numNonzeros) );
              ASSERT_CMR_CALL( CMRreallocBlockArray(cmr, &matrix->entryValues, 2*matrix->numNonzeros) );
              matrix->numNonzeros *= 2;
            }
            matrix->entryColumns[entry] = column;
            matrix->entryValues[entry] = 1;
            ++entry;
          }
        }
      }
      matrix->rowSlice[numRows] = entry;
      matrix->numNonzeros = entry;
    }

    ASSERT_CMR_CALL( CMRcamionComputeSigns(cmr, matrix, NULL, NULL, NULL, DBL_MAX) );

    printf("TU.Random matrix:\n");
    ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );
    fflush(stdout);

    bool isTU;
    CMR_MATROID_DEC* dec = NULL;
    CMR_SUBMAT* violatorSubmatrix = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.regular.threeSumStrategy = CMR_MATROID_DEC_THREESUM_FLAG_SEYMOUR;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    bool matroidDecompositionCorrect = true;
    ASSERT_CMR_CALL( checkDecompositionTreePartition(cmr, dec, &matroidDecompositionCorrect) );
    if (!matroidDecompositionCorrect)
    {
      printf("Matroid decomposition tree:\n");
      ASSERT_CMR_CALL( CMRmatroiddecPrint(cmr, dec, stdout, 2, true, true, true, false, true, true) );
    }
    ASSERT_TRUE( matroidDecompositionCorrect );
    ASSERT_CMR_CALL( CMRmatroiddecFree(cmr, &dec) );

    printf("-> %stotally unimodular.\n", isTU ? "" : "NOT ");

    if (!isTU) /* TODO: Implement submatrix extraction. */
    {
//       ASSERT_TRUE( violatorSubmatrix );

      if (violatorSubmatrix)
      {
        CMR_CHRMAT* violatorMatrix = NULL;
        ASSERT_CMR_CALL( CMRchrmatZoomSubmat(cmr, matrix, violatorSubmatrix, &violatorMatrix) );

        int64_t determinant;
        ASSERT_CMR_CALL( CMRchrmatDeterminant(cmr, violatorMatrix, &determinant) );
        ASSERT_TRUE( (determinant > 1) || (determinant < -1) );

        ASSERT_CMR_CALL( CMRchrmatFree(cmr, &violatorMatrix) );
        ASSERT_CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
      }
    }

    if (numRows <= maxSizePartition || numColumns <= maxSizePartition)
    {
      params.algorithm = CMR_TU_ALGORITHM_PARTITION;
      bool partitionIsTU;
      ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &partitionIsTU, NULL, NULL, &params, NULL, DBL_MAX) );
      ASSERT_EQ( isTU, partitionIsTU );
    }

    /* Cleanup. */
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

#endif /* MASSIVE_RANDOM */
