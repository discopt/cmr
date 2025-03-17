#define MASSIVE_RANDOM /* Uncomment to test a large number of random matrices. */

#define MASSIVE_RANDOM_REPETITIONS 1000000
#define MASSIVE_RANDOM_WIDTH 12
#define MASSIVE_RANDOM_HEIGHT 12
#define MASSIVE_RANDOM_PROBABILITY 0.3
#define MASSIVE_RANDOM_DECOMPOSE_STRATEGY CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_YSUM | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_THREESUM

#include <gtest/gtest.h>

#include "common.h"

#include <cmr/tu.h>
#include <cmr/separation.h>
#include <cmr/graphic.h>
#include <cmr/linear_algebra.h>

TEST(TU, EulerianAlgorithm)
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
    size_t specials[4] = { 0, SIZE_MAX, SIZE_MAX, 0 };
    ASSERT_CMR_CALL( CMRtwosumCompose(cmr, K_3_3, K_3_3_dual, &specials[0], NULL, NULL, &specials[3], 3, &twoSum) );

    size_t rowPermutations[] = { 4, 6, 5, 7, 0, 1, 2, 3 };
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRchrmatPermute(cmr, twoSum, rowPermutations, NULL, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &twoSum) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.algorithm = CMR_TU_ALGORITHM_EULERIAN;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );
    ASSERT_EQ( dec, (CMR_SEYMOUR_NODE*) NULL );

    ASSERT_TRUE( isTU );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &K_3_3_dual) );
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "12 12 "
      "1 1 1 0 1 0 1 1 1 1 1 1 "
      "1 0 1 0 1 0 1 1 1 1 1 0 "
      "0 1 1 0 0 0 0 0 0 0 0 0 "
      "0 1 1 1 0 0 0 0 0 0 0 0 "
      "0 1 1 1 1 0 0 0 0 0 0 0 "
      "0 1 1 1 1 1 0 0 0 0 0 0 "
      "0 1 1 1 1 1 1 0 0 0 0 0 "
      "0 0 0 0 0 0 1 1 0 0 0 0 "
      "0 0 0 0 0 0 0 1 1 0 0 0 "
      "0 0 0 0 0 0 0 0 1 1 0 0 "
      "0 0 0 0 0 0 0 0 0 1 1 0 "
      "0 0 0 0 0 0 0 0 0 0 1 1 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.algorithm = CMR_TU_ALGORITHM_EULERIAN;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, NULL, NULL, &params, NULL, DBL_MAX) );
    ASSERT_FALSE( isTU );
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
    size_t specials[4] = { 1, SIZE_MAX, SIZE_MAX, 1 };
    ASSERT_CMR_CALL( CMRtwosumCompose(cmr, K_3_3, K_3_3_dual, &specials[0], &specials[1], &specials[2], &specials[3], 3,
      &twoSum) );

    size_t rowPermutations[] = { 4, 6, 5, 7, 0, 1, 2, 3 };
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRchrmatPermute(cmr, twoSum, rowPermutations, NULL, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &twoSum) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.algorithm = CMR_TU_ALGORITHM_PARTITION;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );
    ASSERT_EQ( dec, (CMR_SEYMOUR_NODE*) NULL );

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

#if defined(MASSIVE_RANDOM)

static
CMR_ERROR checkDecompositionTreePartition(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEYMOUR_NODE* node, /**< Root of decomposition tree. */
  bool *psuccess          /**< Pointer for storing success. */
)
{
  assert(cmr);
  assert(node);
  assert(psuccess);

  if (!CMRseymourIsTernary(node))
    return CMR_OKAY;

  if (CMRseymourRegularity(node) != 0)
  {
    CMR_TU_PARAMS params;
    CMR_CALL( CMRtuParamsInit(&params) );
    params.algorithm = CMR_TU_ALGORITHM_PARTITION;
    bool isTU;
    CMR_CALL( CMRtuTest(cmr, CMRseymourGetMatrix(node), &isTU, NULL, NULL, &params, NULL, DBL_MAX) );

    if (isTU != (CMRseymourRegularity(node) > 0))
    {
      printf("========== The following node has wrong regularity status! ==========\n");
      CMR_CALL( CMRseymourPrint(cmr, node, stdout, false, true, true, true, true, true) );
      *psuccess = false;
    }
  }

  for (size_t c = 0; c < CMRseymourNumChildren(node); ++c)
  {
    CMR_CALL( checkDecompositionTreePartition(cmr, CMRseymourChild(node, c), psuccess) );
  }

  return CMR_OKAY;
}

#endif /* MASSIVE_RANDOM */

TEST(TU, Onesum)
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
    CMR_CHRMAT* matrices[2] = { K_3_3, K_3_3_dual };
    ASSERT_CMR_CALL( CMRonesumCompose(cmr, 2, matrices, &matrix) );

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.planarityCheck = true;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
    ASSERT_TRUE( isTU );
    ASSERT_FALSE( CMRseymourHasTranspose(dec) ); /* Default settings should mean that the transpose is never computed. */
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    ASSERT_EQ( CMRseymourNumChildren(dec), 2UL );
    ASSERT_LT( CMRseymourGraphicness(CMRseymourChild(dec, 0))
      * CMRseymourGraphicness(CMRseymourChild(dec, 1)), 0);
    ASSERT_LT( CMRseymourCographicness(CMRseymourChild(dec, 0))
      * CMRseymourCographicness(CMRseymourChild(dec, 1)), 0);

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
    size_t specials[4] = { 1, SIZE_MAX, SIZE_MAX, 1 };
    ASSERT_CMR_CALL( CMRtwosumCompose(cmr, K_3_3, K_3_3_dual, &specials[0], &specials[1], &specials[2], &specials[3], 3,
      &twoSum) );

    size_t rowPermutations[] = { 4, 6, 5, 7, 0, 1, 2, 3 };
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRchrmatPermute(cmr, twoSum, rowPermutations, NULL, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &twoSum) );

    CMRchrmatPrintDense(cmr, K_3_3, stdout, '0', true);
    CMRchrmatPrintDense(cmr, K_3_3_dual, stdout, '0', true);
    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
    ASSERT_TRUE( isTU );
    ASSERT_TRUE( CMRseymourHasTranspose(dec) ); /* As we test for graphicness, the transpose is constructed. */
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_TWOSUM );
    ASSERT_EQ( CMRseymourNumChildren(dec), 2UL );
    ASSERT_GT( CMRseymourGraphicness(CMRseymourChild(dec, 0)), 0 );
    ASSERT_LT( CMRseymourCographicness(CMRseymourChild(dec, 0)), 1 );
    ASSERT_LT( CMRseymourGraphicness(CMRseymourChild(dec, 1)), 0 );
    ASSERT_GT( CMRseymourCographicness(CMRseymourChild(dec, 1)), 0 );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
    size_t specials[4] = { 1, SIZE_MAX, SIZE_MAX, 1 };
    ASSERT_CMR_CALL( CMRtwosumCompose(cmr, K_3_3, K_3_3_dual, &specials[0], &specials[1], &specials[2], &specials[3], 3,
      &matrix) );

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
    ASSERT_TRUE( isTU );
    ASSERT_TRUE( CMRseymourHasTranspose(dec) ); /* As we test for graphicness, the transpose is constructed. */
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_TWOSUM );
    ASSERT_EQ( CMRseymourNumChildren(dec), 2UL );
    ASSERT_GT( CMRseymourGraphicness(CMRseymourChild(dec, 0)), 0 );
    ASSERT_LE( CMRseymourCographicness(CMRseymourChild(dec, 0)), 0 );
    ASSERT_LT( CMRseymourGraphicness(CMRseymourChild(dec, 1)), 0 );
    ASSERT_GT( CMRseymourCographicness(CMRseymourChild(dec, 1)), 0 );
    
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
    CMR_SEYMOUR_NODE* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
    ASSERT_FALSE( isTU );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_IRREGULAR );
    ASSERT_EQ( CMRseymourNumMinors(dec), 1UL );
    ASSERT_EQ( CMRminorType(CMRseymourMinor(dec, 0)), CMR_MINOR_TYPE_DETERMINANT );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
  CMR_SEYMOUR_NODE* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
  
  ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
  CMR_SEYMOUR_NODE* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  params.seymour.directGraphicness = false;
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
  
  ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
  CMR_SEYMOUR_NODE* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  params.seymour.directGraphicness = false;
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
  
  ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
  CMR_SEYMOUR_NODE* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  params.seymour.directGraphicness = false;
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
  
  ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
  CMR_SEYMOUR_NODE* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  params.seymour.directGraphicness = false;
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  if (knowNetwork)
  {
    ASSERT_TRUE( isTU );
    ASSERT_GT( CMRseymourGraphicness(dec), 0 );
    CMR_CHRMAT* networkMatrix = NULL;
    bool isForest;
    ASSERT_CMR_CALL( CMRnetworkComputeMatrix(cmr, CMRseymourGraph(dec), &networkMatrix, NULL,
      CMRseymourGraphArcsReversed(dec), CMRseymourGraphSizeForest(dec), CMRseymourGraphForest(dec),
      CMRseymourGraphSizeCoforest(dec), CMRseymourGraphCoforest(dec), &isForest) );

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
    ASSERT_LT( CMRseymourGraphicness(dec), 0 );
  }

//   ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, 0, true, true, true) );

  ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
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
    CMR_SEYMOUR_NODE* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );
    ASSERT_GT( CMRseymourRegularity(dec), 0 );
    ASSERT_EQ( CMRseymourNumChildren(dec), 0UL );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
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
    CMR_SEYMOUR_NODE* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );
    ASSERT_TRUE( dec );
    ASSERT_GT( CMRseymourRegularity(dec), 0 );
    ASSERT_EQ( CMRseymourNumChildren(dec), 0UL );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
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
  CMR_SEYMOUR_NODE* dec = NULL;
  CMR_TU_PARAMS params;
  ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
  params.seymour.directGraphicness = false;
  ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

  if (knowRegular)
  {
    ASSERT_TRUE( isTU );
  }
  else
  {
    ASSERT_FALSE( isTU );
  }

  ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

  ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
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

TEST(TU, DeltasumR12)
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
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

    ASSERT_GT( CMRseymourRegularity(dec), 0 );
    ASSERT_LT( CMRseymourGraphicness(dec), 0 );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_PIVOTS );
    ASSERT_EQ( CMRseymourNumPivots(dec), 1UL );
    ASSERT_EQ( CMRseymourNumChildren(dec), 1UL );

    CMR_SEYMOUR_NODE* child = CMRseymourChild(dec, 0);

    ASSERT_EQ( CMRseymourType(child), CMR_SEYMOUR_NODE_TYPE_DELTASUM );
    ASSERT_EQ( CMRseymourNumChildren(child), 2UL );

    CMR_SEYMOUR_NODE* grandChild1 = CMRseymourChild(child, 0);
    CMR_SEYMOUR_NODE* grandChild2 = CMRseymourChild(child, 1);

    ASSERT_LT( CMRseymourGraphicness(grandChild1), 0 );
    ASSERT_GT( CMRseymourCographicness(grandChild1), 0 );

    ASSERT_GT( CMRseymourGraphicness(grandChild2), 0 );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, YsumR12)
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
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_YSUM
      | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_PIVOT;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

    ASSERT_GT( CMRseymourRegularity(dec), 0 );
    ASSERT_LT( CMRseymourGraphicness(dec), 0 );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_PIVOTS );
    ASSERT_EQ( CMRseymourNumPivots(dec), 1UL );
    ASSERT_EQ( CMRseymourNumChildren(dec), 1UL );

    CMR_SEYMOUR_NODE* child = CMRseymourChild(dec, 0);

    ASSERT_EQ( CMRseymourType(child), CMR_SEYMOUR_NODE_TYPE_YSUM );
    ASSERT_EQ( CMRseymourNumChildren(child), 2UL );

    CMR_SEYMOUR_NODE* grandChild1 = CMRseymourChild(child, 0);
    CMR_SEYMOUR_NODE* grandChild2 = CMRseymourChild(child, 1);

    ASSERT_GT( CMRseymourGraphicness(grandChild1), 0 );
    ASSERT_GT( CMRseymourGraphicness(grandChild2), 0 );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, ThreesumR12)
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
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_TRUEMPER;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

    ASSERT_GT( CMRseymourRegularity(dec), 0 );
    ASSERT_LT( CMRseymourGraphicness(dec), 0 );

    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_THREESUM );
    ASSERT_EQ( CMRseymourNumChildren(dec), 2UL );

    CMR_SEYMOUR_NODE* child1 = CMRseymourChild(dec, 0);
    CMR_SEYMOUR_NODE* child2 = CMRseymourChild(dec, 1);

    ASSERT_LT( CMRseymourGraphicness(child1), 0 );
    ASSERT_GT( CMRseymourCographicness(child1), 0 );

    ASSERT_GT( CMRseymourGraphicness(child2), 0 );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, DeltasumSigns)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "10 10 "
      "0  0  1  0  1  1   0  0  0  0 "
      "1  0  0  1  1  0   0  0  0  0 "
      "0  0  0  1  0  1   0  0  0  0 "
      "1  1  1  0  0  0  -1  0  0  0 "
      "                              "
      "0  0  0 -1 -1 -1   0  1  0  1 "
      "0  0  0 -1 -1 -1   0  0  1  1 "
      "0  0  0  0  0  0   1  0  0  1 "
      "0  0  0  0  0  0   0  1  1  0 "
      "0  0  0  0  0  0   1  1  0  0 "
      "0  0  0  0  0  0   1  0  1  0 "
    ) );
    bool isTU;
    CMR_SEYMOUR_NODE* root;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &root, NULL, &params, NULL, DBL_MAX) );

    // TODO: The decomposition of this matrix is not checked manually, yet. assert(false);

    if (root)
      ASSERT_CMR_CALL( CMRseymourRelease(cmr, &root) );


    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, YsumSigns)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "10 10 "
      "0  0  1  0  1  1   0  0  0  0 "
      "1  0  0  1  1  0   0  0  0  0 "
      "0  0  0  1  0  1   0  0  0  0 "
      "1  1  1  0  0  0  -1  0  0  0 "
      "                              "
      "0  0  0 -1 -1 -1   0  1  0  1 "
      "0  0  0 -1 -1 -1   0  0  1  1 "
      "0  0  0  0  0  0   1  0  0  1 "
      "0  0  0  0  0  0   0  1  1  0 "
      "0  0  0  0  0  0   1  1  0  0 "
      "0  0  0  0  0  0   1  0  1  0 "
    ) );
    bool isTU;
    CMR_SEYMOUR_NODE* root;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_YSUM
      | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_PIVOT;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &root, NULL, &params, NULL, DBL_MAX) );

    // TODO: The decomposition of this matrix is not checked manually, yet. assert(false);

    if (root)
      ASSERT_CMR_CALL( CMRseymourRelease(cmr, &root) );


    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}



TEST(TU, ThreesumForbiddenSubmatrix)
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
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_SUBMAT* forbiddenSubmatrix = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, &forbiddenSubmatrix, &params, NULL, DBL_MAX) );

    ASSERT_LT( CMRseymourRegularity(dec), 0 );
    ASSERT_EQ( forbiddenSubmatrix->numRows, 8UL );
    ASSERT_EQ( forbiddenSubmatrix->numColumns, 8UL );
    // TODO: Compute determinant once implemented.
    ASSERT_CMR_CALL( CMRsubmatFree(cmr, &forbiddenSubmatrix) );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, YsumForbiddenSubmatrix)
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
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_SUBMAT* forbiddenSubmatrix = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_YSUM
      | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_PIVOT;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, &forbiddenSubmatrix, &params, NULL, DBL_MAX) );

    ASSERT_LT( CMRseymourRegularity(dec), 0 );
    ASSERT_EQ( forbiddenSubmatrix->numRows, 8UL );
    ASSERT_EQ( forbiddenSubmatrix->numColumns, 8UL );
    // TODO: Compute determinant once implemented.
    ASSERT_CMR_CALL( CMRsubmatFree(cmr, &forbiddenSubmatrix) );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, DeltasumForbiddenSubmatrix)
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
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_SUBMAT* forbiddenSubmatrix = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_TRUEMPER;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, &forbiddenSubmatrix, &params, NULL, DBL_MAX) );

    ASSERT_LT( CMRseymourRegularity(dec), 0 );
    ASSERT_EQ( forbiddenSubmatrix->numRows, 8UL );
    ASSERT_EQ( forbiddenSubmatrix->numColumns, 8UL );
    // TODO: Compute determinant once implemented.
    ASSERT_CMR_CALL( CMRsubmatFree(cmr, &forbiddenSubmatrix) );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
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
    CMR_SEYMOUR_NODE* dec = NULL;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_IRREGULAR );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
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
    CMR_SEYMOUR_NODE* dec = NULL;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_IRREGULAR );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, ThreesumPivotHighRank)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    /* Here, pivoting in a distributed rank case yields +/- 2 entries. */
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "7 7 "
      "1 1 1 0 0 0 0 "
      "0 0 1 1 1 1 0 "
      "0 1 0 0 0 -1 1 "
      "0 1 0 0 0 0 1 "
      "0 1 0 -1 -1 0 0 "
      "1 1 0 -1 0 0 1 "
      "0 0 0 1 1 0 0 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_TRUEMPER;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_PIVOTS );
    ASSERT_EQ( CMRseymourNumChildren(dec), 1 );
    CMR_SEYMOUR_NODE* threeSumNode = CMRseymourChild(dec, 0);
    ASSERT_EQ( CMRseymourType(threeSumNode), CMR_SEYMOUR_NODE_TYPE_THREESUM );
    ASSERT_EQ( CMRseymourNumChildren(threeSumNode), 2 );
    CMR_SEYMOUR_NODE* node = CMRseymourChild(threeSumNode, 0);
    ASSERT_EQ( CMRseymourType(node), CMR_SEYMOUR_NODE_TYPE_IRREGULAR );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "10 10 "
      "0 0 0 0 0 0 0 0 1 0 "
      "0 0 1 0 0 0 0 0 0 -1 "
      "0 1 0 0 0 0 0 0 0 0 "
      "0 -1 1 1 0 0 0 0 0 0 "
      "0 1 0 0 0 1 1 1 1 1 "
      "1 0 1 1 0 0 0 1 1 0 "
      "1 0 0 0 0 0 0 1 0 0 "
      "0 1 0 0 0 0 0 0 0 0 "
      "0 -1 0 1 0 0 0 1 0 0 "
      "0 0 0 1 1 -1 0 1 0 0 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_TRUEMPER;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_SERIES_PARALLEL );
    ASSERT_EQ( CMRseymourNumChildren(dec), 1 );
    CMR_SEYMOUR_NODE* pivotNode = CMRseymourChild(dec, 0);
    ASSERT_EQ( CMRseymourType(pivotNode), CMR_SEYMOUR_NODE_TYPE_PIVOTS );
    ASSERT_EQ( CMRseymourNumPivots(pivotNode), 1 );
    ASSERT_EQ( CMRseymourNumChildren(pivotNode), 1 );
    CMR_SEYMOUR_NODE* threeSumNode = CMRseymourChild(pivotNode, 0);
    ASSERT_EQ( CMRseymourType(threeSumNode), CMR_SEYMOUR_NODE_TYPE_THREESUM );
    ASSERT_EQ( CMRseymourNumChildren(threeSumNode), 2 );
    CMR_SEYMOUR_NODE* node = CMRseymourChild(threeSumNode, 0);
    ASSERT_EQ( CMRseymourType(node), CMR_SEYMOUR_NODE_TYPE_IRREGULAR );
    node = CMRseymourChild(threeSumNode, 1);
    ASSERT_EQ( CMRseymourType(node), CMR_SEYMOUR_NODE_TYPE_IRREGULAR );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(TU, CompleteTree)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "8 8 "
      "1 1  0 1 0 0  0 0 "
      "0 1  1 1 0 0  0 0 "
      "1 0 -1 1 0 0  0 0 "
      "1 1  1 0 0 0  0 0 "
      "0 0  0 0 1 0  1 1 "
      "0 0  0 0 1 1  0 1 "
      "0 0  0 0 0 1 -1 1 "
      "0 0  0 0 1 1  1 0 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* root = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &root, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, root, stdout, true, true, true, true, true, true) );

    ASSERT_EQ( CMRseymourRegularity(root), -1 );
    ASSERT_EQ( CMRseymourType(root), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    CMR_SEYMOUR_NODE* child0 = CMRseymourChild(root, 0);
    CMR_SEYMOUR_NODE* child1 = CMRseymourChild(root, 1);
    ASSERT_EQ( CMRseymourType(child0), CMR_SEYMOUR_NODE_TYPE_IRREGULAR );
    ASSERT_EQ( CMRseymourType(child1), CMR_SEYMOUR_NODE_TYPE_IRREGULAR );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &root) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  /* 1-sum of two irregular ones, the top-left is detected in direct network test. */
  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "8 8 "
      "1 1  0 1 0 0  0 0 "
      "0 1  1 0 0 0  0 0 "
      "1 0  1 0 0 0  0 0 "
      "1 0  0 0 0 0  0 0 "
      "0 0  0 0 1 0  1 1 "
      "0 0  0 0 1 1  0 1 "
      "0 0  0 0 0 1 -1 1 "
      "0 0  0 0 1 1  1 0 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.stopWhenIrregular = true;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

    ASSERT_EQ( CMRseymourRegularity(dec), -1 );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    CMR_SEYMOUR_NODE* child0 = CMRseymourChild(dec, 0);
    CMR_SEYMOUR_NODE* child1 = CMRseymourChild(dec, 1);

    int numIrregular = 0;
    numIrregular += (CMRseymourRegularity(child0) < 0) ? 1 : 0;
    numIrregular += (CMRseymourRegularity(child1) < 0) ? 1 : 0;
    ASSERT_EQ( numIrregular, 1);

    int numUnknown = 0;
    numUnknown += CMRseymourType(child0) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    numUnknown += CMRseymourType(child1) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    ASSERT_EQ( numUnknown, 1);

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  /* 1-sum of two irregular ones, the top-left is detected in sequence network test. */
  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "8 8 "
      "1 1  0 0 0 0  0 0 "
      "0 1  1 0 0 0  0 0 "
      "0 0 -1 1 0 0  0 0 "
      "1 0  0 1 0 0  0 0 "
      "0 0  0 0 1 0  1 1 "
      "0 0  0 0 1 1  0 1 "
      "0 0  0 0 0 1 -1 1 "
      "0 0  0 0 1 1  1 0 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.directGraphicness = false;
    params.seymour.stopWhenIrregular = true;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

    ASSERT_EQ( CMRseymourRegularity(dec), -1 );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    CMR_SEYMOUR_NODE* child0 = CMRseymourChild(dec, 0);
    CMR_SEYMOUR_NODE* child1 = CMRseymourChild(dec, 1);

    int numIrregular = 0;
    numIrregular += (CMRseymourRegularity(child0) < 0) ? 1 : 0;
    numIrregular += (CMRseymourRegularity(child1) < 0) ? 1 : 0;
    ASSERT_EQ( numIrregular, 1);

    int numUnknown = 0;
    numUnknown += CMRseymourType(child0) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    numUnknown += CMRseymourType(child1) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    ASSERT_EQ( numUnknown, 1);

    ASSERT_EQ( CMRseymourGetUsed(dec), 1 );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  /* 1-sum of two irregular ones, the top-left is binary series-parallel but not TU. */
  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "8 8 "
      "1 0  1 1 0 0  0 0 "
      "1 1  0 0 0 0  0 0 "
      "0 1 -1 1 0 0  0 0 "
      "1 1  1 1 0 0  0 0 "
      "0 0  0 0 1 0  1 1 "
      "0 0  0 0 1 1  0 1 "
      "0 0  0 0 0 1 -1 1 "
      "0 0  0 0 1 1  1 0 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.stopWhenIrregular = true;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

    ASSERT_EQ( CMRseymourRegularity(dec), -1 );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    CMR_SEYMOUR_NODE* child0 = CMRseymourChild(dec, 0);
    CMR_SEYMOUR_NODE* child1 = CMRseymourChild(dec, 1);

    int numIrregular = 0;
    numIrregular += (CMRseymourRegularity(child0) < 0) ? 1 : 0;
    numIrregular += (CMRseymourRegularity(child1) < 0) ? 1 : 0;
    ASSERT_EQ( numIrregular, 1);

    int numUnknown = 0;
    numUnknown += CMRseymourType(child0) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    numUnknown += CMRseymourType(child1) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    ASSERT_EQ( numUnknown, 1);

    ASSERT_EQ( CMRseymourGetUsed(dec), 1 );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  /* 1-sum of two irregular ones, triggers something I forgot :-) */
  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "12 10 "
      " 1 1  0  0 0  0 0  0  0 0 "
      " 1 0  0 -1 0  0 0  0  0 0 "
      " 0 1  1  1 0  0 0  0  0 0 "
      " 0 0  1  1 0  0 0  0  0 0 "
      " 1 1  1  0 1  0 0  0  0 0 "
      " 1 1  1  0 1  1 0  0  0 0 "
      " 1 1 -1  0 0  1 0  0  0 0 "
      " 0 0  0  0 0  0 0 -1  1 1 "
      " 0 0  0  0 0  0 1  0  1 1 "
      " 0 0  0  0 0  0 1  1  0 1 "
      " 0 0  0  0 0  0 0  1 -1 1 "
      " 0 0  0  0 0  0 1  1  1 0 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.directGraphicness = false;
    params.seymour.stopWhenIrregular = true;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

    ASSERT_EQ( CMRseymourRegularity(dec), -1 );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    CMR_SEYMOUR_NODE* child0 = CMRseymourChild(dec, 0);
    CMR_SEYMOUR_NODE* child1 = CMRseymourChild(dec, 1);

    int numIrregular = 0;
    numIrregular += (CMRseymourRegularity(child0) < 0) ? 1 : 0;
    numIrregular += (CMRseymourRegularity(child1) < 0) ? 1 : 0;
    ASSERT_EQ( numIrregular, 1);

    int numUnknown = 0;
    numUnknown += CMRseymourType(child0) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    numUnknown += CMRseymourType(child1) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    ASSERT_EQ( numUnknown, 1);

    ASSERT_EQ( CMRseymourGetUsed(dec), 1 );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  /* 1-sum of two irregular ones, the top-left is 4-connected but irregular. */
  if (false)
  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "8 8 "
      "1 1  0 1 0 0  0 0 "
      "0 1  1 1 0 0  0 0 "
      "1 0  1 1 0 0  0 0 "
      "1 1  1 0 0 0  0 0 "
      "0 0  0 0 1 0  1 1 "
      "0 0  0 0 1 1  0 1 "
      "0 0  0 0 0 1 -1 1 "
      "0 0  0 0 1 1  1 0 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

    ASSERT_EQ( CMRseymourRegularity(dec), -1 );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    CMR_SEYMOUR_NODE* child0 = CMRseymourChild(dec, 0);
    CMR_SEYMOUR_NODE* child1 = CMRseymourChild(dec, 1);

    int numIrregular = 0;
    numIrregular += (CMRseymourRegularity(child0) < 0) ? 1 : 0;
    numIrregular += (CMRseymourRegularity(child1) < 0) ? 1 : 0;
    ASSERT_EQ( numIrregular, 1);

    int numUnknown = 0;
    numUnknown += CMRseymourType(child0) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    numUnknown += CMRseymourType(child1) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    ASSERT_EQ( numUnknown, 1);

    ASSERT_EQ( CMRseymourGetUsed(dec), 1 );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  /* 1-sum of two irregular ones, triggers something I forgot :-) */
  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "9 9 "
      "-1 0 0 1 1 0 0  0 0 "
      " 1 1 0 0 1 0 0  0 0 "
      " 0 1 1 0 1 0 0  0 0 "
      " 0 0 1 1 1 0 0  0 0 "
      " 1 1 1 1 1 0 0  0 0 "
      " 0 0 0 0 0 1 0  1 1 "
      " 0 0 0 0 0 1 1  0 1 "
      " 0 0 0 0 0 0 1 -1 1 "
      " 0 0 0 0 0 1 1  1 0 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.stopWhenIrregular = true;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

    ASSERT_EQ( CMRseymourRegularity(dec), -1 );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    CMR_SEYMOUR_NODE* child0 = CMRseymourChild(dec, 0);
    CMR_SEYMOUR_NODE* child1 = CMRseymourChild(dec, 1);

    int numIrregular = 0;
    numIrregular += (CMRseymourRegularity(child0) < 0) ? 1 : 0;
    numIrregular += (CMRseymourRegularity(child1) < 0) ? 1 : 0;
    ASSERT_EQ( numIrregular, 1);

    int numUnknown = 0;
    numUnknown += CMRseymourType(child0) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    numUnknown += CMRseymourType(child1) == CMR_SEYMOUR_NODE_TYPE_UNKNOWN ? 1 : 0;
    ASSERT_EQ( numUnknown, 1);

    ASSERT_EQ( CMRseymourGetUsed(dec), 1 );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  /* 1-sum of two irregular ones, the top-left is detected in direct network test, and the second will is completed
   * later. */
  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "8 8 "
      "1 1  0 1 0 0  0 0 "
      "0 1  1 0 0 0  0 0 "
      "1 0  1 0 0 0  0 0 "
      "1 0  0 0 0 0  0 0 "
      "0 0  0 0 1 0  1 1 "
      "0 0  0 0 1 1  0 1 "
      "0 0  0 0 0 1  1 1 "
      "0 0  0 0 1 1  1 0 "
    ) );

    printf("TU.Complete for 1-sum with completion:\n");
    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.stopWhenIrregular = true;

    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

    ASSERT_EQ( CMRseymourRegularity(dec), -1 );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    CMR_SEYMOUR_NODE* child0 = CMRseymourChild(dec, 0);
    CMR_SEYMOUR_NODE* child1 = CMRseymourChild(dec, 1);

    ASSERT_LT( CMRseymourRegularity(child0), 0 );
    ASSERT_EQ( CMRseymourRegularity(child1), 0 );

    params.seymour.stopWhenIrregular = false;
    ASSERT_CMR_CALL( CMRtuCompleteDecomposition(cmr, child1, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

    ASSERT_LT( CMRseymourRegularity(child1), 0 );
    ASSERT_EQ( CMRseymourType(child1), CMR_SEYMOUR_NODE_TYPE_IRREGULAR );

    ASSERT_EQ( CMRseymourGetUsed(dec), 1 );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
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
"7 7 "
"1 1 1 0 0 0 0 "
"0 0 1 1 1 1 0 "
"0 1 0 0 0 -1 1 "
"0 1 0 0 0 0 1 "
"0 1 0 -1 -1 0 0 "
"1 1 0 -1 0 0 1 "
"0 0 0 1 1 0 0 "
  ) );

  size_t repetitions = MASSIVE_RANDOM_REPETITIONS;
  size_t numRows = MASSIVE_RANDOM_HEIGHT;
  size_t numColumns = MASSIVE_RANDOM_WIDTH;
  double probability1 = MASSIVE_RANDOM_PROBABILITY;
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

    bool isTU;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_SUBMAT* violatorSubmatrix = NULL;
    CMR_TU_PARAMS params;
    ASSERT_CMR_CALL( CMRtuParamsInit(&params) );
    params.seymour.stopWhenIrregular = false;
    params.seymour.decomposeStrategy = MASSIVE_RANDOM_DECOMPOSE_STRATEGY;
    ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &isTU, &dec, NULL, &params, NULL, DBL_MAX) );

    bool matroidDecompositionCorrect = true;
    ASSERT_CMR_CALL( checkDecompositionTreePartition(cmr, dec, &matroidDecompositionCorrect) );
    if (!matroidDecompositionCorrect)
    {
      printf("Matroid decomposition tree:\n");
      ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, false, true, true) );
    }
    ASSERT_TRUE( matroidDecompositionCorrect );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

    printf("-> %stotally unimodular.\n", isTU ? "" : "NOT ");

    if (!isTU) /* TODO: Implement submatrix extraction. */
    {
//       ASSERT_TRUE( violatorSubmatrix );

      if (violatorSubmatrix)
      {
        CMR_CHRMAT* violatorMatrix = NULL;
        ASSERT_CMR_CALL( CMRchrmatSlice(cmr, matrix, violatorSubmatrix, &violatorMatrix) );

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
