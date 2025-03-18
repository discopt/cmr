// #define MASSIVE_RANDOM /* Uncomment to test a large number of random matrices. */

#define MASSIVE_RANDOM_REPETITIONS 1000000
#define MASSIVE_RANDOM_WIDTH 12
#define MASSIVE_RANDOM_HEIGHT 12
#define MASSIVE_RANDOM_PROBABILITY 0.3
#define MASSIVE_RANDOM_DECOMPOSE_STRATEGY CMR_SEYMOUR_DECOMPOSE_FLAG_DISTRIBUTED_YSUM | CMR_SEYMOUR_DECOMPOSE_FLAG_CONCENTRATED_THREESUM

#include <gtest/gtest.h>

#include "common.h"

#include <cmr/regular.h>
#include <cmr/separation.h>
#include <cmr/graphic.h>

TEST(Regular, Onesum)
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
    CMR_CHRMAT* matrices[2] = { K_3_3, K_3_3_dual };
    ASSERT_CMR_CALL( CMRonesumCompose(cmr, 2, matrices, &matrix) );

    bool isRegular;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_REGULAR_PARAMS params;
    ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
    params.seymour.planarityCheck = true;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
    ASSERT_TRUE( isRegular );
    ASSERT_FALSE( CMRseymourHasTranspose(dec) ); /* Default settings should mean that the transpose is never computed. */
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    ASSERT_EQ( CMRseymourNumChildren(dec), 2UL );
    ASSERT_LT( CMRseymourGraphicness(CMRseymourChild(dec, 0))
      * CMRseymourGraphicness(CMRseymourChild(dec, 1)), 0 );
    ASSERT_LT( CMRseymourCographicness(CMRseymourChild(dec, 0))
      * CMRseymourCographicness(CMRseymourChild(dec, 1)), 0 );
    
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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

    CMR_CHRMAT* twosum = NULL;
    size_t specials[4] = { 1, SIZE_MAX, SIZE_MAX, 1 };
    ASSERT_CMR_CALL( CMRtwosumCompose(cmr, K_3_3, K_3_3_dual, &specials[0], &specials[1], &specials[2], &specials[3], 3,
      &twosum) );

    size_t rowPermutations[] = { 4, 6, 5, 7, 0, 1, 2, 3 };
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRchrmatPermute(cmr, twosum, rowPermutations, NULL, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &twosum) );

    CMRchrmatPrintDense(cmr, K_3_3, stdout, '0', true);
    CMRchrmatPrintDense(cmr, K_3_3_dual, stdout, '0', true);
    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isRegular;
    CMR_SEYMOUR_NODE* dec = NULL;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
    ASSERT_TRUE( isRegular );
    ASSERT_TRUE( CMRseymourHasTranspose(dec) ); /* As we test for graphicness, the transpose is constructed. */
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_TWOSUM );
    ASSERT_EQ( CMRseymourNumChildren(dec), 2UL );
    ASSERT_EQ( CMRseymourType(CMRseymourChild(dec, 0)), CMR_SEYMOUR_NODE_TYPE_GRAPH );
    ASSERT_EQ( CMRseymourType(CMRseymourChild(dec, 1)), CMR_SEYMOUR_NODE_TYPE_COGRAPH );
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
    size_t specials[4] = { 1, SIZE_MAX, SIZE_MAX, 1 };
    ASSERT_CMR_CALL( CMRtwosumCompose(cmr, K_3_3, K_3_3_dual, &specials[0], &specials[1], &specials[2], &specials[3], 3,
      &matrix) );

//     CMRchrmatPrintDense(cmr, K_3_3, stdout, '0', true);
//     CMRchrmatPrintDense(cmr, K_3_3_dual, stdout, '0', true);
//     CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isRegular;
    CMR_SEYMOUR_NODE* dec = NULL;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, NULL, NULL, DBL_MAX) );

    ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
    ASSERT_TRUE( isRegular );
    ASSERT_TRUE( CMRseymourHasTranspose(dec) ); /* As we test for graphicness, the transpose is constructed. */
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_TWOSUM );
    ASSERT_EQ( CMRseymourNumChildren(dec), 2UL );
    ASSERT_GT( CMRseymourGraphicness(CMRseymourChild(dec, 0)), 0);
    ASSERT_LE( CMRseymourCographicness(CMRseymourChild(dec, 0)), 0);
    ASSERT_LT( CMRseymourGraphicness(CMRseymourChild(dec, 1)), 0);
    ASSERT_GT( CMRseymourCographicness(CMRseymourChild(dec, 1)), 0);

//     int graphic = () ? 2 : 0) + (CMRseymourGraphicness(CMRseymourChild(dec, 1)) ? 1 : 0);
//     int cographic = (CMRseymourCographicness(CMRseymourChild(dec, 0)) ? 1 : 0) + (CMRseymourCographicness(CMRseymourChild(dec, 1)) ? 2 : 0);
//     ASSERT_EQ( graphic, cographic );
//     ASSERT_NE( graphic, 0);
//     ASSERT_NE( graphic, 3);

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
    CMR_SEYMOUR_NODE* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.seymour.directGraphicness = false;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
  
  ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
    CMR_SEYMOUR_NODE* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.seymour.directGraphicness = false;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
  
  ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
    CMR_SEYMOUR_NODE* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.seymour.directGraphicness = false;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );
  
  ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

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
    CMR_SEYMOUR_NODE* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.seymour.directGraphicness = false;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

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
  bool knowGraphic    /**< Whether the matrix is graphic. */
)
{
  printf("Testing matrix for graphicness via nested minor sequence:\n");
  ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );

  bool isRegular;
    CMR_SEYMOUR_NODE* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.seymour.directGraphicness = false;
  params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

  ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, false, false, false, false, false) );

  if (knowGraphic)
  {
    ASSERT_TRUE( isRegular );
    ASSERT_TRUE( CMRseymourGraphicness(dec) );
    CMR_CHRMAT* graphicMatrix = NULL;
    bool isForest;
    ASSERT_CMR_CALL( CMRgraphicComputeMatrix(cmr, CMRseymourGraph(dec), &graphicMatrix, NULL,
      CMRseymourGraphSizeForest(dec), CMRseymourGraphForest(dec), CMRseymourGraphSizeCoforest(dec),
      CMRseymourGraphCoforest(dec), &isForest) );
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
    ASSERT_LT( CMRseymourGraphicness(dec), 0 );
  }

  ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

  ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
}

TEST(Regular, SequenceGraphicnessWheel)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 6 "
      "1 1 0 0 0 0 "
      "0 1 1 0 0 0 "
      "1 0 0 1 0 0 "
      "0 0 1 0 0 1 "
      "0 0 0 0 1 1 "
      "0 0 0 1 1 0 "
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
    printf("Matrix 1 of 5\n");
    fflush(stdout);

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
    printf("Matrix 2 of 5\n");
    fflush(stdout);

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
    printf("Matrix 3 of 5\n");
    fflush(stdout);

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
    printf("Matrix 4 of 5\n");
    fflush(stdout);

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
    printf("Matrix 5 of 5\n");
    fflush(stdout);

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
    CMR_SEYMOUR_NODE* dec = NULL;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, NULL, NULL, DBL_MAX) );
    ASSERT_GT( CMRseymourRegularity(dec), 0 );
    ASSERT_EQ( CMRseymourNumChildren(dec), 0UL );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
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
    CMR_SEYMOUR_NODE* dec = NULL;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, NULL, NULL, DBL_MAX) );
    ASSERT_GT( CMRseymourRegularity(dec), 0 );
    ASSERT_EQ( CMRseymourNumChildren(dec), 0 );
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
  CMR* cmr,                                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,                         /**< Matrix to test. */
  CMR_SEYMOUR_DECOMPOSE_FLAG threeSumStrategy, /**< Strategy for the 3-sum. */
  bool knowRegular                            /**< Whether the matrix is regular. */
)
{
  printf("Testing matrix for regularity with a 3-separation:\n");
  ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );

  bool isRegular;
    CMR_SEYMOUR_NODE* dec = NULL;
  CMR_REGULAR_PARAMS params;
  ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
  params.seymour.directGraphicness = false;
  params.seymour.decomposeStrategy = threeSumStrategy;
  ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

  if (knowRegular)
  {
    ASSERT_TRUE( isRegular );
  }
  else
  {
    ASSERT_FALSE( isRegular );
  }

  ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, true, true, true) );

  ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
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
    testEnumerate(cmr, matrix, CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR, false);
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
    testEnumerate(cmr, matrix, CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR, false);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Regular, EnumerateConcentratedRankForcePivot)
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
    testEnumerate(cmr, matrix, CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR, false);
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

    ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );

    bool isRegular;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_REGULAR_PARAMS params;
    ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
    params.seymour.decomposeStrategy = CMR_SEYMOUR_DECOMPOSE_FLAG_SEYMOUR;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );
    ASSERT_GT( CMRseymourRegularity(dec), 0 );
//     ASSERT_EQ( CMRseymourNumChildren(dec), 2UL );
//     size_t graphicChildren = (CMRseymourGraphicness(CMRseymourChild(dec, 0)) ? 2 : 0)
//       + (CMRseymourGraphicness(CMRseymourChild(dec, 1)) ? 1 : 0);
//     size_t cographicChildren = (CMRseymourCographicness(CMRseymourChild(dec, 0)) ? 2 : 0)
//       + (CMRseymourCographicness(CMRseymourChild(dec, 1)) ? 1 : 0);
//     ASSERT_EQ( CMRseymourNumChildren(CMRseymourChild(dec, 0)), 0UL );
//     ASSERT_EQ( CMRseymourNumChildren(CMRseymourChild(dec, 1)), 0UL );
//     ASSERT_EQ( graphicChildren + cographicChildren, 3UL );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Regular, TreeFlagsStopNoncographic)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "9 9 "
      " 1 1 0 0 0 0 0 0 0 "
      " 1 1 1 0 0 0 0 0 0 "
      " 1 0 0 1 0 0 0 0 0 "
      " 0 1 1 1 0 0 0 0 0 "
      " 0 0 1 1 0 0 0 0 0 "
      " 0 0 0 0 1 1 1 0 0 "
      " 0 0 0 0 1 1 0 1 0 "
      " 0 0 0 0 0 1 0 1 1 "
      " 0 0 0 0 0 0 1 1 1 "
    ) );

    ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );

    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_REGULAR_PARAMS params;
    ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
    params.seymour.stopWhenNoncographic = true;
    params.seymour.planarityCheck = true;

    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, NULL, &dec, NULL, &params, NULL, DBL_MAX) );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    ASSERT_EQ( CMRseymourType(CMRseymourChild(dec, 0)), CMR_SEYMOUR_NODE_TYPE_GRAPH );
    ASSERT_EQ( CMRseymourType(CMRseymourChild(dec, 1)), CMR_SEYMOUR_NODE_TYPE_UNKNOWN );
    ASSERT_LT( CMRseymourCographicness(dec), 0 );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


TEST(Regular, TreeFlagsStopNongraphic)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "9 9 "
      " 1 1 1 0 0 0 0 0 0 "
      " 1 1 0 1 0 0 0 0 0 "
      " 0 1 0 1 1 0 0 0 0 "
      " 0 0 1 1 1 0 0 0 0 "
      " 0 0 0 0 0 1 1 0 0 "
      " 0 0 0 0 0 1 1 1 0 "
      " 0 0 0 0 0 1 0 0 1 "
      " 0 0 0 0 0 0 1 1 1 "
      " 0 0 0 0 0 0 0 1 1 "
    ) );

    ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );

    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_REGULAR_PARAMS params;
    ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
    params.seymour.stopWhenNongraphic = true;
    params.seymour.planarityCheck = true;

    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, NULL, &dec, NULL, &params, NULL, DBL_MAX) );
    ASSERT_EQ( CMRseymourType(dec), CMR_SEYMOUR_NODE_TYPE_ONESUM );
    ASSERT_LT( CMRseymourGraphicness(dec), 0 );

    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


#if defined(MASSIVE_RANDOM)

TEST(Regular, Random)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix,
"7 7 "
"1 1 1 0 0 0 0 "
"0 0 1 1 1 1 0 "
"0 1 0 0 0 1 1 "
"0 1 0 0 0 0 1 "
"0 1 0 1 1 0 0 "
"1 1 0 1 0 0 1 "
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

    printf("Regular.Random matrix:\n");
    ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );

    bool isRegular;
    CMR_SEYMOUR_NODE* dec = NULL;
    CMR_REGULAR_PARAMS params;
    ASSERT_CMR_CALL( CMRregularParamsInit(&params) );
    params.seymour.stopWhenIrregular = false;
    params.seymour.decomposeStrategy = MASSIVE_RANDOM_DECOMPOSE_STRATEGY;
    ASSERT_CMR_CALL( CMRregularTest(cmr, matrix, &isRegular, &dec, NULL, &params, NULL, DBL_MAX) );

    bool matroidDecompositionCorrect = true;
    // ASSERT_CMR_CALL( checkDecompositionTreePartition(cmr, dec, &matroidDecompositionCorrect) );
    if (!matroidDecompositionCorrect)
    {
      printf("Matroid decomposition tree:\n");
      ASSERT_CMR_CALL( CMRseymourPrint(cmr, dec, stdout, true, true, true, false, true, true) );
    }
    ASSERT_TRUE( matroidDecompositionCorrect );
    ASSERT_CMR_CALL( CMRseymourRelease(cmr, &dec) );

    printf("-> %sregular.\n", isRegular ? "" : "NOT ");

    if (numRows <= maxSizePartition || numColumns <= maxSizePartition)
    {
      // params.algorithm = CMR_TU_ALGORITHM_PARTITION;
      // bool partitionIsTU;
      // ASSERT_CMR_CALL( CMRtuTest(cmr, matrix, &partitionIsTU, NULL, NULL, &params, NULL, DBL_MAX) );
      // ASSERT_EQ( isRegular, partitionIsTU );
    }

    /* Cleanup. */
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

#endif /* MASSIVE_RANDOM */
