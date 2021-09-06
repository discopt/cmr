  #include <gtest/gtest.h>

#include "common.h"

#include <cmr/regular.h>
#include <cmr/separation.h>

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
    CMR_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtestBinaryRegular(cmr, matrix, &isRegular, &dec, NULL, NULL) );

    ASSERT_CMR_CALL( CMRdecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_TRUE( isRegular );
    ASSERT_FALSE( CMRdecHasMatrix(dec) ); /* Default settings should mean that the matrix is not copied. */
    ASSERT_FALSE( CMRdecHasTranspose(dec) ); /* Default settings should mean that the transpose is never computed. */
    ASSERT_EQ( CMRdecIsSum(dec, NULL, NULL), 1 );
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

TEST(Regular, TwoSumDuringSeriesParallel)
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
    ASSERT_CMR_CALL( CMRtwoSum(cmr, K_3_3, K_3_3_dual, CMRrowToElement(1), CMRcolumnToElement(1), &twoSum) );

    size_t rowPermutations[] = { 4, 6, 5, 7, 0, 1, 2, 3 };
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( CMRchrmatPermute(cmr, twoSum, rowPermutations, NULL, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &twoSum) );

    CMRchrmatPrintDense(cmr, K_3_3, stdout, '0', true);
    CMRchrmatPrintDense(cmr, K_3_3_dual, stdout, '0', true);
    CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

    bool isRegular;
    CMR_DEC* dec = NULL;
    ASSERT_CMR_CALL( CMRtestBinaryRegular(cmr, matrix, &isRegular, &dec, NULL, NULL) );

    ASSERT_CMR_CALL( CMRdecPrint(cmr, dec, stdout, 0, true, true, true) );
    ASSERT_TRUE( isRegular );
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
