#include <gtest/gtest.h>

#include "common.h"

#include <cmr/balanced.h>

TEST(Balanced, Submatrix)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
      " 1 0 0 1 "
      " 1 1 0 0 "
      " 0 1 1 0 "
      " 0 0 1 1 "
    ) );

    bool isBalanced;
    CMR_BALANCED_PARAMS params;
    ASSERT_CMR_CALL( CMRbalancedParamsInit(&params) );
    params.algorithm = CMR_BALANCED_ALGORITHM_SUBMATRIX;
    ASSERT_CMR_CALL( CMRbalancedTest(cmr, matrix, &isBalanced, NULL, &params, NULL, DBL_MAX) );

    ASSERT_TRUE( isBalanced );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
      " 1 0 0  1 "
      " 1 1 0  0 "
      " 0 1 1  0 "
      " 0 0 1 -1 "
    ) );

    bool isBalanced;
    CMR_BALANCED_PARAMS params;
    ASSERT_CMR_CALL( CMRbalancedParamsInit(&params) );
    params.algorithm = CMR_BALANCED_ALGORITHM_SUBMATRIX;
    ASSERT_CMR_CALL( CMRbalancedTest(cmr, matrix, &isBalanced, NULL, &params, NULL, DBL_MAX) );

    ASSERT_FALSE( isBalanced );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
      " 1 0 0 0 1 "
      " 1 1 0 0 0 "
      " 0 1 1 0 0 "
      " 0 0 1 1 0 "
      " 0 0 0 1 1 "
    ) );

    bool isBalanced;
    CMR_BALANCED_PARAMS params;
    ASSERT_CMR_CALL( CMRbalancedParamsInit(&params) );
    params.algorithm = CMR_BALANCED_ALGORITHM_SUBMATRIX;
    ASSERT_CMR_CALL( CMRbalancedTest(cmr, matrix, &isBalanced, NULL, &params, NULL, DBL_MAX) );

    ASSERT_FALSE( isBalanced );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 5 "
      " 1 0 0 0 0 1 "
      " 1 1 0 0 0 0 "
      " 0 1 1 0 0 0 "
      " 0 0 1 1 0 0 "
      " 0 0 0 1 1 0 "
      " 0 0 0 0 1 1 "
    ) );

    bool isBalanced;
    CMR_BALANCED_PARAMS params;
    ASSERT_CMR_CALL( CMRbalancedParamsInit(&params) );
    params.algorithm = CMR_BALANCED_ALGORITHM_SUBMATRIX;
    ASSERT_CMR_CALL( CMRbalancedTest(cmr, matrix, &isBalanced, NULL, &params, NULL, DBL_MAX) );

    ASSERT_TRUE( isBalanced );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "5 4 "
      " 1 1 0  0 "
      " 1 1 1  0 "
      " 1 0 0 -1 "
      " 0 1 1  1 "
      " 0 0 1  0 "
    ) );

    bool isBalanced;
    CMR_BALANCED_PARAMS params;
    ASSERT_CMR_CALL( CMRbalancedParamsInit(&params) );
    params.algorithm = CMR_BALANCED_ALGORITHM_SUBMATRIX;
    ASSERT_CMR_CALL( CMRbalancedTest(cmr, matrix, &isBalanced, NULL, &params, NULL, DBL_MAX) );

    ASSERT_TRUE( isBalanced );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

