  #include <gtest/gtest.h>

#include "common.h"
#include <cmr/ctu.h>

TEST(ComplementTotalUnimodularity, Examples)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;

  {
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 4 5 "
      "1 0 0 1 1 "
      "1 1 0 0 1 "
      "0 1 1 0 1 "
      "0 0 1 1 1 "
    ) );

    bool isCTU;
    ASSERT_CMR_CALL( CMRtestComplementTotalUnimodularity(cmr, matrix, &isCTU, NULL, NULL, NULL, DBL_MAX) );
    
    ASSERT_TRUE(isCTU);
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 5 5 "
      "0 1 0 1 0 "
      "1 1 1 0 1 "
      "0 1 0 0 1 "
      "1 0 0 0 1 "
      "0 1 1 1 1 "
    ) );

    bool isCTU;
    size_t complementRow;
    size_t complementColumn;
    ASSERT_CMR_CALL( CMRtestComplementTotalUnimodularity(cmr, matrix, &isCTU, &complementRow, &complementColumn, NULL,
      DBL_MAX) );
    ASSERT_FALSE(isCTU);
    ASSERT_EQ(complementRow, 0UL);
    ASSERT_EQ(complementColumn, 0UL);
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
