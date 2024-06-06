#include <gtest/gtest.h>

#include "common.h"
#include <cmr/ctu.h>

TEST(CTU, ExamplesWideWide)
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
    CMR_CTU_PARAMS params;
    ASSERT_CMR_CALL( CMRctuParamsInit(&params) );
    params.tu.seymour.threeSumStrategy = CMR_SEYMOUR_THREESUM_FLAG_DISTRIBUTED_RANKS
      | CMR_SEYMOUR_THREESUM_FLAG_FIRST_WIDE | CMR_SEYMOUR_THREESUM_FLAG_SECOND_WIDE;

    ASSERT_CMR_CALL( CMRctuTest(cmr, matrix, &isCTU, NULL, NULL, &params, NULL, DBL_MAX) );
    
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
    CMR_CTU_PARAMS params;
    ASSERT_CMR_CALL( CMRctuParamsInit(&params) );
    params.tu.seymour.threeSumStrategy = CMR_SEYMOUR_THREESUM_FLAG_DISTRIBUTED_RANKS
      | CMR_SEYMOUR_THREESUM_FLAG_FIRST_WIDE | CMR_SEYMOUR_THREESUM_FLAG_SECOND_WIDE;

    ASSERT_CMR_CALL( CMRctuTest(cmr, matrix, &isCTU, &complementRow, &complementColumn, &params, NULL, DBL_MAX) );
    ASSERT_FALSE(isCTU);
    ASSERT_EQ(complementRow, 0UL);
    ASSERT_EQ(complementColumn, 0UL);
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
