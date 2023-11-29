#include <gtest/gtest.h>

#include "common.h"

#include <cmr/equimodular.h>
#include <cmr/separation.h>
#include <cmr/graphic.h>

TEST(Equimodular, Examples)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_INTMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToIntMatrix(cmr, &matrix, "5 9 "
      " -2 -1  1  5 -3 -2  -7  5  6 "
      "  1  0  1  7  1  2  -6  8  8 "
      " -3 -1  1  4 -4 -3  -7  4  5 "
      "  6  2  0  4  8  8   2  6  4 "
      "  4  2 -2 -9  6  4  13 -9 -11 "
    ) );


    bool isEquimodular;
    int64_t k = 0;
    CMR_EQUIMODULAR_STATISTICS stats;
    CMRstatsEquimodularityInit(&stats);
    ASSERT_CMR_CALL( CMRtestEquimodularity(cmr, matrix, &isEquimodular, &k, NULL, &stats, DBL_MAX) );

    ASSERT_TRUE(isEquimodular);
    ASSERT_EQ(k, 1);
    ASSERT_EQ(stats.totalCount, 1);

    ASSERT_CMR_CALL( CMRintmatFree(cmr, &matrix) );
  }

  {
    CMR_INTMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToIntMatrix(cmr, &matrix, "5 9 "
      " -2 -1  1  5 -3 -2  -7  5  6 "
      "  1  0  1  7  1  2  -6  8  8 "
      " -3 -1  1  4 -4 -3  -7  4  5 "
      "  6  2  0  4  8  8   2  6  4 "
      "  4  2 -2 -9  6  4  13 -9 -12 "
    ) );

    bool isEquimodular;
    int64_t k = 0;
    CMR_EQUIMODULAR_STATISTICS stats;
    CMRstatsEquimodularityInit(&stats);
    ASSERT_CMR_CALL( CMRtestEquimodularity(cmr, matrix, &isEquimodular, &k, NULL, &stats, DBL_MAX) );

    ASSERT_FALSE(isEquimodular);
    ASSERT_EQ(k, 0);
    ASSERT_EQ(stats.totalCount, 1);

    ASSERT_CMR_CALL( CMRintmatFree(cmr, &matrix) );
  }

  {
    CMR_INTMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToIntMatrix(cmr, &matrix, "4 4 "
      " 1 1 1 1 "
      " 1 1 0 0 "
      " 1 0 1 0 "
      " 1 0 0 1 "
    ) );


    bool isEquimodular;
    int64_t k = 0;
    CMR_EQUIMODULAR_STATISTICS stats;
    CMRstatsEquimodularityInit(&stats);
    ASSERT_CMR_CALL( CMRtestEquimodularity(cmr, matrix, &isEquimodular, &k, NULL, &stats, DBL_MAX) );

    ASSERT_TRUE(isEquimodular);
    ASSERT_EQ(k, 2);
    ASSERT_EQ(stats.totalCount, 1);

    ASSERT_CMR_CALL( CMRintmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
