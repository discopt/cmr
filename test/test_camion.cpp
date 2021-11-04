#include <gtest/gtest.h>

#include "common.h"
#include <cmr/camion.h>

TEST(Camion, Change)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "10 10 "
    "+1 -1  0  0  0  0  0  0  0  0 "
    "-1 +1  0  0  0  0  0  0  0  0 "
    "0   0 +1  0  0  0  0 -1  0  0 "
    "0   0  0  0 -1  0 +1  0  0  0 "
    "0   0  0  0  0  0 -1 -1 -1  0 "
    "0   0  0 -1  0 +1  0  0  0  0 "
    "0   0  0  0 +1 -1  0  0  0  0 "
    "0   0 -1 +1  0  0  0  0  0  0 "
    "0   0  0  0  0  0  0  0 +1 -1 "
    "0   0  0  0  0  0  0 -1  0 +1 "
  ) );
  CMR_CHRMAT* check = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &check, "10 10 "
    "+1 -1  0  0  0  0  0  0  0  0 "
    "-1 +1  0  0  0  0  0  0  0  0 "
    "0   0 +1  0  0  0  0 -1  0  0 "
    "0   0  0  0 -1  0 +1  0  0  0 "
    "0   0  0  0  0  0 -1 -1 -1  0 "
    "0   0  0 -1  0 +1  0  0  0  0 "
    "0   0  0  0 +1 -1  0  0  0  0 "
    "0   0 -1 -1  0  0  0  0  0  0 "
    "0   0  0  0  0  0  0  0 -1 -1 "
    "0   0  0  0  0  0  0 -1  0 +1 "
  ) );
  CMR_CHRMAT* checkViolator = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkViolator, "3 3 "
    "-1 -1 0 "
    "0 1 -1 "
    "-1 0 1 "
  ) );

  CMR_SUBMAT* submatrix = NULL;
  CMR_CHRMAT* violator = NULL;

  bool alreadySigned;
  ASSERT_CMR_CALL( CMRtestCamionSigned(cmr, matrix, &alreadySigned, &submatrix) );
  ASSERT_FALSE(alreadySigned);
  ASSERT_TRUE(submatrix != NULL);
  ASSERT_CMR_CALL( CMRchrmatFilterSubmat(cmr, matrix, submatrix, &violator) );
  ASSERT_TRUE(CMRchrmatCheckEqual(violator, checkViolator));
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &violator) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &submatrix) );

  ASSERT_CMR_CALL( CMRcomputeCamionSigned(cmr, matrix, &alreadySigned, NULL) );
  ASSERT_FALSE(alreadySigned);
  ASSERT_TRUE(CMRchrmatCheckEqual(matrix, check));

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkViolator) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


TEST(Camion, Issue10)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "4 4 "
    "1 1 0 0 "
    "1 0 1 1 "
    "0 1 1 0 "
    "0 0 1 1 "
  ) );
  CMR_CHRMAT* check = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &check, "4 4 "
    "1 1  0  0 "
    "1 0 -1  1 "
    "0 1  1  0 "
    "0 0  1 -1 "
  ) );


  bool alreadySigned;
  ASSERT_CMR_CALL( CMRcomputeCamionSigned(cmr, matrix, &alreadySigned, NULL) );
  ASSERT_FALSE(alreadySigned);
  ASSERT_TRUE(CMRchrmatCheckEqual(matrix, check));

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
