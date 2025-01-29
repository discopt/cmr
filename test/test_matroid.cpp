#include <gtest/gtest.h>

#include <stdio.h>

#include "common.h"
#include <cmr/matroid.h>

TEST(Matroid, BinaryPivot)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Single pivot on a 1-entry. */
  {
    CMR_CHRMAT* matrix = NULL;
    stringToCharMatrix(cmr, &matrix, "10 10 "
      "1 1 0 0 0 1 0 1 0 0 "
      "1 0 0 0 0 1 1 0 1 0 "
      "0 0 0 0 1 1 0 0 0 0 "
      "0 0 0 1 1 0 0 0 0 0 "
      "0 0 1 1 0 0 0 0 0 0 "
      "0 1 1 0 0 0 0 0 0 0 "
      "0 0 0 0 0 0 0 0 1 0 "
      "0 0 0 0 0 1 0 0 0 1 "
      "1 0 0 0 0 0 1 0 1 1 "
      "1 1 0 0 0 1 0 0 0 0 "
    );
    CMR_CHRMAT* check = NULL;
    stringToCharMatrix(cmr, &check, "10 10 "
      "1 1 0 0 0 1 0 1 0 0 "
      "1 1 0 0 0 0 1 1 1 0 "
      "0 0 0 0 1 1 0 0 0 0 "
      "0 0 0 1 1 0 0 0 0 0 "
      "0 0 1 1 0 0 0 0 0 0 "
      "0 1 1 0 0 0 0 0 0 0 "
      "0 0 0 0 0 0 0 0 1 0 "
      "0 0 0 0 0 1 0 0 0 1 "
      "1 1 0 0 0 1 1 1 1 1 "
      "1 0 0 0 0 0 0 1 0 0 "
    );

    CMR_CHRMAT* result = NULL;
    ASSERT_CMR_CALL( CMRchrmatBinaryPivot(cmr, matrix, 0, 0, &result) );
    ASSERT_TRUE( CMRchrmatCheckEqual(result, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &result) );
  }

  /* Severals pivots. */
  {
    CMR_CHRMAT* matrix = NULL;
    stringToCharMatrix(cmr, &matrix, "10 10 "
      "1 1 0 0 0 1 0 1 0 0 "
      "1 0 0 0 0 1 1 0 1 0 "
      "0 0 0 0 1 1 0 0 0 0 "
      "0 0 0 1 1 0 0 0 0 0 "
      "0 0 1 1 0 0 0 0 0 0 "
      "0 1 1 0 0 0 0 0 0 0 "
      "0 0 0 0 0 0 0 0 1 0 "
      "0 0 0 0 0 1 0 0 0 1 "
      "1 0 0 0 0 0 1 0 1 1 "
      "1 1 0 0 0 1 0 0 0 0 "
    );
    CMR_CHRMAT* check = NULL;
    stringToCharMatrix(cmr, &check, "10 10 "
      "1 0 1 1 1 1 0 1 0 0 "
      "1 1 1 1 1 1 1 0 1 0 "
      "0 1 1 1 1 1 0 0 0 0 "
      "0 1 1 1 1 0 0 0 0 0 "
      "0 1 1 1 0 0 0 0 0 0 "
      "0 1 1 0 0 0 0 0 0 0 "
      "0 0 0 0 0 0 0 0 1 0 "
      "0 1 1 1 1 1 0 0 0 1 "
      "1 0 0 0 0 0 1 0 1 1 "
      "1 0 1 1 1 1 0 0 0 0 "
    );

    CMR_CHRMAT* result = NULL;
    size_t pivotRows[4] = {5, 4, 3, 2};
    size_t pivotColumns[4] = {2, 3, 4, 5};
    ASSERT_CMR_CALL( CMRchrmatBinaryPivots(cmr, matrix, 4, pivotRows, pivotColumns, &result) );

    ASSERT_TRUE( CMRchrmatCheckEqual(result, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &result) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Matroid, TernaryPivot)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Single pivot on a 1-entry. */
  {
    CMR_CHRMAT* matrix = NULL;
    stringToCharMatrix(cmr, &matrix, "10 10 "
      " 1 1 0 0 0 -1 0 1 0 0 "
      "-1 0 0 0 0  1 1 0 1 0 "
      " 0 0 0 0 1  1 0 0 0 0 "
      " 0 0 0 1 1  0 0 0 0 0 "
      " 0 0 1 1 0  0 0 0 0 0 "
      " 0 1 1 0 0  0 0 0 0 0 "
      " 0 0 0 0 0  0 0 0 1 0 "
      " 0 0 0 0 0  1 0 0 0 1 "
      " 1 0 0 0 0  0 1 0 1 1 "
      " 1 1 0 0 0  1 0 0 0 0 "
    );
    CMR_CHRMAT* check = NULL;
    stringToCharMatrix(cmr, &check, "10 10 "
      "-1  1 0 0 0 -1 0  1 0 0 "
      "-1  1 0 0 0  0 1  1 1 0 "
      " 0  0 0 0 1  1 0  0 0 0 "
      " 0  0 0 1 1  0 0  0 0 0 "
      " 0  0 1 1 0  0 0  0 0 0 "
      " 0  1 1 0 0  0 0  0 0 0 "
      " 0  0 0 0 0  0 0  0 1 0 "
      " 0  0 0 0 0  1 0  0 0 1 "
      " 1 -1 0 0 0  1 1 -1 1 1 "
      " 1  0 0 0 0 -1 0 -1 0 0 "
    );

    CMR_CHRMAT* result = NULL;
    ASSERT_CMR_CALL( CMRchrmatTernaryPivot(cmr, matrix, 0, 0, &result) );


    if (!CMRchrmatCheckEqual(result, check))
    {
      printf("Input\n");
      ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );
      printf("After pivoting on r%zu,c%zu:\n", 0UL, 0UL);
      ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, result, stdout, '0', false) );
      printf("Expected result:\n");
      ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, check, stdout, '0', false) );
    }

    ASSERT_TRUE( CMRchrmatCheckEqual(result, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &result) );
  }

  /* Single pivot on a -1-entry. */
  {
    CMR_CHRMAT* matrix = NULL;
    stringToCharMatrix(cmr, &matrix, "10 10 "
      "-1 1 0 0 0 -1 0 1 0 0 "
      "-1 0 0 0 0  1 1 0 1 0 "
      " 0 0 0 0 1  1 0 0 0 0 "
      " 0 0 0 1 1  0 0 0 0 0 "
      " 0 0 1 1 0  0 0 0 0 0 "
      " 0 1 1 0 0  0 0 0 0 0 "
      " 0 0 0 0 0  0 0 0 1 0 "
      " 0 0 0 0 0  1 0 0 0 1 "
      " 1 0 0 0 0  0 1 0 1 1 "
      " 1 1 0 0 0  1 0 0 0 0 "
    );
    CMR_CHRMAT* check = NULL;
    stringToCharMatrix(cmr, &check, "10 10 "
      " 1 -1 0 0 0  1 0 -1 0 0 "
      " 1 -1 0 0 0 -1 1 -1 1 0 "
      " 0  0 0 0 1  1 0  0 0 0 "
      " 0  0 0 1 1  0 0  0 0 0 "
      " 0  0 1 1 0  0 0  0 0 0 "
      " 0  1 1 0 0  0 0  0 0 0 "
      " 0  0 0 0 0  0 0  0 1 0 "
      " 0  0 0 0 0  1 0  0 0 1 "
      "-1  1 0 0 0 -1 1  1 1 1 "
      "-1 -1 0 0 0  0 0  1 0 0 "
    );

    CMR_CHRMAT* result = NULL;
    ASSERT_CMR_CALL( CMRchrmatTernaryPivot(cmr, matrix, 0, 0, &result) );


    if (!CMRchrmatCheckEqual(result, check))
    {
      printf("Input\n");
      ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );
      printf("After pivoting on r%zu,c%zu:\n", 0UL, 0UL);
      ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, result, stdout, '0', false) );
      printf("Expected result:\n");
      ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, check, stdout, '0', false) );
    }

    ASSERT_TRUE( CMRchrmatCheckEqual(result, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &result) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
