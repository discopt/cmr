#include <gtest/gtest.h>

#include "common.h"
#include <tu/sign.h>

TEST(Sign, Change)
{
  TU* tu;
  TU_MATRIX_CHAR matrix = stringToMatrixChar("10 10 "
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
  );
  TU_MATRIX_CHAR check = stringToMatrixChar("10 10 "
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
  );
  TU_MATRIX_CHAR checkViolator = stringToMatrixChar("3 3 "
    "-1 -1 0 "
    "0 1 -1 "
    "-1 0 1 "
  );

  TUinit(&tu);
  TU_SUBMATRIX* submatrix = NULL;
  TU_MATRIX_CHAR violator;

  ASSERT_FALSE(TUtestSignChar(tu, &matrix, &submatrix));
  ASSERT_TRUE(submatrix != NULL);
  TUfilterSubmatrixChar(&matrix, submatrix, &violator);
  ASSERT_TRUE(TUcheckMatrixEqualChar(&violator, &checkViolator));
  TUclearMatrixChar(&violator);
  TUfreeSubmatrix(&submatrix);

  ASSERT_FALSE(TUcorrectSignChar(tu, &matrix, NULL));
  ASSERT_TRUE(TUcheckMatrixEqualChar(&matrix, &check));

  TUclearMatrixChar(&checkViolator);
  TUclearMatrixChar(&check);
  TUclearMatrixChar(&matrix);

  TUfree(&tu);
}
