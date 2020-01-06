#include <gtest/gtest.h>

#include "common.h"
#include <tu/sign.h>

TEST(Sign, Change)
{
  TU* tu;
  TUcreateEnvironment(&tu);
  TU_CHAR_MATRIX* matrix = NULL;
  stringToCharMatrix(tu, &matrix, "10 10 "
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
  TU_CHAR_MATRIX* check = NULL;
  stringToCharMatrix(tu, &check, "10 10 "
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
  TU_CHAR_MATRIX* checkViolator = NULL;
  stringToCharMatrix(tu, &checkViolator, "3 3 "
    "-1 -1 0 "
    "0 1 -1 "
    "-1 0 1 "
  );

  TU_SUBMATRIX* submatrix = NULL;
  TU_CHAR_MATRIX* violator = NULL;

  ASSERT_FALSE(TUtestSignChar(tu, matrix, &submatrix));
  ASSERT_TRUE(submatrix != NULL);
  TUfilterCharSubmatrix(tu, matrix, submatrix, &violator);
  ASSERT_TRUE(TUcheckCharMatrixEqual(violator, checkViolator));
  TUfreeCharMatrix(tu, &violator);
  TUfreeSubmatrix(tu, &submatrix);

  ASSERT_FALSE(TUcorrectSignChar(tu, matrix, NULL));
  ASSERT_TRUE(TUcheckCharMatrixEqual(matrix, check));

  TUfreeCharMatrix(tu, &checkViolator);
  TUfreeCharMatrix(tu, &check);
  TUfreeCharMatrix(tu, &matrix);

  TUfreeEnvironment(&tu);
}
