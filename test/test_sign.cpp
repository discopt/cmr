#include <gtest/gtest.h>

#include "common.h"
#include <tu/sign.h>

TEST(Sign, Change)
{
  TU* tu;
  TUcreateEnvironment(&tu);
  TU_CHRMAT* matrix = NULL;
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
  TU_CHRMAT* check = NULL;
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
  TU_CHRMAT* checkViolator = NULL;
  stringToCharMatrix(tu, &checkViolator, "3 3 "
    "-1 -1 0 "
    "0 1 -1 "
    "-1 0 1 "
  );

  TU_SUBMATRIX* submatrix = NULL;
  TU_CHRMAT* violator = NULL;

  ASSERT_FALSE(TUtestSignChar(tu, matrix, &submatrix));
  ASSERT_TRUE(submatrix != NULL);
  TUfilterCharSubmatrix(tu, matrix, submatrix, &violator);
  ASSERT_TRUE(TUchrmatCheckEqual(violator, checkViolator));
  TUchrmatFree(tu, &violator);
  TUfreeSubmatrix(tu, &submatrix);

  ASSERT_FALSE(TUcorrectSignChar(tu, matrix, NULL));
  ASSERT_TRUE(TUchrmatCheckEqual(matrix, check));

  TUchrmatFree(tu, &checkViolator);
  TUchrmatFree(tu, &check);
  TUchrmatFree(tu, &matrix);

  TUfreeEnvironment(&tu);
}
