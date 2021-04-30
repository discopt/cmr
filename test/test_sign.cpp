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

  TU_SUBMAT* submatrix = NULL;
  TU_CHRMAT* violator = NULL;

  bool alreadySigned;
  ASSERT_TU_CALL( TUtestSignChr(tu, matrix, &alreadySigned, &submatrix) );
  ASSERT_FALSE(alreadySigned);
  ASSERT_TRUE(submatrix != NULL);
  TUchrsubmatFilter(tu, matrix, submatrix, &violator);
  ASSERT_TRUE(TUchrmatCheckEqual(violator, checkViolator));
  TUchrmatFree(tu, &violator);
  TUsubmatFree(tu, &submatrix);

  ASSERT_TU_CALL( TUcorrectSignChr(tu, matrix, &alreadySigned, NULL) );
  ASSERT_FALSE(alreadySigned);
  ASSERT_TRUE(TUchrmatCheckEqual(matrix, check));

  TUchrmatFree(tu, &checkViolator);
  TUchrmatFree(tu, &check);
  TUchrmatFree(tu, &matrix);

  TUfreeEnvironment(&tu);
}
