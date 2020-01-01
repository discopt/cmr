#include <gtest/gtest.h>

#include "common.h"
#include <tu/sign.h>

TEST(Sign, Change)
{
  TU* tu;
  TU_SPARSE_CHAR matrix = stringToSparseChar("10 10 "
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
  TU_SPARSE_CHAR check = stringToSparseChar("10 10 "
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
  TU_SPARSE_CHAR checkViolator = stringToSparseChar("3 3 "
    "-1 -1 0 "
    "0 1 -1 "
    "-1 0 1 "
  );

  TUinit(&tu);
  TU_SUBMATRIX* submatrix = NULL;
  TU_SPARSE_CHAR violator;

  ASSERT_FALSE(TUtestSignChar(tu, &matrix, &submatrix));
  ASSERT_TRUE(submatrix != NULL);
  TUfilterSubmatrixChar(&matrix, submatrix, &violator);
  ASSERT_TRUE(TUcheckSparseEqualChar(&violator, &checkViolator));
  TUclearSparseChar(&violator);
  TUfreeSubmatrix(&submatrix);

  ASSERT_FALSE(TUcorrectSignChar(tu, &matrix, NULL));
  ASSERT_TRUE(TUcheckSparseEqualChar(&matrix, &check));

  TUclearSparseChar(&checkViolator);
  TUclearSparseChar(&check);
  TUclearSparseChar(&matrix);

  TUfree(&tu);
}
