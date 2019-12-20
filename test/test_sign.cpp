#include <gtest/gtest.h>

#include "common.h"
#include <tu/sign.h>

TEST(Sign, Change)
{
  TU* tu;
  TU_SPARSE_CHAR matrix = stringToSparseChar("10 10 "
    "+1 -1  0  0  0  0  0  0  0  0 "
    "-1 -1  0  0  0  0  0  0  0  0 "
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

  TUinit(&tu);

  ASSERT_FALSE(TUtestSignChar(tu, &matrix));

//   TUprintSparseAsDenseChar(stdout, &matrix, ' ', true);

  ASSERT_FALSE(TUcorrectSignChar(tu, &matrix));

//   TUprintSparseAsDenseChar(stdout, &matrix, ' ', true);

  ASSERT_TRUE(TUcheckSparseEqualChar(&matrix, &check));

  TUclearSparseChar(&check);
  TUclearSparseChar(&matrix);

  TUfree(&tu);
}
