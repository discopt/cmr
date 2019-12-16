#include <gtest/gtest.h>

#include "common.h"
#include "../src/tu/one_sum.h"

TEST(OneSum, IntToInt)
{
  TU* tu;
  TU_SPARSE_INT matrix;

  TUinit(&tu);

  matrix = stringToSparseInt("10 10 "
    "0 1 0 2 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 1 0 0 "
    "0 0 2 0 0 0 0 3 0 0 "
    "0 0 0 0 0 0 0 0 0 0 "
    "1 0 0 0 0 0 0 0 0 0 "
    "0 3 0 4 0 0 0 0 0 0 "
    "0 0 0 5 0 0 0 0 6 0 "
    "0 0 0 0 1 0 0 0 0 0 "
    "0 0 0 0 2 3 4 0 0 0 "
    "0 0 0 0 0 0 5 0 0 0 "
  );

  TUprintSparseAsDenseInt(stdout, &matrix, ' ', true);

  int numComponents;
  TU_SPARSE_INT* compMatrices;
  TU_SPARSE_INT* compTransposes;
  TU_SPARSE_INT check;
  int** rowMapping;
  int** columnMapping;
  decomposeOneSumIntToInt(tu, &matrix, &numComponents, &compMatrices, &compTransposes, &rowMapping,
    &columnMapping);

  ASSERT_EQ(numComponents, 6);

  check = stringToSparseInt("3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &compMatrices[0]));
  ASSERT_EQ(rowMapping[0][0], 0);
  ASSERT_EQ(rowMapping[0][1], 5);
  ASSERT_EQ(rowMapping[0][2], 6);
  ASSERT_EQ(columnMapping[0][0], 1);
  ASSERT_EQ(columnMapping[0][1], 3);
  ASSERT_EQ(columnMapping[0][2], 8);
  TUclearSparseInt(&check);

  check = stringToSparseInt("2 2 "
    "1 0 "
    "3 2 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &compMatrices[1]));
  ASSERT_EQ(rowMapping[1][0], 1);
  ASSERT_EQ(rowMapping[1][1], 2);
  ASSERT_EQ(columnMapping[1][0], 7);
  ASSERT_EQ(columnMapping[1][1], 2);
  TUclearSparseInt(&check);

  check = stringToSparseInt("1 0 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &compMatrices[2]));
  ASSERT_EQ(rowMapping[2][0], 3);
  TUclearSparseInt(&check);

  check = stringToSparseInt("1 1 "
    "1 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &compMatrices[3]));
  ASSERT_EQ(rowMapping[3][0], 4);
  ASSERT_EQ(columnMapping[3][0], 0);
  TUclearSparseInt(&check);

  check = stringToSparseInt("3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &compMatrices[4]));
  ASSERT_EQ(rowMapping[4][0], 7);
  ASSERT_EQ(rowMapping[4][1], 8);
  ASSERT_EQ(rowMapping[4][2], 9);
  ASSERT_EQ(columnMapping[4][0], 4);
  ASSERT_EQ(columnMapping[4][1], 5);
  ASSERT_EQ(columnMapping[4][2], 6);
  TUclearSparseInt(&check);

  check = stringToSparseInt("0 1 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &compMatrices[5]));
  ASSERT_EQ(columnMapping[5][0], 9);
  TUclearSparseInt(&check);

  for (int c = 0; c < numComponents; ++c)
  {
    TUclearSparseInt(&compMatrices[c]);
    TUclearSparseInt(&compTransposes[c]);
    free(rowMapping[c]);
    free(columnMapping[c]);
  }
  free(compMatrices);
  free(compTransposes);
  free(rowMapping);
  free(columnMapping);

  TUclearSparseInt(&matrix);
  TUfree(&tu);
}

TEST(OneSum, CharToChar)
{
  TU* tu;
  TU_SPARSE_CHAR matrix;

  TUinit(&tu);

  matrix = stringToSparseChar("10 10 "
    "0 1 0 2 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 1 0 0 "
    "0 0 2 0 0 0 0 3 0 0 "
    "0 0 0 0 0 0 0 0 0 0 "
    "1 0 0 0 0 0 0 0 0 0 "
    "0 3 0 4 0 0 0 0 0 0 "
    "0 0 0 5 0 0 0 0 6 0 "
    "0 0 0 0 1 0 0 0 0 0 "
    "0 0 0 0 2 3 4 0 0 0 "
    "0 0 0 0 0 0 5 0 0 0 "
  );

  TUprintSparseAsDenseChar(stdout, &matrix, ' ', true);

  int numComponents;
  TU_SPARSE_CHAR* compMatrices;
  TU_SPARSE_CHAR* compTransposes;
  TU_SPARSE_CHAR check;
  int** rowMapping;
  int** columnMapping;
  decomposeOneSumCharToChar(tu, &matrix, &numComponents, &compMatrices, &compTransposes,
    &rowMapping, &columnMapping);

  ASSERT_EQ(numComponents, 6);

  check = stringToSparseChar("3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &compMatrices[0]));
  ASSERT_EQ(rowMapping[0][0], 0);
  ASSERT_EQ(rowMapping[0][1], 5);
  ASSERT_EQ(rowMapping[0][2], 6);
  ASSERT_EQ(columnMapping[0][0], 1);
  ASSERT_EQ(columnMapping[0][1], 3);
  ASSERT_EQ(columnMapping[0][2], 8);
  TUclearSparseChar(&check);

  check = stringToSparseChar("2 2 "
    "1 0 "
    "3 2 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &compMatrices[1]));
  ASSERT_EQ(rowMapping[1][0], 1);
  ASSERT_EQ(rowMapping[1][1], 2);
  ASSERT_EQ(columnMapping[1][0], 7);
  ASSERT_EQ(columnMapping[1][1], 2);
  TUclearSparseChar(&check);

  check = stringToSparseChar("1 0 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &compMatrices[2]));
  ASSERT_EQ(rowMapping[2][0], 3);
  TUclearSparseChar(&check);

  check = stringToSparseChar("1 1 "
    "1 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &compMatrices[3]));
  ASSERT_EQ(rowMapping[3][0], 4);
  ASSERT_EQ(columnMapping[3][0], 0);
  TUclearSparseChar(&check);

  check = stringToSparseChar("3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &compMatrices[4]));
  ASSERT_EQ(rowMapping[4][0], 7);
  ASSERT_EQ(rowMapping[4][1], 8);
  ASSERT_EQ(rowMapping[4][2], 9);
  ASSERT_EQ(columnMapping[4][0], 4);
  ASSERT_EQ(columnMapping[4][1], 5);
  ASSERT_EQ(columnMapping[4][2], 6);
  TUclearSparseChar(&check);

  check = stringToSparseChar("0 1 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &compMatrices[5]));
  ASSERT_EQ(columnMapping[5][0], 9);
  TUclearSparseChar(&check);

  for (int c = 0; c < numComponents; ++c)
  {
    TUclearSparseChar(&compMatrices[c]);
    TUclearSparseChar(&compTransposes[c]);
    free(rowMapping[c]);
    free(columnMapping[c]);
  }
  free(compMatrices);
  free(compTransposes);
  free(rowMapping);
  free(columnMapping);

  TUclearSparseChar(&matrix);
  TUfree(&tu);
}
