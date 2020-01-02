#include <gtest/gtest.h>

#include "common.h"
#include "../src/tu/one_sum.h"

TEST(OneSum, IntToInt)
{
  TU* tu;
  TU_MATRIX_INT matrix;

  TUinit(&tu);

  matrix = stringToMatrixInt("10 10 "
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

  int numComponents;
  TU_ONESUM_COMPONENT* components;
  int rowsToComponents[10];
  int columnsToComponents[10];
  int rowsToComponentRows[10];
  int columnsToComponentColumns[10];

  TU_MATRIX_INT check, checkTranspose;
  decomposeOneSum(tu, (TU_MATRIX*) &matrix, sizeof(int), sizeof(int), &numComponents, &components,
    rowsToComponents, columnsToComponents, rowsToComponentRows, columnsToComponentColumns);

  ASSERT_EQ(numComponents, 6);

  check = stringToMatrixInt("3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  checkTranspose = stringToMatrixInt("3 3 "
    "1 3 0 "
    "2 4 5 "
    "0 0 6 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualInt(&check, (TU_MATRIX_INT*) &components[0].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualInt(&checkTranspose, (TU_MATRIX_INT*) &components[0].transpose));
  ASSERT_EQ(components[0].rowsToOriginal[0], 0);
  ASSERT_EQ(components[0].rowsToOriginal[1], 5);
  ASSERT_EQ(components[0].rowsToOriginal[2], 6);
  ASSERT_EQ(components[0].columnsToOriginal[0], 1);
  ASSERT_EQ(components[0].columnsToOriginal[1], 3);
  ASSERT_EQ(components[0].columnsToOriginal[2], 8);
  ASSERT_EQ(rowsToComponents[0], 0);
  ASSERT_EQ(rowsToComponents[5], 0);
  ASSERT_EQ(rowsToComponents[6], 0);
  ASSERT_EQ(columnsToComponents[1], 0);
  ASSERT_EQ(columnsToComponents[3], 0);
  ASSERT_EQ(columnsToComponents[8], 0);
  ASSERT_EQ(rowsToComponentRows[0], 0);
  ASSERT_EQ(rowsToComponentRows[5], 1);
  ASSERT_EQ(rowsToComponentRows[6], 2);
  ASSERT_EQ(columnsToComponentColumns[1], 0);
  ASSERT_EQ(columnsToComponentColumns[3], 1);
  ASSERT_EQ(columnsToComponentColumns[8], 2);
  TUclearMatrixInt(&check);
  TUclearMatrixInt(&checkTranspose);

  check = stringToMatrixInt("2 2 "
    "1 0 "
    "3 2 "
  );
  checkTranspose = stringToMatrixInt("2 2 "
    "1 3 "
    "0 2 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualInt(&check, (TU_MATRIX_INT*) &components[1].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualInt(&checkTranspose, (TU_MATRIX_INT*) &components[1].transpose));
  ASSERT_EQ(components[1].rowsToOriginal[0], 1);
  ASSERT_EQ(components[1].rowsToOriginal[1], 2);
  ASSERT_EQ(components[1].columnsToOriginal[0], 7);
  ASSERT_EQ(components[1].columnsToOriginal[1], 2);
  ASSERT_EQ(rowsToComponents[1], 1);
  ASSERT_EQ(rowsToComponents[2], 1);
  ASSERT_EQ(columnsToComponents[7], 1);
  ASSERT_EQ(columnsToComponents[2], 1);
  ASSERT_EQ(rowsToComponentRows[1], 0);
  ASSERT_EQ(rowsToComponentRows[2], 1);
  ASSERT_EQ(columnsToComponentColumns[7], 0);
  ASSERT_EQ(columnsToComponentColumns[2], 1);
  TUclearMatrixInt(&check);
  TUclearMatrixInt(&checkTranspose);

  check = stringToMatrixInt("1 0 "
  );
  checkTranspose = stringToMatrixInt("0 1 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualInt(&check, (TU_MATRIX_INT*) &components[2].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualInt(&checkTranspose, (TU_MATRIX_INT*) &components[2].transpose));
  ASSERT_EQ(components[2].rowsToOriginal[0], 3);
  ASSERT_EQ(rowsToComponents[3], 2);
  ASSERT_EQ(rowsToComponentRows[3], 0);
  TUclearMatrixInt(&check);
  TUclearMatrixInt(&checkTranspose);

  check = stringToMatrixInt("1 1 "
    "1 "
  );
  checkTranspose = stringToMatrixInt("1 1 "
    "1 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualInt(&check, (TU_MATRIX_INT*) &components[3].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualInt(&checkTranspose, (TU_MATRIX_INT*) &components[3].transpose));
  ASSERT_EQ(components[3].rowsToOriginal[0], 4);
  ASSERT_EQ(components[3].columnsToOriginal[0], 0);
  ASSERT_EQ(rowsToComponents[4], 3);
  ASSERT_EQ(columnsToComponents[0], 3);
  ASSERT_EQ(rowsToComponentRows[4], 0);
  ASSERT_EQ(columnsToComponentColumns[0], 0);
  TUclearMatrixInt(&check);
  TUclearMatrixInt(&checkTranspose);

  check = stringToMatrixInt("3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  checkTranspose = stringToMatrixInt("3 3 "
    "1 2 0 "
    "0 3 0 "
    "0 4 5 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualInt(&check, (TU_MATRIX_INT*) &components[4].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualInt(&checkTranspose, (TU_MATRIX_INT*) &components[4].transpose));
  ASSERT_EQ(components[4].rowsToOriginal[0], 7);
  ASSERT_EQ(components[4].rowsToOriginal[1], 8);
  ASSERT_EQ(components[4].rowsToOriginal[2], 9);
  ASSERT_EQ(components[4].columnsToOriginal[0], 4);
  ASSERT_EQ(components[4].columnsToOriginal[1], 5);
  ASSERT_EQ(components[4].columnsToOriginal[2], 6);
  ASSERT_EQ(rowsToComponents[7], 4);
  ASSERT_EQ(rowsToComponents[8], 4);
  ASSERT_EQ(rowsToComponents[9], 4);
  ASSERT_EQ(columnsToComponents[4], 4);
  ASSERT_EQ(columnsToComponents[5], 4);
  ASSERT_EQ(columnsToComponents[6], 4);
  ASSERT_EQ(rowsToComponentRows[7], 0);
  ASSERT_EQ(rowsToComponentRows[8], 1);
  ASSERT_EQ(rowsToComponentRows[9], 2);
  ASSERT_EQ(columnsToComponentColumns[4], 0);
  ASSERT_EQ(columnsToComponentColumns[5], 1);
  ASSERT_EQ(columnsToComponentColumns[6], 2);
  TUclearMatrixInt(&check);
  TUclearMatrixInt(&checkTranspose);

  check = stringToMatrixInt("0 1 "
  );
  checkTranspose = stringToMatrixInt("1 0 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualInt(&check, (TU_MATRIX_INT*) &components[5].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualInt(&checkTranspose, (TU_MATRIX_INT*) &components[5].transpose));
  ASSERT_EQ(components[5].columnsToOriginal[0], 9);
  ASSERT_EQ(columnsToComponents[9], 5);
  ASSERT_EQ(columnsToComponentColumns[9], 0);
  TUclearMatrixInt(&check);
  TUclearMatrixInt(&checkTranspose);

  for (int c = 0; c < numComponents; ++c)
  {
    TUclearMatrixInt((TU_MATRIX_INT*) &components[c].matrix);
    TUclearMatrixInt((TU_MATRIX_INT*) &components[c].transpose);
    free(components[c].rowsToOriginal);
    free(components[c].columnsToOriginal);
  }
  free(components);

  TUclearMatrixInt(&matrix);
  TUfree(&tu);
}

TEST(OneSum, CharToChar)
{
  TU* tu;
  TU_MATRIX_CHAR matrix;

  TUinit(&tu);

  matrix = stringToMatrixChar("10 10 "
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

  int numComponents;
  TU_ONESUM_COMPONENT* components;
  int rowsToComponents[10];
  int columnsToComponents[10];
  int rowsToComponentRows[10];
  int columnsToComponentColumns[10];

  TU_MATRIX_CHAR check, checkTranspose;
  decomposeOneSum(tu, (TU_MATRIX*) &matrix, sizeof(char), sizeof(char), &numComponents, &components,
    rowsToComponents, columnsToComponents, rowsToComponentRows, columnsToComponentColumns);

  ASSERT_EQ(numComponents, 6);

  check = stringToMatrixChar("3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  checkTranspose = stringToMatrixChar("3 3 "
    "1 3 0 "
    "2 4 5 "
    "0 0 6 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualChar(&check, (TU_MATRIX_CHAR*) &components[0].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualChar(&checkTranspose, (TU_MATRIX_CHAR*) &components[0].transpose));
  ASSERT_EQ(components[0].rowsToOriginal[0], 0);
  ASSERT_EQ(components[0].rowsToOriginal[1], 5);
  ASSERT_EQ(components[0].rowsToOriginal[2], 6);
  ASSERT_EQ(components[0].columnsToOriginal[0], 1);
  ASSERT_EQ(components[0].columnsToOriginal[1], 3);
  ASSERT_EQ(components[0].columnsToOriginal[2], 8);
  ASSERT_EQ(rowsToComponents[0], 0);
  ASSERT_EQ(rowsToComponents[5], 0);
  ASSERT_EQ(rowsToComponents[6], 0);
  ASSERT_EQ(columnsToComponents[1], 0);
  ASSERT_EQ(columnsToComponents[3], 0);
  ASSERT_EQ(columnsToComponents[8], 0);
  ASSERT_EQ(rowsToComponentRows[0], 0);
  ASSERT_EQ(rowsToComponentRows[5], 1);
  ASSERT_EQ(rowsToComponentRows[6], 2);
  ASSERT_EQ(columnsToComponentColumns[1], 0);
  ASSERT_EQ(columnsToComponentColumns[3], 1);
  ASSERT_EQ(columnsToComponentColumns[8], 2);
  TUclearMatrixChar(&check);
  TUclearMatrixChar(&checkTranspose);

  check = stringToMatrixChar("2 2 "
    "1 0 "
    "3 2 "
  );
  checkTranspose = stringToMatrixChar("2 2 "
    "1 3 "
    "0 2 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualChar(&check, (TU_MATRIX_CHAR*) &components[1].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualChar(&checkTranspose, (TU_MATRIX_CHAR*) &components[1].transpose));
  ASSERT_EQ(components[1].rowsToOriginal[0], 1);
  ASSERT_EQ(components[1].rowsToOriginal[1], 2);
  ASSERT_EQ(components[1].columnsToOriginal[0], 7);
  ASSERT_EQ(components[1].columnsToOriginal[1], 2);
  ASSERT_EQ(rowsToComponents[1], 1);
  ASSERT_EQ(rowsToComponents[2], 1);
  ASSERT_EQ(columnsToComponents[7], 1);
  ASSERT_EQ(columnsToComponents[2], 1);
  ASSERT_EQ(rowsToComponentRows[1], 0);
  ASSERT_EQ(rowsToComponentRows[2], 1);
  ASSERT_EQ(columnsToComponentColumns[7], 0);
  ASSERT_EQ(columnsToComponentColumns[2], 1);
  TUclearMatrixChar(&check);
  TUclearMatrixChar(&checkTranspose);

  check = stringToMatrixChar("1 0 "
  );
  checkTranspose = stringToMatrixChar("0 1 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualChar(&check, (TU_MATRIX_CHAR*) &components[2].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualChar(&checkTranspose, (TU_MATRIX_CHAR*) &components[2].transpose));
  ASSERT_EQ(components[2].rowsToOriginal[0], 3);
  ASSERT_EQ(rowsToComponents[3], 2);
  ASSERT_EQ(rowsToComponentRows[3], 0);
  TUclearMatrixChar(&check);
  TUclearMatrixChar(&checkTranspose);

  check = stringToMatrixChar("1 1 "
    "1 "
  );
  checkTranspose = stringToMatrixChar("1 1 "
    "1 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualChar(&check, (TU_MATRIX_CHAR*) &components[3].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualChar(&checkTranspose, (TU_MATRIX_CHAR*) &components[3].transpose));
  ASSERT_EQ(components[3].rowsToOriginal[0], 4);
  ASSERT_EQ(components[3].columnsToOriginal[0], 0);
  ASSERT_EQ(rowsToComponents[4], 3);
  ASSERT_EQ(columnsToComponents[0], 3);
  ASSERT_EQ(rowsToComponentRows[4], 0);
  ASSERT_EQ(columnsToComponentColumns[0], 0);
  TUclearMatrixChar(&check);
  TUclearMatrixChar(&checkTranspose);

  check = stringToMatrixChar("3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  checkTranspose = stringToMatrixChar("3 3 "
    "1 2 0 "
    "0 3 0 "
    "0 4 5 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualChar(&check, (TU_MATRIX_CHAR*) &components[4].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualChar(&checkTranspose, (TU_MATRIX_CHAR*) &components[4].transpose));
  ASSERT_EQ(components[4].rowsToOriginal[0], 7);
  ASSERT_EQ(components[4].rowsToOriginal[1], 8);
  ASSERT_EQ(components[4].rowsToOriginal[2], 9);
  ASSERT_EQ(components[4].columnsToOriginal[0], 4);
  ASSERT_EQ(components[4].columnsToOriginal[1], 5);
  ASSERT_EQ(components[4].columnsToOriginal[2], 6);
  ASSERT_EQ(rowsToComponents[7], 4);
  ASSERT_EQ(rowsToComponents[8], 4);
  ASSERT_EQ(rowsToComponents[9], 4);
  ASSERT_EQ(columnsToComponents[4], 4);
  ASSERT_EQ(columnsToComponents[5], 4);
  ASSERT_EQ(columnsToComponents[6], 4);
  ASSERT_EQ(rowsToComponentRows[7], 0);
  ASSERT_EQ(rowsToComponentRows[8], 1);
  ASSERT_EQ(rowsToComponentRows[9], 2);
  ASSERT_EQ(columnsToComponentColumns[4], 0);
  ASSERT_EQ(columnsToComponentColumns[5], 1);
  ASSERT_EQ(columnsToComponentColumns[6], 2);
  TUclearMatrixChar(&check);
  TUclearMatrixChar(&checkTranspose);

  check = stringToMatrixChar("0 1 "
  );
  checkTranspose = stringToMatrixChar("1 0 "
  );
  ASSERT_TRUE(TUcheckMatrixEqualChar(&check, (TU_MATRIX_CHAR*) &components[5].matrix));
  ASSERT_TRUE(TUcheckMatrixEqualChar(&checkTranspose, (TU_MATRIX_CHAR*) &components[5].transpose));
  ASSERT_EQ(components[5].columnsToOriginal[0], 9);
  ASSERT_EQ(columnsToComponents[9], 5);
  ASSERT_EQ(columnsToComponentColumns[9], 0);
  TUclearMatrixChar(&check);
  TUclearMatrixChar(&checkTranspose);

  for (int c = 0; c < numComponents; ++c)
  {
    TUclearMatrixChar((TU_MATRIX_CHAR*) &components[c].matrix);
    TUclearMatrixChar((TU_MATRIX_CHAR*) &components[c].transpose);
    free(components[c].rowsToOriginal);
    free(components[c].columnsToOriginal);
  }
  free(components);

  TUclearMatrixChar(&matrix);
  TUfree(&tu);
}
