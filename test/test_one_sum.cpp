#include <gtest/gtest.h>

#include "common.h"
#include "../src/tu/one_sum.h"

TEST(OneSum, DoubleToDouble)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);

  TU_DBLMAT* matrix = NULL;
  stringToDoubleMatrix(tu, &matrix, "10 10 "
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

  size_t numComponents;
  TU_ONESUM_COMPONENT* components = NULL;
  size_t rowsToComponents[10];
  size_t columnsToComponents[10];
  size_t rowsToComponentRows[10];
  size_t columnsToComponentColumns[10];

  TU_DBLMAT* check = NULL;
  TU_DBLMAT* checkTranspose = NULL;
  ASSERT_TU_CALL( decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(double), sizeof(double), &numComponents, &components,
    rowsToComponents, columnsToComponents, rowsToComponentRows, columnsToComponentColumns) );

  ASSERT_EQ(numComponents, 6);
  stringToDoubleMatrix(tu, &check, "3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  stringToDoubleMatrix(tu, &checkTranspose, "3 3 "
    "1 3 0 "
    "2 4 5 "
    "0 0 6 "
  );
  ASSERT_TRUE(TUdblmatCheckEqual(check, (TU_DBLMAT*) components[0].matrix));
  ASSERT_TRUE(TUdblmatCheckEqual(checkTranspose, (TU_DBLMAT*) components[0].transpose));
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
  TUdblmatFree(tu, &check);
  TUdblmatFree(tu, &checkTranspose);

  stringToDoubleMatrix(tu, &check, "2 2 "
    "1 0 "
    "3 2 "
  );
  stringToDoubleMatrix(tu, &checkTranspose, "2 2 "
    "1 3 "
    "0 2 "
  );
  ASSERT_TRUE(TUdblmatCheckEqual(check, (TU_DBLMAT*) components[1].matrix));
  ASSERT_TRUE(TUdblmatCheckEqual(checkTranspose, (TU_DBLMAT*) components[1].transpose));
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
  TUdblmatFree(tu, &check);
  TUdblmatFree(tu, &checkTranspose);

  stringToDoubleMatrix(tu, &check, "1 0 "
  );
  stringToDoubleMatrix(tu, &checkTranspose, "0 1 "
  );
  ASSERT_TRUE(TUdblmatCheckEqual(check, (TU_DBLMAT*) components[2].matrix));
  ASSERT_TRUE(TUdblmatCheckEqual(checkTranspose, (TU_DBLMAT*) components[2].transpose));
  ASSERT_EQ(components[2].rowsToOriginal[0], 3);
  ASSERT_EQ(rowsToComponents[3], 2);
  ASSERT_EQ(rowsToComponentRows[3], 0);
  TUdblmatFree(tu, &check);
  TUdblmatFree(tu, &checkTranspose);

  stringToDoubleMatrix(tu, &check, "1 1 "
    "1 "
  );
  stringToDoubleMatrix(tu, &checkTranspose, "1 1 "
    "1 "
  );
  ASSERT_TRUE(TUdblmatCheckEqual(check, (TU_DBLMAT*) components[3].matrix));
  ASSERT_TRUE(TUdblmatCheckEqual(checkTranspose, (TU_DBLMAT*) components[3].transpose));
  ASSERT_EQ(components[3].rowsToOriginal[0], 4);
  ASSERT_EQ(components[3].columnsToOriginal[0], 0);
  ASSERT_EQ(rowsToComponents[4], 3);
  ASSERT_EQ(columnsToComponents[0], 3);
  ASSERT_EQ(rowsToComponentRows[4], 0);
  ASSERT_EQ(columnsToComponentColumns[0], 0);
  TUdblmatFree(tu, &check);
  TUdblmatFree(tu, &checkTranspose);

  stringToDoubleMatrix(tu, &check, "3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  stringToDoubleMatrix(tu, &checkTranspose, "3 3 "
    "1 2 0 "
    "0 3 0 "
    "0 4 5 "
  );
  ASSERT_TRUE(TUdblmatCheckEqual(check, (TU_DBLMAT*) components[4].matrix));
  ASSERT_TRUE(TUdblmatCheckEqual(checkTranspose, (TU_DBLMAT*) components[4].transpose));
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
  TUdblmatFree(tu, &check);
  TUdblmatFree(tu, &checkTranspose);

  stringToDoubleMatrix(tu, &check, "0 1 "
  );
  stringToDoubleMatrix(tu, &checkTranspose, "1 0 "
  );
  ASSERT_TRUE(TUdblmatCheckEqual(check, (TU_DBLMAT*) components[5].matrix));
  ASSERT_TRUE(TUdblmatCheckEqual(checkTranspose, (TU_DBLMAT*) components[5].transpose));
  ASSERT_EQ(components[5].columnsToOriginal[0], 9);
  ASSERT_EQ(columnsToComponents[9], 5);
  ASSERT_EQ(columnsToComponentColumns[9], 0);
  TUdblmatFree(tu, &check);
  TUdblmatFree(tu, &checkTranspose);

  for (int c = 0; c < numComponents; ++c)
  {
    TUdblmatFree(tu, (TU_DBLMAT**) &components[c].matrix);
    TUdblmatFree(tu, (TU_DBLMAT**) &components[c].transpose);
    TUfreeBlockArray(tu, &components[c].rowsToOriginal);
    TUfreeBlockArray(tu, &components[c].columnsToOriginal);
  }
  TUfreeBlockArray(tu, &components);

  TUdblmatFree(tu, &matrix);
  TUfreeEnvironment(&tu);
}

TEST(OneSum, IntToInt)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);

  TU_INTMAT* matrix = NULL;
  stringToIntMatrix(tu, &matrix, "10 10 "
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

  size_t numComponents;
  TU_ONESUM_COMPONENT* components = NULL;
  size_t rowsToComponents[10];
  size_t columnsToComponents[10];
  size_t rowsToComponentRows[10];
  size_t columnsToComponentColumns[10];

  TU_INTMAT* check = NULL;
  TU_INTMAT* checkTranspose = NULL;
  ASSERT_TU_CALL( decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(int), sizeof(int), &numComponents, &components,
    rowsToComponents, columnsToComponents, rowsToComponentRows, columnsToComponentColumns) );

  ASSERT_EQ(numComponents, 6);
  stringToIntMatrix(tu, &check, "3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  stringToIntMatrix(tu, &checkTranspose, "3 3 "
    "1 3 0 "
    "2 4 5 "
    "0 0 6 "
  );
  ASSERT_TRUE(TUintmatCheckEqual(check, (TU_INTMAT*) components[0].matrix));
  ASSERT_TRUE(TUintmatCheckEqual(checkTranspose, (TU_INTMAT*) components[0].transpose));
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
  TUintmatFree(tu, &check);
  TUintmatFree(tu, &checkTranspose);

  stringToIntMatrix(tu, &check, "2 2 "
    "1 0 "
    "3 2 "
  );
  stringToIntMatrix(tu, &checkTranspose, "2 2 "
    "1 3 "
    "0 2 "
  );
  ASSERT_TRUE(TUintmatCheckEqual(check, (TU_INTMAT*) components[1].matrix));
  ASSERT_TRUE(TUintmatCheckEqual(checkTranspose, (TU_INTMAT*) components[1].transpose));
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
  TUintmatFree(tu, &check);
  TUintmatFree(tu, &checkTranspose);

  stringToIntMatrix(tu, &check, "1 0 "
  );
  stringToIntMatrix(tu, &checkTranspose, "0 1 "
  );
  ASSERT_TRUE(TUintmatCheckEqual(check, (TU_INTMAT*) components[2].matrix));
  ASSERT_TRUE(TUintmatCheckEqual(checkTranspose, (TU_INTMAT*) components[2].transpose));
  ASSERT_EQ(components[2].rowsToOriginal[0], 3);
  ASSERT_EQ(rowsToComponents[3], 2);
  ASSERT_EQ(rowsToComponentRows[3], 0);
  TUintmatFree(tu, &check);
  TUintmatFree(tu, &checkTranspose);

  stringToIntMatrix(tu, &check, "1 1 "
    "1 "
  );
  stringToIntMatrix(tu, &checkTranspose, "1 1 "
    "1 "
  );
  ASSERT_TRUE(TUintmatCheckEqual(check, (TU_INTMAT*) components[3].matrix));
  ASSERT_TRUE(TUintmatCheckEqual(checkTranspose, (TU_INTMAT*) components[3].transpose));
  ASSERT_EQ(components[3].rowsToOriginal[0], 4);
  ASSERT_EQ(components[3].columnsToOriginal[0], 0);
  ASSERT_EQ(rowsToComponents[4], 3);
  ASSERT_EQ(columnsToComponents[0], 3);
  ASSERT_EQ(rowsToComponentRows[4], 0);
  ASSERT_EQ(columnsToComponentColumns[0], 0);
  TUintmatFree(tu, &check);
  TUintmatFree(tu, &checkTranspose);

  stringToIntMatrix(tu, &check, "3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  stringToIntMatrix(tu, &checkTranspose, "3 3 "
    "1 2 0 "
    "0 3 0 "
    "0 4 5 "
  );
  ASSERT_TRUE(TUintmatCheckEqual(check, (TU_INTMAT*) components[4].matrix));
  ASSERT_TRUE(TUintmatCheckEqual(checkTranspose, (TU_INTMAT*) components[4].transpose));
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
  TUintmatFree(tu, &check);
  TUintmatFree(tu, &checkTranspose);

  stringToIntMatrix(tu, &check, "0 1 "
  );
  stringToIntMatrix(tu, &checkTranspose, "1 0 "
  );
  ASSERT_TRUE(TUintmatCheckEqual(check, (TU_INTMAT*) components[5].matrix));
  ASSERT_TRUE(TUintmatCheckEqual(checkTranspose, (TU_INTMAT*) components[5].transpose));
  ASSERT_EQ(components[5].columnsToOriginal[0], 9);
  ASSERT_EQ(columnsToComponents[9], 5);
  ASSERT_EQ(columnsToComponentColumns[9], 0);
  TUintmatFree(tu, &check);
  TUintmatFree(tu, &checkTranspose);

  for (int c = 0; c < numComponents; ++c)
  {
    TUintmatFree(tu, (TU_INTMAT**) &components[c].matrix);
    TUintmatFree(tu, (TU_INTMAT**) &components[c].transpose);
    TUfreeBlockArray(tu, &components[c].rowsToOriginal);
    TUfreeBlockArray(tu, &components[c].columnsToOriginal);
  }
  TUfreeBlockArray(tu, &components);

  TUintmatFree(tu, &matrix);
  TUfreeEnvironment(&tu);
}

TEST(OneSum, CharToChar)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);

  TU_CHRMAT* matrix = NULL;
  stringToCharMatrix(tu, &matrix, "10 10 "
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

  size_t numComponents;
  TU_ONESUM_COMPONENT* components = NULL;
  size_t rowsToComponents[10];
  size_t columnsToComponents[10];
  size_t rowsToComponentRows[10];
  size_t columnsToComponentColumns[10];

  TU_CHRMAT* check = NULL;
  TU_CHRMAT* checkTranspose = NULL;
  ASSERT_TU_CALL( decomposeOneSum(tu, (TU_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components,
    rowsToComponents, columnsToComponents, rowsToComponentRows, columnsToComponentColumns) );

  ASSERT_EQ(numComponents, 6);
  stringToCharMatrix(tu, &check, "3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  stringToCharMatrix(tu, &checkTranspose, "3 3 "
    "1 3 0 "
    "2 4 5 "
    "0 0 6 "
  );
  ASSERT_TRUE(TUchrmatCheckEqual(check, (TU_CHRMAT*) components[0].matrix));
  ASSERT_TRUE(TUchrmatCheckEqual(checkTranspose, (TU_CHRMAT*) components[0].transpose));
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
  TUchrmatFree(tu, &check);
  TUchrmatFree(tu, &checkTranspose);

  stringToCharMatrix(tu, &check, "2 2 "
    "1 0 "
    "3 2 "
  );
  stringToCharMatrix(tu, &checkTranspose, "2 2 "
    "1 3 "
    "0 2 "
  );
  ASSERT_TRUE(TUchrmatCheckEqual(check, (TU_CHRMAT*) components[1].matrix));
  ASSERT_TRUE(TUchrmatCheckEqual(checkTranspose, (TU_CHRMAT*) components[1].transpose));
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
  TUchrmatFree(tu, &check);
  TUchrmatFree(tu, &checkTranspose);

  stringToCharMatrix(tu, &check, "1 0 "
  );
  stringToCharMatrix(tu, &checkTranspose, "0 1 "
  );
  ASSERT_TRUE(TUchrmatCheckEqual(check, (TU_CHRMAT*) components[2].matrix));
  ASSERT_TRUE(TUchrmatCheckEqual(checkTranspose, (TU_CHRMAT*) components[2].transpose));
  ASSERT_EQ(components[2].rowsToOriginal[0], 3);
  ASSERT_EQ(rowsToComponents[3], 2);
  ASSERT_EQ(rowsToComponentRows[3], 0);
  TUchrmatFree(tu, &check);
  TUchrmatFree(tu, &checkTranspose);

  stringToCharMatrix(tu, &check, "1 1 "
    "1 "
  );
  stringToCharMatrix(tu, &checkTranspose, "1 1 "
    "1 "
  );
  ASSERT_TRUE(TUchrmatCheckEqual(check, (TU_CHRMAT*) components[3].matrix));
  ASSERT_TRUE(TUchrmatCheckEqual(checkTranspose, (TU_CHRMAT*) components[3].transpose));
  ASSERT_EQ(components[3].rowsToOriginal[0], 4);
  ASSERT_EQ(components[3].columnsToOriginal[0], 0);
  ASSERT_EQ(rowsToComponents[4], 3);
  ASSERT_EQ(columnsToComponents[0], 3);
  ASSERT_EQ(rowsToComponentRows[4], 0);
  ASSERT_EQ(columnsToComponentColumns[0], 0);
  TUchrmatFree(tu, &check);
  TUchrmatFree(tu, &checkTranspose);

  stringToCharMatrix(tu, &check, "3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  stringToCharMatrix(tu, &checkTranspose, "3 3 "
    "1 2 0 "
    "0 3 0 "
    "0 4 5 "
  );
  ASSERT_TRUE(TUchrmatCheckEqual(check, (TU_CHRMAT*) components[4].matrix));
  ASSERT_TRUE(TUchrmatCheckEqual(checkTranspose, (TU_CHRMAT*) components[4].transpose));
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
  TUchrmatFree(tu, &check);
  TUchrmatFree(tu, &checkTranspose);

  stringToCharMatrix(tu, &check, "0 1 "
  );
  stringToCharMatrix(tu, &checkTranspose, "1 0 "
  );
  ASSERT_TRUE(TUchrmatCheckEqual(check, (TU_CHRMAT*) components[5].matrix));
  ASSERT_TRUE(TUchrmatCheckEqual(checkTranspose, (TU_CHRMAT*) components[5].transpose));
  ASSERT_EQ(components[5].columnsToOriginal[0], 9);
  ASSERT_EQ(columnsToComponents[9], 5);
  ASSERT_EQ(columnsToComponentColumns[9], 0);
  TUchrmatFree(tu, &check);
  TUchrmatFree(tu, &checkTranspose);

  for (int c = 0; c < numComponents; ++c)
  {
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].matrix);
    TUchrmatFree(tu, (TU_CHRMAT**) &components[c].transpose);
    TUfreeBlockArray(tu, &components[c].rowsToOriginal);
    TUfreeBlockArray(tu, &components[c].columnsToOriginal);
  }
  TUfreeBlockArray(tu, &components);

  TUchrmatFree(tu, &matrix);
  TUfreeEnvironment(&tu);
}
