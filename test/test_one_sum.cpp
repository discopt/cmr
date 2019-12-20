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

  int numComponents;
  TU_ONESUM_COMPONENT_INT* components;
  int rowsToComponents[10];
  int columnsToComponents[10];
  int rowsToComponentRows[10];
  int columnsToComponentColumns[10];

  TU_SPARSE_INT check, checkTranspose;
  decomposeOneSumIntToInt(tu, &matrix, &numComponents, &components,  rowsToComponents,
    columnsToComponents, rowsToComponentRows, columnsToComponentColumns);

  ASSERT_EQ(numComponents, 6);

  check = stringToSparseInt("3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  checkTranspose = stringToSparseInt("3 3 "
    "1 3 0 "
    "2 4 5 "
    "0 0 6 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &components[0].matrix));
  ASSERT_TRUE(TUcheckSparseEqualInt(&checkTranspose, &components[0].transpose));
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
  TUclearSparseInt(&check);
  TUclearSparseInt(&checkTranspose);

  check = stringToSparseInt("2 2 "
    "1 0 "
    "3 2 "
  );
  checkTranspose = stringToSparseInt("2 2 "
    "1 3 "
    "0 2 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &components[1].matrix));
  ASSERT_TRUE(TUcheckSparseEqualInt(&checkTranspose, &components[1].transpose));
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
  TUclearSparseInt(&check);
  TUclearSparseInt(&checkTranspose);

  check = stringToSparseInt("1 0 "
  );
  checkTranspose = stringToSparseInt("0 1 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &components[2].matrix));
  ASSERT_TRUE(TUcheckSparseEqualInt(&checkTranspose, &components[2].transpose));
  ASSERT_EQ(components[2].rowsToOriginal[0], 3);
  ASSERT_EQ(rowsToComponents[3], 2);
  ASSERT_EQ(rowsToComponentRows[3], 0);
  TUclearSparseInt(&check);
  TUclearSparseInt(&checkTranspose);

  check = stringToSparseInt("1 1 "
    "1 "
  );
  checkTranspose = stringToSparseInt("1 1 "
    "1 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &components[3].matrix));
  ASSERT_TRUE(TUcheckSparseEqualInt(&checkTranspose, &components[3].transpose));
  ASSERT_EQ(components[3].rowsToOriginal[0], 4);
  ASSERT_EQ(components[3].columnsToOriginal[0], 0);
  ASSERT_EQ(rowsToComponents[4], 3);
  ASSERT_EQ(columnsToComponents[0], 3);
  ASSERT_EQ(rowsToComponentRows[4], 0);
  ASSERT_EQ(columnsToComponentColumns[0], 0);
  TUclearSparseInt(&check);
  TUclearSparseInt(&checkTranspose);

  check = stringToSparseInt("3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  checkTranspose = stringToSparseInt("3 3 "
    "1 2 0 "
    "0 3 0 "
    "0 4 5 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &components[4].matrix));
  ASSERT_TRUE(TUcheckSparseEqualInt(&checkTranspose, &components[4].transpose));
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
  TUclearSparseInt(&check);
  TUclearSparseInt(&checkTranspose);

  check = stringToSparseInt("0 1 "
  );
  checkTranspose = stringToSparseInt("1 0 "
  );
  ASSERT_TRUE(TUcheckSparseEqualInt(&check, &components[5].matrix));
  ASSERT_TRUE(TUcheckSparseEqualInt(&checkTranspose, &components[5].transpose));
  ASSERT_EQ(components[5].columnsToOriginal[0], 9);
  ASSERT_EQ(columnsToComponents[9], 5);
  ASSERT_EQ(columnsToComponentColumns[9], 0);
  TUclearSparseInt(&check);
  TUclearSparseInt(&checkTranspose);

  for (int c = 0; c < numComponents; ++c)
  {
    TUclearSparseInt(&components[c].matrix);
    TUclearSparseInt(&components[c].transpose);
    free(components[c].rowsToOriginal);
    free(components[c].columnsToOriginal);
  }
  free(components);

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

  int numComponents;
  TU_ONESUM_COMPONENT_CHAR* components;
  int rowsToComponents[10];
  int columnsToComponents[10];
  int rowsToComponentRows[10];
  int columnsToComponentColumns[10];

  TU_SPARSE_CHAR check, checkTranspose;
  decomposeOneSumCharToChar(tu, &matrix, &numComponents, &components,  rowsToComponents,
    columnsToComponents, rowsToComponentRows, columnsToComponentColumns);

  ASSERT_EQ(numComponents, 6);

  check = stringToSparseChar("3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  checkTranspose = stringToSparseChar("3 3 "
    "1 3 0 "
    "2 4 5 "
    "0 0 6 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &components[0].matrix));
  ASSERT_TRUE(TUcheckSparseEqualChar(&checkTranspose, &components[0].transpose));
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
  TUclearSparseChar(&check);
  TUclearSparseChar(&checkTranspose);

  check = stringToSparseChar("2 2 "
    "1 0 "
    "3 2 "
  );
  checkTranspose = stringToSparseChar("2 2 "
    "1 3 "
    "0 2 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &components[1].matrix));
  ASSERT_TRUE(TUcheckSparseEqualChar(&checkTranspose, &components[1].transpose));
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
  TUclearSparseChar(&check);
  TUclearSparseChar(&checkTranspose);

  check = stringToSparseChar("1 0 "
  );
  checkTranspose = stringToSparseChar("0 1 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &components[2].matrix));
  ASSERT_TRUE(TUcheckSparseEqualChar(&checkTranspose, &components[2].transpose));
  ASSERT_EQ(components[2].rowsToOriginal[0], 3);
  ASSERT_EQ(rowsToComponents[3], 2);
  ASSERT_EQ(rowsToComponentRows[3], 0);
  TUclearSparseChar(&check);
  TUclearSparseChar(&checkTranspose);

  check = stringToSparseChar("1 1 "
    "1 "
  );
  checkTranspose = stringToSparseChar("1 1 "
    "1 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &components[3].matrix));
  ASSERT_TRUE(TUcheckSparseEqualChar(&checkTranspose, &components[3].transpose));
  ASSERT_EQ(components[3].rowsToOriginal[0], 4);
  ASSERT_EQ(components[3].columnsToOriginal[0], 0);
  ASSERT_EQ(rowsToComponents[4], 3);
  ASSERT_EQ(columnsToComponents[0], 3);
  ASSERT_EQ(rowsToComponentRows[4], 0);
  ASSERT_EQ(columnsToComponentColumns[0], 0);
  TUclearSparseChar(&check);
  TUclearSparseChar(&checkTranspose);

  check = stringToSparseChar("3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  checkTranspose = stringToSparseChar("3 3 "
    "1 2 0 "
    "0 3 0 "
    "0 4 5 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &components[4].matrix));
  ASSERT_TRUE(TUcheckSparseEqualChar(&checkTranspose, &components[4].transpose));
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
  TUclearSparseChar(&check);
  TUclearSparseChar(&checkTranspose);

  check = stringToSparseChar("0 1 "
  );
  checkTranspose = stringToSparseChar("1 0 "
  );
  ASSERT_TRUE(TUcheckSparseEqualChar(&check, &components[5].matrix));
  ASSERT_TRUE(TUcheckSparseEqualChar(&checkTranspose, &components[5].transpose));
  ASSERT_EQ(components[5].columnsToOriginal[0], 9);
  ASSERT_EQ(columnsToComponents[9], 5);
  ASSERT_EQ(columnsToComponentColumns[9], 0);
  TUclearSparseChar(&check);
  TUclearSparseChar(&checkTranspose);

  for (int c = 0; c < numComponents; ++c)
  {
    TUclearSparseChar(&components[c].matrix);
    TUclearSparseChar(&components[c].transpose);
    free(components[c].rowsToOriginal);
    free(components[c].columnsToOriginal);
  }
  free(components);

  TUclearSparseChar(&matrix);
  TUfree(&tu);
}
