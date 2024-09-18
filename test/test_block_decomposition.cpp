#include <gtest/gtest.h>

#include "common.h"
#include "../src/cmr/block_decomposition.h"

TEST(BlockDecomposition, DoubleToDouble)
{
  CMR* cmr = NULL;
  CMRcreateEnvironment(&cmr);

  CMR_DBLMAT* matrix = NULL;
  stringToDoubleMatrix(cmr, &matrix, "10 10 "
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
  CMR_BLOCK* components = NULL;
  size_t rowsToComponents[10];
  size_t columnsToComponents[10];
  size_t rowsToComponentRows[10];
  size_t columnsToComponentColumns[10];

  CMR_DBLMAT* check = NULL;
  CMR_DBLMAT* checkTranspose = NULL;
  ASSERT_CMR_CALL( CMRdecomposeBlocks(cmr, (CMR_MATRIX*) matrix, sizeof(double), sizeof(double), &numComponents, &components,
    rowsToComponents, columnsToComponents, rowsToComponentRows, columnsToComponentColumns) );

  ASSERT_EQ(numComponents, 6);
  stringToDoubleMatrix(cmr, &check, "3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  stringToDoubleMatrix(cmr, &checkTranspose, "3 3 "
    "1 3 0 "
    "2 4 5 "
    "0 0 6 "
  );
  ASSERT_TRUE(CMRdblmatCheckEqual(check, (CMR_DBLMAT*) components[0].matrix));
  ASSERT_TRUE(CMRdblmatCheckEqual(checkTranspose, (CMR_DBLMAT*) components[0].transpose));
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
  CMRdblmatFree(cmr, &check);
  CMRdblmatFree(cmr, &checkTranspose);

  stringToDoubleMatrix(cmr, &check, "2 2 "
    "1 0 "
    "3 2 "
  );
  stringToDoubleMatrix(cmr, &checkTranspose, "2 2 "
    "1 3 "
    "0 2 "
  );
  ASSERT_TRUE(CMRdblmatCheckEqual(check, (CMR_DBLMAT*) components[1].matrix));
  ASSERT_TRUE(CMRdblmatCheckEqual(checkTranspose, (CMR_DBLMAT*) components[1].transpose));
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
  CMRdblmatFree(cmr, &check);
  CMRdblmatFree(cmr, &checkTranspose);

  stringToDoubleMatrix(cmr, &check, "1 0 "
  );
  stringToDoubleMatrix(cmr, &checkTranspose, "0 1 "
  );
  ASSERT_TRUE(CMRdblmatCheckEqual(check, (CMR_DBLMAT*) components[2].matrix));
  ASSERT_TRUE(CMRdblmatCheckEqual(checkTranspose, (CMR_DBLMAT*) components[2].transpose));
  ASSERT_EQ(components[2].rowsToOriginal[0], 3);
  ASSERT_EQ(rowsToComponents[3], 2);
  ASSERT_EQ(rowsToComponentRows[3], 0);
  CMRdblmatFree(cmr, &check);
  CMRdblmatFree(cmr, &checkTranspose);

  stringToDoubleMatrix(cmr, &check, "1 1 "
    "1 "
  );
  stringToDoubleMatrix(cmr, &checkTranspose, "1 1 "
    "1 "
  );
  ASSERT_TRUE(CMRdblmatCheckEqual(check, (CMR_DBLMAT*) components[3].matrix));
  ASSERT_TRUE(CMRdblmatCheckEqual(checkTranspose, (CMR_DBLMAT*) components[3].transpose));
  ASSERT_EQ(components[3].rowsToOriginal[0], 4);
  ASSERT_EQ(components[3].columnsToOriginal[0], 0);
  ASSERT_EQ(rowsToComponents[4], 3);
  ASSERT_EQ(columnsToComponents[0], 3);
  ASSERT_EQ(rowsToComponentRows[4], 0);
  ASSERT_EQ(columnsToComponentColumns[0], 0);
  CMRdblmatFree(cmr, &check);
  CMRdblmatFree(cmr, &checkTranspose);

  stringToDoubleMatrix(cmr, &check, "3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  stringToDoubleMatrix(cmr, &checkTranspose, "3 3 "
    "1 2 0 "
    "0 3 0 "
    "0 4 5 "
  );
  ASSERT_TRUE(CMRdblmatCheckEqual(check, (CMR_DBLMAT*) components[4].matrix));
  ASSERT_TRUE(CMRdblmatCheckEqual(checkTranspose, (CMR_DBLMAT*) components[4].transpose));
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
  CMRdblmatFree(cmr, &check);
  CMRdblmatFree(cmr, &checkTranspose);

  stringToDoubleMatrix(cmr, &check, "0 1 "
  );
  stringToDoubleMatrix(cmr, &checkTranspose, "1 0 "
  );
  ASSERT_TRUE(CMRdblmatCheckEqual(check, (CMR_DBLMAT*) components[5].matrix));
  ASSERT_TRUE(CMRdblmatCheckEqual(checkTranspose, (CMR_DBLMAT*) components[5].transpose));
  ASSERT_EQ(components[5].columnsToOriginal[0], 9);
  ASSERT_EQ(columnsToComponents[9], 5);
  ASSERT_EQ(columnsToComponentColumns[9], 0);
  CMRdblmatFree(cmr, &check);
  CMRdblmatFree(cmr, &checkTranspose);

  for (size_t c = 0; c < numComponents; ++c)
  {
    CMRdblmatFree(cmr, (CMR_DBLMAT**) &components[c].matrix);
    CMRdblmatFree(cmr, (CMR_DBLMAT**) &components[c].transpose);
    CMRfreeBlockArray(cmr, &components[c].rowsToOriginal);
    CMRfreeBlockArray(cmr, &components[c].columnsToOriginal);
  }
  CMRfreeBlockArray(cmr, &components);

  CMRdblmatFree(cmr, &matrix);
  CMRfreeEnvironment(&cmr);
}

TEST(BlockDecomposition, IntToInt)
{
  CMR* cmr = NULL;
  CMRcreateEnvironment(&cmr);

  CMR_INTMAT* matrix = NULL;
  stringToIntMatrix(cmr, &matrix, "10 10 "
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
  CMR_BLOCK* components = NULL;
  size_t rowsToComponents[10];
  size_t columnsToComponents[10];
  size_t rowsToComponentRows[10];
  size_t columnsToComponentColumns[10];

  CMR_INTMAT* check = NULL;
  CMR_INTMAT* checkTranspose = NULL;
  ASSERT_CMR_CALL( CMRdecomposeBlocks(cmr, (CMR_MATRIX*) matrix, sizeof(int), sizeof(int), &numComponents, &components,
    rowsToComponents, columnsToComponents, rowsToComponentRows, columnsToComponentColumns) );

  ASSERT_EQ(numComponents, 6);
  stringToIntMatrix(cmr, &check, "3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  stringToIntMatrix(cmr, &checkTranspose, "3 3 "
    "1 3 0 "
    "2 4 5 "
    "0 0 6 "
  );
  ASSERT_TRUE(CMRintmatCheckEqual(check, (CMR_INTMAT*) components[0].matrix));
  ASSERT_TRUE(CMRintmatCheckEqual(checkTranspose, (CMR_INTMAT*) components[0].transpose));
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
  CMRintmatFree(cmr, &check);
  CMRintmatFree(cmr, &checkTranspose);

  stringToIntMatrix(cmr, &check, "2 2 "
    "1 0 "
    "3 2 "
  );
  stringToIntMatrix(cmr, &checkTranspose, "2 2 "
    "1 3 "
    "0 2 "
  );
  ASSERT_TRUE(CMRintmatCheckEqual(check, (CMR_INTMAT*) components[1].matrix));
  ASSERT_TRUE(CMRintmatCheckEqual(checkTranspose, (CMR_INTMAT*) components[1].transpose));
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
  CMRintmatFree(cmr, &check);
  CMRintmatFree(cmr, &checkTranspose);

  stringToIntMatrix(cmr, &check, "1 0 "
  );
  stringToIntMatrix(cmr, &checkTranspose, "0 1 "
  );
  ASSERT_TRUE(CMRintmatCheckEqual(check, (CMR_INTMAT*) components[2].matrix));
  ASSERT_TRUE(CMRintmatCheckEqual(checkTranspose, (CMR_INTMAT*) components[2].transpose));
  ASSERT_EQ(components[2].rowsToOriginal[0], 3);
  ASSERT_EQ(rowsToComponents[3], 2);
  ASSERT_EQ(rowsToComponentRows[3], 0);
  CMRintmatFree(cmr, &check);
  CMRintmatFree(cmr, &checkTranspose);

  stringToIntMatrix(cmr, &check, "1 1 "
    "1 "
  );
  stringToIntMatrix(cmr, &checkTranspose, "1 1 "
    "1 "
  );
  ASSERT_TRUE(CMRintmatCheckEqual(check, (CMR_INTMAT*) components[3].matrix));
  ASSERT_TRUE(CMRintmatCheckEqual(checkTranspose, (CMR_INTMAT*) components[3].transpose));
  ASSERT_EQ(components[3].rowsToOriginal[0], 4);
  ASSERT_EQ(components[3].columnsToOriginal[0], 0);
  ASSERT_EQ(rowsToComponents[4], 3);
  ASSERT_EQ(columnsToComponents[0], 3);
  ASSERT_EQ(rowsToComponentRows[4], 0);
  ASSERT_EQ(columnsToComponentColumns[0], 0);
  CMRintmatFree(cmr, &check);
  CMRintmatFree(cmr, &checkTranspose);

  stringToIntMatrix(cmr, &check, "3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  stringToIntMatrix(cmr, &checkTranspose, "3 3 "
    "1 2 0 "
    "0 3 0 "
    "0 4 5 "
  );
  ASSERT_TRUE(CMRintmatCheckEqual(check, (CMR_INTMAT*) components[4].matrix));
  ASSERT_TRUE(CMRintmatCheckEqual(checkTranspose, (CMR_INTMAT*) components[4].transpose));
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
  CMRintmatFree(cmr, &check);
  CMRintmatFree(cmr, &checkTranspose);

  stringToIntMatrix(cmr, &check, "0 1 "
  );
  stringToIntMatrix(cmr, &checkTranspose, "1 0 "
  );
  ASSERT_TRUE(CMRintmatCheckEqual(check, (CMR_INTMAT*) components[5].matrix));
  ASSERT_TRUE(CMRintmatCheckEqual(checkTranspose, (CMR_INTMAT*) components[5].transpose));
  ASSERT_EQ(components[5].columnsToOriginal[0], 9);
  ASSERT_EQ(columnsToComponents[9], 5);
  ASSERT_EQ(columnsToComponentColumns[9], 0);
  CMRintmatFree(cmr, &check);
  CMRintmatFree(cmr, &checkTranspose);

  for (size_t c = 0; c < numComponents; ++c)
  {
    CMRintmatFree(cmr, (CMR_INTMAT**) &components[c].matrix);
    CMRintmatFree(cmr, (CMR_INTMAT**) &components[c].transpose);
    CMRfreeBlockArray(cmr, &components[c].rowsToOriginal);
    CMRfreeBlockArray(cmr, &components[c].columnsToOriginal);
  }
  CMRfreeBlockArray(cmr, &components);

  CMRintmatFree(cmr, &matrix);
  CMRfreeEnvironment(&cmr);
}

TEST(BlockDecomposition, CharToChar)
{
  CMR* cmr = NULL;
  CMRcreateEnvironment(&cmr);

  CMR_CHRMAT* matrix = NULL;
  stringToCharMatrix(cmr, &matrix, "10 10 "
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
  CMR_BLOCK* components = NULL;
  size_t rowsToComponents[10];
  size_t columnsToComponents[10];
  size_t rowsToComponentRows[10];
  size_t columnsToComponentColumns[10];

  CMR_CHRMAT* check = NULL;
  CMR_CHRMAT* checkTranspose = NULL;
  ASSERT_CMR_CALL( CMRdecomposeBlocks(cmr, (CMR_MATRIX*) matrix, sizeof(char), sizeof(char), &numComponents, &components,
    rowsToComponents, columnsToComponents, rowsToComponentRows, columnsToComponentColumns) );

  ASSERT_EQ(numComponents, 6);
  stringToCharMatrix(cmr, &check, "3 3 "
    "1 2 0 "
    "3 4 0 "
    "0 5 6 "
  );
  stringToCharMatrix(cmr, &checkTranspose, "3 3 "
    "1 3 0 "
    "2 4 5 "
    "0 0 6 "
  );
  ASSERT_TRUE(CMRchrmatCheckEqual(check, (CMR_CHRMAT*) components[0].matrix));
  ASSERT_TRUE(CMRchrmatCheckEqual(checkTranspose, (CMR_CHRMAT*) components[0].transpose));
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
  CMRchrmatFree(cmr, &check);
  CMRchrmatFree(cmr, &checkTranspose);

  stringToCharMatrix(cmr, &check, "2 2 "
    "1 0 "
    "3 2 "
  );
  stringToCharMatrix(cmr, &checkTranspose, "2 2 "
    "1 3 "
    "0 2 "
  );
  ASSERT_TRUE(CMRchrmatCheckEqual(check, (CMR_CHRMAT*) components[1].matrix));
  ASSERT_TRUE(CMRchrmatCheckEqual(checkTranspose, (CMR_CHRMAT*) components[1].transpose));
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
  CMRchrmatFree(cmr, &check);
  CMRchrmatFree(cmr, &checkTranspose);

  stringToCharMatrix(cmr, &check, "1 0 "
  );
  stringToCharMatrix(cmr, &checkTranspose, "0 1 "
  );
  ASSERT_TRUE(CMRchrmatCheckEqual(check, (CMR_CHRMAT*) components[2].matrix));
  ASSERT_TRUE(CMRchrmatCheckEqual(checkTranspose, (CMR_CHRMAT*) components[2].transpose));
  ASSERT_EQ(components[2].rowsToOriginal[0], 3);
  ASSERT_EQ(rowsToComponents[3], 2);
  ASSERT_EQ(rowsToComponentRows[3], 0);
  CMRchrmatFree(cmr, &check);
  CMRchrmatFree(cmr, &checkTranspose);

  stringToCharMatrix(cmr, &check, "1 1 "
    "1 "
  );
  stringToCharMatrix(cmr, &checkTranspose, "1 1 "
    "1 "
  );
  ASSERT_TRUE(CMRchrmatCheckEqual(check, (CMR_CHRMAT*) components[3].matrix));
  ASSERT_TRUE(CMRchrmatCheckEqual(checkTranspose, (CMR_CHRMAT*) components[3].transpose));
  ASSERT_EQ(components[3].rowsToOriginal[0], 4);
  ASSERT_EQ(components[3].columnsToOriginal[0], 0);
  ASSERT_EQ(rowsToComponents[4], 3);
  ASSERT_EQ(columnsToComponents[0], 3);
  ASSERT_EQ(rowsToComponentRows[4], 0);
  ASSERT_EQ(columnsToComponentColumns[0], 0);
  CMRchrmatFree(cmr, &check);
  CMRchrmatFree(cmr, &checkTranspose);

  stringToCharMatrix(cmr, &check, "3 3 "
    "1 0 0 "
    "2 3 4 "
    "0 0 5 "
  );
  stringToCharMatrix(cmr, &checkTranspose, "3 3 "
    "1 2 0 "
    "0 3 0 "
    "0 4 5 "
  );
  ASSERT_TRUE(CMRchrmatCheckEqual(check, (CMR_CHRMAT*) components[4].matrix));
  ASSERT_TRUE(CMRchrmatCheckEqual(checkTranspose, (CMR_CHRMAT*) components[4].transpose));
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
  CMRchrmatFree(cmr, &check);
  CMRchrmatFree(cmr, &checkTranspose);

  stringToCharMatrix(cmr, &check, "0 1 "
  );
  stringToCharMatrix(cmr, &checkTranspose, "1 0 "
  );
  ASSERT_TRUE(CMRchrmatCheckEqual(check, (CMR_CHRMAT*) components[5].matrix));
  ASSERT_TRUE(CMRchrmatCheckEqual(checkTranspose, (CMR_CHRMAT*) components[5].transpose));
  ASSERT_EQ(components[5].columnsToOriginal[0], 9);
  ASSERT_EQ(columnsToComponents[9], 5);
  ASSERT_EQ(columnsToComponentColumns[9], 0);
  CMRchrmatFree(cmr, &check);
  CMRchrmatFree(cmr, &checkTranspose);

  for (size_t c = 0; c < numComponents; ++c)
  {
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].matrix);
    CMRchrmatFree(cmr, (CMR_CHRMAT**) &components[c].transpose);
    CMRfreeBlockArray(cmr, &components[c].rowsToOriginal);
    CMRfreeBlockArray(cmr, &components[c].columnsToOriginal);
  }
  CMRfreeBlockArray(cmr, &components);

  CMRchrmatFree(cmr, &matrix);
  CMRfreeEnvironment(&cmr);
}
