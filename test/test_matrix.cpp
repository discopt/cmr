#include <gtest/gtest.h>

#include "common.h"
#include <tu/matrix.h>

TEST(Matrix, Submatrix)
{
  TU* tu = NULL;
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

  TU_SUBMATRIX* submatrix = NULL;
  TUcreateSubmatrix(tu, &submatrix, 3, 3);
  submatrix->rows[0] = 1;
  submatrix->rows[1] = 3;
  submatrix->rows[2] = 4;
  submatrix->columns[0] = 1;
  submatrix->columns[1] = 4;
  submatrix->columns[2] = 6;

  TU_CHAR_MATRIX* result = NULL;
  TUfilterCharSubmatrix(tu, matrix, submatrix, &result);

  TU_CHAR_MATRIX* check = NULL;
  stringToCharMatrix(tu, &check, "3 3 "
    "+1  0   0"
    " 0 -1  +1"
    " 0  0  -1"
  );
  ASSERT_TRUE(TUcheckCharMatrixEqual(result, check));
  TUfreeCharMatrix(tu, &check);

  TUfreeCharMatrix(tu, &result);
  TUfreeSubmatrix(tu, &submatrix);
  TUfreeCharMatrix(tu, &matrix);

  TUfreeEnvironment(&tu);
}
