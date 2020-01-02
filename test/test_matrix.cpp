#include <gtest/gtest.h>

#include "common.h"
#include <tu/matrix.h>

TEST(Matrix, Submatrix)
{
  TU* tu;
  TU_MATRIX_CHAR matrix = stringToMatrixChar("10 10 "
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

  TU_SUBMATRIX* submatrix;
  TUcreateSubmatrix(&submatrix, 3, 3);
  submatrix->rows[0] = 1;
  submatrix->rows[1] = 3;
  submatrix->rows[2] = 4;
  submatrix->columns[0] = 1;
  submatrix->columns[1] = 4;
  submatrix->columns[2] = 6;

  TU_MATRIX_CHAR result;
  TUfilterSubmatrixChar(&matrix, submatrix, &result);

  TU_MATRIX_CHAR check = stringToMatrixChar("3 3 "
    "+1  0   0"
    " 0 -1  +1"
    " 0  0  -1"
  );
  ASSERT_TRUE(TUcheckMatrixEqualChar(&result, &check));
  TUclearMatrixChar(&check);

  TUclearMatrixChar(&result);
  TUfreeSubmatrix(&submatrix);
  TUclearMatrixChar(&matrix);
}
