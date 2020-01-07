#include <gtest/gtest.h>

#include "common.h"
#include <tu/matrix.h>

TEST(Matrix, Transpose)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);

  /* Double matrices. */
  {
    TU_DOUBLE_MATRIX* A = NULL;
    stringToDoubleMatrix(tu, &A, "4 5 "
      "1 2 3 0 0 "
      "0 4 5 0 6 "
      "7 0 0 8 0 "
      "0 0 9 0 0 "
    );

    TU_DOUBLE_MATRIX* B = NULL;
    stringToDoubleMatrix(tu, &B, "5 4 "
      "1 0 7 0 "
      "2 4 0 0 "
      "3 5 0 9 "
      "0 0 8 0 "
      "0 6 0 0 "
    );

    ASSERT_TRUE(TUcheckDoubleMatrixTranspose(A, B));

    TUfreeDoubleMatrix(tu, &B);
    TUfreeDoubleMatrix(tu, &A);
  }

  /* Int matrices. */
  {
    TU_INT_MATRIX* A = NULL;
    stringToIntMatrix(tu, &A, "4 5 "
      "1 2 3 0 0 "
      "0 4 5 0 6 "
      "7 0 0 8 0 "
      "0 0 9 0 0 "
    );

    TU_INT_MATRIX* B = NULL;
    stringToIntMatrix(tu, &B, "5 4 "
      "1 0 7 0 "
      "2 4 0 0 "
      "3 5 0 9 "
      "0 0 8 0 "
      "0 6 0 0 "
    );

    ASSERT_TRUE(TUcheckIntMatrixTranspose(A, B));

    TUfreeIntMatrix(tu, &B);
    TUfreeIntMatrix(tu, &A);
  }

  /* Char matrices. */
  {
    TU_CHAR_MATRIX* A = NULL;
    stringToCharMatrix(tu, &A, "4 5 "
      "1 2 3 0 0 "
      "0 4 5 0 6 "
      "7 0 0 8 0 "
      "0 0 9 0 0 "
    );

    TU_CHAR_MATRIX* B = NULL;
    stringToCharMatrix(tu, &B, "5 4 "
      "1 0 7 0 "
      "2 4 0 0 "
      "3 5 0 9 "
      "0 0 8 0 "
      "0 6 0 0 "
    );

    ASSERT_TRUE(TUcheckCharMatrixTranspose(A, B));

    TUfreeCharMatrix(tu, &B);
    TUfreeCharMatrix(tu, &A);
  }

  TUfreeEnvironment(&tu);
}

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
