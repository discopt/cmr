#include <gtest/gtest.h>

#include "common.h"
#include <tu/matrix.h>

TEST(Matrix, Transpose)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);

  /* Double matrices. */
  {
    TU_DBLMAT* A = NULL;
    stringToDoubleMatrix(tu, &A, "4 5 "
      "1 2 3 0 0 "
      "0 4 5 0 6 "
      "7 0 0 8 0 "
      "0 0 9 0 0 "
    );

    TU_DBLMAT* B = NULL;
    stringToDoubleMatrix(tu, &B, "5 4 "
      "1 0 7 0 "
      "2 4 0 0 "
      "3 5 0 9 "
      "0 0 8 0 "
      "0 6 0 0 "
    );

    ASSERT_TRUE(TUdblmatCheckTranspose(A, B));

    TUdblmatFree(tu, &B);
    TUdblmatFree(tu, &A);
  }

  /* Int matrices. */
  {
    TU_INTMAT* A = NULL;
    stringToIntMatrix(tu, &A, "4 5 "
      "1 2 3 0 0 "
      "0 4 5 0 6 "
      "7 0 0 8 0 "
      "0 0 9 0 0 "
    );

    TU_INTMAT* B = NULL;
    stringToIntMatrix(tu, &B, "5 4 "
      "1 0 7 0 "
      "2 4 0 0 "
      "3 5 0 9 "
      "0 0 8 0 "
      "0 6 0 0 "
    );

    ASSERT_TRUE(TUintmatCheckTranspose(A, B));

    TUintmatFree(tu, &B);
    TUintmatFree(tu, &A);
  }

  /* Char matrices. */
  {
    TU_CHRMAT* A = NULL;
    stringToCharMatrix(tu, &A, "4 5 "
      "1 2 3 0 0 "
      "0 4 5 0 6 "
      "7 0 0 8 0 "
      "0 0 9 0 0 "
    );

    TU_CHRMAT* B = NULL;
    stringToCharMatrix(tu, &B, "5 4 "
      "1 0 7 0 "
      "2 4 0 0 "
      "3 5 0 9 "
      "0 0 8 0 "
      "0 6 0 0 "
    );

    ASSERT_TRUE(TUchrmatCheckTranspose(A, B));

    TUchrmatFree(tu, &B);
    TUchrmatFree(tu, &A);
  }

  TUfreeEnvironment(&tu);
}

TEST(Matrix, Submatrix)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);

  TU_CHRMAT* matrix = NULL;
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

  TU_SUBMAT* submatrix = NULL;
  TUsubmatCreate(tu, &submatrix, 3, 3);
  submatrix->rows[0] = 1;
  submatrix->rows[1] = 3;
  submatrix->rows[2] = 4;
  submatrix->columns[0] = 1;
  submatrix->columns[1] = 4;
  submatrix->columns[2] = 6;

  TU_CHRMAT* result = NULL;
  ASSERT_TU_CALL( TUchrmatFilterSubmat(tu, matrix, submatrix, &result) );

  TU_CHRMAT* check = NULL;
  stringToCharMatrix(tu, &check, "3 3 "
    "+1  0   0"
    " 0 -1  +1"
    " 0  0  -1"
  );
  ASSERT_TRUE(TUchrmatCheckEqual(result, check));
  TUchrmatFree(tu, &check);

  TUchrmatFree(tu, &result);
  TUsubmatFree(tu, &submatrix);
  TUchrmatFree(tu, &matrix);

  TUfreeEnvironment(&tu);
}
