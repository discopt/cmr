#include <gtest/gtest.h>

#include <stdio.h>

#include "common.h"
#include <tu/matrix.h>

TEST(Matrix, Read)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  /* Double matrices. */
  {
    const char* denseInput = "3 3 "
      "1 0 0 "
      "0 1 1 "
      "-1 0 2 ";
    FILE* stream = fmemopen((char*) denseInput, strlen(denseInput), "r");
    TU_DBLMAT* dense = NULL;
    ASSERT_TU_CALL( TUdblmatCreateFromDenseStream(tu, &dense, stream) );
    fclose(stream);

    const char* sparseInput = "3 3 5 "
      "0 0 1 "
      "1 1 1 "
      "1 2 1 "
      "2 0 -1 "
      "2 2 2 ";
    stream = fmemopen((char*) sparseInput, strlen(sparseInput), "r");
    TU_DBLMAT* sparse = NULL;
    ASSERT_TU_CALL( TUdblmatCreateFromSparseStream(tu, &sparse, stream) );
    fclose(stream);

    ASSERT_TRUE( TUdblmatCheckEqual(dense, sparse) );

    ASSERT_TU_CALL( TUdblmatFree(tu, &dense) );
    ASSERT_TU_CALL( TUdblmatFree(tu, &sparse) );
  }

  /* Int matrices. */
  {
    const char* denseInput = "3 3 "
      "1 0 0 "
      "0 1 1 "
      "-1 0 2 ";
    FILE* stream = fmemopen((char*) denseInput, strlen(denseInput), "r");
    TU_INTMAT* dense = NULL;
    ASSERT_TU_CALL( TUintmatCreateFromDenseStream(tu, &dense, stream) );
    fclose(stream);

    const char* sparseInput = "3 3 5 "
      "0 0 1 "
      "1 1 1 "
      "1 2 1 "
      "2 0 -1 "
      "2 2 2 ";
    stream = fmemopen((char*) sparseInput, strlen(sparseInput), "r");
    TU_INTMAT* sparse = NULL;
    ASSERT_TU_CALL( TUintmatCreateFromSparseStream(tu, &sparse, stream) );
    fclose(stream);

    ASSERT_TRUE( TUintmatCheckEqual(dense, sparse) );

    ASSERT_TU_CALL( TUintmatFree(tu, &dense) );
    ASSERT_TU_CALL( TUintmatFree(tu, &sparse) );
  }

  /* Char matrices. */
  {
    const char* denseInput = "3 3 "
      "1 0 0 "
      "0 1 1 "
      "-1 0 2 ";
    FILE* stream = fmemopen((char*) denseInput, strlen(denseInput), "r");
    TU_CHRMAT* dense = NULL;
    ASSERT_TU_CALL( TUchrmatCreateFromDenseStream(tu, &dense, stream) );
    fclose(stream);

    const char* sparseInput = "3 3 5 "
      "0 0 1 "
      "1 1 1 "
      "1 2 1 "
      "2 0 -1 "
      "2 2 2 ";
    stream = fmemopen((char*) sparseInput, strlen(sparseInput), "r");
    TU_CHRMAT* sparse = NULL;
    ASSERT_TU_CALL( TUchrmatCreateFromSparseStream(tu, &sparse, stream) );
    fclose(stream);

    ASSERT_TRUE( TUchrmatCheckEqual(dense, sparse) );

    ASSERT_TU_CALL( TUchrmatFree(tu, &dense) );
    ASSERT_TU_CALL( TUchrmatFree(tu, &sparse) );
  }

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Matrix, Transpose)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

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

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Matrix, Submatrix)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

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

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}
