#include <gtest/gtest.h>

#include <stdio.h>

#include "common.h"
#include <cmr/matrix.h>

TEST(Matrix, Read)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Double matrices. */
  {
    /* Zero matrix. */
    const char* denseInput = "3 3 "
      "0 0 0 "
      "0 0 0 "
      "0 0 0 ";
    FILE* stream = fmemopen((char*) denseInput, strlen(denseInput), "r");
    CMR_DBLMAT* dense = NULL;
    ASSERT_CMR_CALL( CMRdblmatCreateFromDenseStream(cmr, stream, &dense) );
    fclose(stream);
  
    const char* sparseInput = "3 3 0 ";
    stream = fmemopen((char*) sparseInput, strlen(sparseInput), "r");
    CMR_DBLMAT* sparse = NULL;
    ASSERT_CMR_CALL( CMRdblmatCreateFromSparseStream(cmr, stream, &sparse) );
    fclose(stream);

    ASSERT_TRUE( CMRdblmatCheckEqual(dense, sparse) );

    ASSERT_CMR_CALL( CMRdblmatFree(cmr, &dense) );
    ASSERT_CMR_CALL( CMRdblmatFree(cmr, &sparse) );
  }

  {
    /* Arbitrary matrix */
    const char* denseInput = "3 3 "
      "0.25 0 0 "
      "0 0.5 1 "
      "-1 0 2 ";
    FILE* stream = fmemopen((char*) denseInput, strlen(denseInput), "r");
    CMR_DBLMAT* dense = NULL;
    ASSERT_CMR_CALL( CMRdblmatCreateFromDenseStream(cmr, stream, &dense) );
    fclose(stream);

    const char* sparseInput = "3 3 5 "
      "1 1 0.25 "
      "2 2 0.5 "
      "2 3 1 "
      "3 1 -1 "
      "3 3 2 ";
    stream = fmemopen((char*) sparseInput, strlen(sparseInput), "r");
    CMR_DBLMAT* sparse = NULL;
    ASSERT_CMR_CALL( CMRdblmatCreateFromSparseStream(cmr, stream, &sparse) );
    fclose(stream);

    ASSERT_TRUE( CMRdblmatCheckEqual(dense, sparse) );

    ASSERT_CMR_CALL( CMRdblmatFree(cmr, &dense) );
    ASSERT_CMR_CALL( CMRdblmatFree(cmr, &sparse) );
  }

  /* Int matrices. */
  {
    /* Zero matrix. */
    const char* denseInput = "3 3 "
      "0 0 0 "
      "0 0 0 "
      "0 0 0 ";
    FILE* stream = fmemopen((char*) denseInput, strlen(denseInput), "r");
    CMR_INTMAT* dense = NULL;
    ASSERT_CMR_CALL( CMRintmatCreateFromDenseStream(cmr, stream, &dense) );
    fclose(stream);
  
    const char* sparseInput = "3 3 0 ";
    stream = fmemopen((char*) sparseInput, strlen(sparseInput), "r");
    CMR_INTMAT* sparse = NULL;
    ASSERT_CMR_CALL( CMRintmatCreateFromSparseStream(cmr, stream, &sparse) );
    fclose(stream);

    ASSERT_TRUE( CMRintmatCheckEqual(dense, sparse) );

    ASSERT_CMR_CALL( CMRintmatFree(cmr, &dense) );
    ASSERT_CMR_CALL( CMRintmatFree(cmr, &sparse) );
  }

  {
    /* Arbitrary matrix */
    const char* denseInput = "3 3 "
      "1 0 0 "
      "0 1 1 "
      "-1 0 2 ";
    FILE* stream = fmemopen((char*) denseInput, strlen(denseInput), "r");
    CMR_INTMAT* dense = NULL;
    ASSERT_CMR_CALL( CMRintmatCreateFromDenseStream(cmr, stream, &dense) );
    fclose(stream);

    const char* sparseInput = "3 3 5 "
      "1 1 1 "
      "2 2 1 "
      "2 3 1 "
      "3 1 -1 "
      "3 3 2 ";
    stream = fmemopen((char*) sparseInput, strlen(sparseInput), "r");
    CMR_INTMAT* sparse = NULL;
    ASSERT_CMR_CALL( CMRintmatCreateFromSparseStream(cmr, stream, &sparse) );
    fclose(stream);

    ASSERT_TRUE( CMRintmatCheckEqual(dense, sparse) );

    ASSERT_CMR_CALL( CMRintmatFree(cmr, &dense) );
    ASSERT_CMR_CALL( CMRintmatFree(cmr, &sparse) );
  }

  /* Char matrices. */
  {
    /* Zero matrix. */
    const char* denseInput = "3 3 "
      "0 0 0 "
      "0 0 0 "
      "0 0 0 ";
    FILE* stream = fmemopen((char*) denseInput, strlen(denseInput), "r");
    CMR_CHRMAT* dense = NULL;
    ASSERT_CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, stream, &dense) );
    fclose(stream);
  
    const char* sparseInput = "3 3 0 ";
    stream = fmemopen((char*) sparseInput, strlen(sparseInput), "r");
    CMR_CHRMAT* sparse = NULL;
    ASSERT_CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, stream, &sparse) );
    fclose(stream);

    ASSERT_TRUE( CMRchrmatCheckEqual(dense, sparse) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &dense) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &sparse) );
  }

  {
    /* Arbitrary matrix. */
    const char* denseInput = "3 3 "
      "1 0 0 "
      "0 1 1 "
      "-1 0 2 ";
    FILE* stream = fmemopen((char*) denseInput, strlen(denseInput), "r");
    CMR_CHRMAT* dense = NULL;
    ASSERT_CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, stream, &dense) );
    fclose(stream);

    const char* sparseInput = "3 3 5 "
      "1 1 1 "
      "2 2 1 "
      "2 3 1 "
      "3 1 -1 "
      "3 3 2 ";
    stream = fmemopen((char*) sparseInput, strlen(sparseInput), "r");
    CMR_CHRMAT* sparse = NULL;
    ASSERT_CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, stream, &sparse) );
    fclose(stream);

    ASSERT_TRUE( CMRchrmatCheckEqual(dense, sparse) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &dense) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &sparse) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Matrix, Transpose)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Double matrices. */
  {
    CMR_DBLMAT* A = NULL;
    stringToDoubleMatrix(cmr, &A, "4 5 "
      "1 2 3 0 0 "
      "0 4 5 0 6 "
      "7 0 0 8 0 "
      "0 0 9 0 0 "
    );

    CMR_DBLMAT* B = NULL;
    stringToDoubleMatrix(cmr, &B, "5 4 "
      "1 0 7 0 "
      "2 4 0 0 "
      "3 5 0 9 "
      "0 0 8 0 "
      "0 6 0 0 "
    );

    bool transposes;
    ASSERT_CMR_CALL( CMRdblmatCheckTranspose(cmr, A, B, &transposes) );
    ASSERT_TRUE(transposes);

    CMRdblmatFree(cmr, &B);
    CMRdblmatFree(cmr, &A);
  }

  /* Int matrices. */
  {
    CMR_INTMAT* A = NULL;
    stringToIntMatrix(cmr, &A, "4 5 "
      "1 2 3 0 0 "
      "0 4 5 0 6 "
      "7 0 0 8 0 "
      "0 0 9 0 0 "
    );

    CMR_INTMAT* B = NULL;
    stringToIntMatrix(cmr, &B, "5 4 "
      "1 0 7 0 "
      "2 4 0 0 "
      "3 5 0 9 "
      "0 0 8 0 "
      "0 6 0 0 "
    );

    bool transposes;
    ASSERT_CMR_CALL( CMRintmatCheckTranspose(cmr, A, B, &transposes) );
    ASSERT_TRUE(transposes);

    CMRintmatFree(cmr, &B);
    CMRintmatFree(cmr, &A);
  }

  /* Char matrices. */
  {
    CMR_CHRMAT* A = NULL;
    stringToCharMatrix(cmr, &A, "4 5 "
      "1 2 3 0 0 "
      "0 4 5 0 6 "
      "7 0 0 8 0 "
      "0 0 9 0 0 "
    );

    CMR_CHRMAT* B = NULL;
    stringToCharMatrix(cmr, &B, "5 4 "
      "1 0 7 0 "
      "2 4 0 0 "
      "3 5 0 9 "
      "0 0 8 0 "
      "0 6 0 0 "
    );

    bool transposes;
    ASSERT_CMR_CALL( CMRchrmatCheckTranspose(cmr, A, B, &transposes) );
    ASSERT_TRUE(transposes);

    CMRchrmatFree(cmr, &B);
    CMRchrmatFree(cmr, &A);
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Matrix, Submatrix)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  stringToCharMatrix(cmr, &matrix, "10 10 "
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

  CMR_SUBMAT* submatrix = NULL;
  ASSERT_CMR_CALL( CMRsubmatCreate(cmr, 3, 3, &submatrix) );
  submatrix->rows[0] = 1;
  submatrix->rows[1] = 3;
  submatrix->rows[2] = 4;
  submatrix->columns[0] = 1;
  submatrix->columns[1] = 4;
  submatrix->columns[2] = 6;

  CMR_CHRMAT* result = NULL;
  ASSERT_CMR_CALL( CMRchrmatSlice(cmr, matrix, submatrix, &result) );

  CMRchrmatPrintDense(cmr, result, stdout, '0', true);
  
  CMR_CHRMAT* check = NULL;
  stringToCharMatrix(cmr, &check, "3 3 "
    "+1  0   0"
    " 0 -1  +1"
    " 0  0  -1"
  );
  CMRchrmatPrintDense(cmr, check, stdout, '0', true);
  ASSERT_TRUE(CMRchrmatCheckEqual(result, check));
  CMRchrmatFree(cmr, &check);

  CMRchrmatFree(cmr, &result);
  CMRsubmatFree(cmr, &submatrix);
  CMRchrmatFree(cmr, &matrix);

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
