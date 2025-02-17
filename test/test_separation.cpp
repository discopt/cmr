#include <gtest/gtest.h>

#include "common.h"
#include <cmr/separation.h>

TEST(Separation, Onesum)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &first, "5 5 "
      " 1  1  0  0  0 "
      " 1  0  1 -1  1 "
      " 0 -1  1  0 -1 "
      " 0  0 -1  1  0 "
      " 0  1  1  0  1 "
    ) );
    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &second, "5 5 "
      " 1 -1  1  0  0 "
      " 1  1  1  1 -1 "
      " 0  0 -1  0  1 "
      " 1  0  0 -1  0 "
      " 0  1  0  0  1 "
    ) );
    CMR_CHRMAT* third = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &third, "2 1 "
      "  1 "
      " -1 "
    ) );

    CMR_CHRMAT* check = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &check, "12 11 "
      " 1  1  0  0  0  0  0  0  0  0  0 "
      " 1  0  1 -1  1  0  0  0  0  0  0 "
      " 0 -1  1  0 -1  0  0  0  0  0  0 "
      " 0  0 -1  1  0  0  0  0  0  0  0 "
      " 0  1  1  0  1  0  0  0  0  0  0 "
      " 0  0  0  0  0  1 -1  1  0  0  0 "
      " 0  0  0  0  0  1  1  1  1 -1  0 "
      " 0  0  0  0  0  0  0 -1  0  1  0 "
      " 0  0  0  0  0  1  0  0 -1  0  0 "
      " 0  0  0  0  0  0  1  0  0  1  0 "
      " 0  0  0  0  0  0  0  0  0  0  1 "
      " 0  0  0  0  0  0  0  0  0  0 -1 "
    ) );

    CMR_CHRMAT* onesum = NULL;
    CMR_CHRMAT* matrices[3] = { first, second, third };
    ASSERT_CMR_CALL( CMRonesumCompose(cmr, 3, matrices, &onesum) );

    ASSERT_TRUE( CMRchrmatCheckEqual(onesum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &onesum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &third) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );  
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Separation, TwosumComposition)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &first, "5 5 "
      " 1  1  0  0  0 "
      " 1  0  1 -1  1 " /* marker row */
      " 0 -1  1  0 -1 "
      " 0  0 -1  1  0 "
      " 0  1  1  0  1 "
    ) );
    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &second, "5 5 "
      " 1 -1  1  0  0 "
      " 1  1  1  1 -1 "
      " 0  0 -1  0  1 "
      " 1  0  0 -1  0 "
      " 0  1  0  0  1 "
      /*      ^
       * marker column */
    ) );

    CMR_CHRMAT* check = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &check, "9 9 "
      " 1  1  0  0  0  0  0  0  0 "
      " 0 -1  1  0 -1  0  0  0  0 "
      " 0  0 -1  1  0  0  0  0  0 "
      " 0  1  1  0  1  0  0  0  0 "
      " 1  0  1 -1  1  1 -1  0  0 "
      " 1  0  1 -1  1  1  1  1 -1 "
      "-1  0 -1  1 -1  0  0  0  1 "
      " 0  0  0  0  0  1  0 -1  0 "
      " 0  0  0  0  0  0  1  0  1 "
    ) );

    CMR_CHRMAT* twosum = NULL;
    size_t firstSpecialRow = 1;
    size_t secondSpecialColumn = 2;
    ASSERT_CMR_CALL( CMRtwosumCompose(cmr, first, second, &firstSpecialRow, NULL, NULL, &secondSpecialColumn, 3,
      &twosum) );

    CMRchrmatPrintDense(cmr, twosum, stdout, '0', false);

    ASSERT_TRUE( CMRchrmatCheckEqual(twosum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &twosum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
  }

  {
    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &first, "5 5 "
      " 1  1  0  0  0 "
      " 1  0  1 -1  1 "
      " 0 -1  1  0 -1 "
      " 0  0 -1  1  0 "
      " 0  1  1  0  1 "
      /*            ^
       * marker column */
    ) );
    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &second, "5 5 "
      " 1 -1  1  0  0 " /* marker row */
      " 1  1  1  1 -1 "
      " 0  0 -1  0  1 "
      " 1  0  0 -1  0 "
      " 0  1  0  0  1 " 
    ) );

    CMR_CHRMAT* check = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &check, "9 9 "
      " 1  1  0  0  0  0  0  0  0 "
      " 1  0  1 -1  1 -1  1  0  0 "
      " 0 -1  1  0 -1  1 -1  0  0 "
      " 0  0 -1  1  0  0  0  0  0 "
      " 0  1  1  0  1 -1  1  0  0 "
      " 0  0  0  0  1  1  1  1 -1 "
      " 0  0  0  0  0  0 -1  0  1 "
      " 0  0  0  0  1  0  0 -1  0 "
      " 0  0  0  0  0  1  0  0  1 " 
    ) );

    CMR_CHRMAT* twosum = NULL;
    size_t firstSpecialColumn = 4;
    size_t secondSpecialRow = 0;
    ASSERT_CMR_CALL( CMRtwosumCompose(cmr, first, second, NULL, &firstSpecialColumn, &secondSpecialRow, NULL, 3,
      &twosum) );

    CMRchrmatPrintDense(cmr, twosum, stdout, '0', false);

    ASSERT_TRUE( CMRchrmatCheckEqual(twosum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &twosum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Separation, TwosumDecomposition)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "9 9 "
      " 1  1  0  0  0  0  0  0  0 "
      " 1  0  1 -1 -1  1  1  0  0 "
      " 0 -1  1  0  1 -1 -1  0  0 "
      " 0  0 -1  1  0  0  0  0  0 "
      " 0  1  1  0 -1  1  1  0  0 "
      " 0  0  0  0  1  1  1  1 -1 "
      " 0  0  0  0  0  0 -1  0  1 "
      " 0  0  0  0  1  0  0 -1  0 "
      " 0  0  0  0  0  1  0  0  1 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', false);

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 9, 9, &sepa) );
    sepa->type = CMR_SEPA_TYPE_TWO;

    sepa->rowsFlags[0] = CMR_SEPA_FIRST;
    sepa->rowsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[2] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[3] = CMR_SEPA_FIRST;
    sepa->rowsFlags[4] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[5] = CMR_SEPA_SECOND;
    sepa->rowsFlags[6] = CMR_SEPA_SECOND;
    sepa->rowsFlags[7] = CMR_SEPA_SECOND;
    sepa->rowsFlags[8] = CMR_SEPA_SECOND;

    sepa->columnsFlags[0] = CMR_SEPA_FIRST;
    sepa->columnsFlags[1] = CMR_SEPA_FIRST;
    sepa->columnsFlags[2] = CMR_SEPA_FIRST;
    sepa->columnsFlags[3] = CMR_SEPA_FIRST;
    sepa->columnsFlags[4] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[5] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[6] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[7] = CMR_SEPA_SECOND;
    sepa->columnsFlags[8] = CMR_SEPA_SECOND;

    CMR_CHRMAT* first = NULL;
    size_t firstSpecialRow;
    size_t firstSpecialColumn;
    ASSERT_CMR_CALL( CMRtwosumDecomposeFirst(cmr, matrix, sepa, &first, NULL, NULL, NULL, NULL, &firstSpecialRow,
      &firstSpecialColumn) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, "5 5 "
      " 1  1  0  0  0 "
      " 1  0  1 -1 -1 "
      " 0 -1  1  0  1 "
      " 0  0 -1  1  0 "
      " 0  1  1  0 -1 "
    ) );
    CMRchrmatPrintDense(cmr, first, stdout, '0', false);
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    size_t secondSpecialRow;
    size_t secondSpecialColumn;
    ASSERT_CMR_CALL( CMRtwosumDecomposeSecond(cmr, matrix, sepa, &second, NULL, NULL, NULL, NULL, &secondSpecialRow,
      &secondSpecialColumn) );
    CMR_CHRMAT* checkSecond = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkSecond, "5 5 "
      " 1 -1 -1  0  0 "
      " 1  1  1  1 -1 "
      " 0  0 -1  0  1 "
      " 1  0  0 -1  0 "
      " 0  1  0  0  1 "
    ) );
    CMRchrmatPrintDense(cmr, second, stdout, '0', false);
    ASSERT_TRUE( CMRchrmatCheckEqual(second, checkSecond) );

    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );

    /* Compose again. */
    CMR_CHRMAT* check = NULL;
    ASSERT_CMR_CALL( CMRtwosumCompose(cmr, first, second, &firstSpecialRow, &firstSpecialColumn, &secondSpecialRow,
      &secondSpecialColumn, 0, &check) );

    ASSERT_TRUE( CMRchrmatCheckEqual(matrix, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkFirst) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkSecond) );
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "9 9 "
      " 1  1  0  0  0  0  0  0  0 "
      " 1  0  1 -1 -1  1  1  0  0 "
      " 0 -1  1  0  1 -1 -1  0  0 "
      " 0  0 -1  1  0  0  0  0  0 "
      " 0  1  1  0 -1  1  1  0  0 "
      " 0  0  0  0  1  1  1  1 -1 "
      " 0  0  0  0  0  0 -1  0  1 "
      " 0  0  0  0  1  0  0 -1  0 "
      " 0  0  0  0  0  1  0  0  1 "
    ) );

    CMRchrmatPrintDense(cmr, matrix, stdout, '0', false);

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 9, 9, &sepa) );
    sepa->type = CMR_SEPA_TYPE_TWO;

    sepa->rowsFlags[0] = CMR_SEPA_SECOND;
    sepa->rowsFlags[1] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[3] = CMR_SEPA_SECOND;
    sepa->rowsFlags[4] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[5] = CMR_SEPA_FIRST;
    sepa->rowsFlags[6] = CMR_SEPA_FIRST;
    sepa->rowsFlags[7] = CMR_SEPA_FIRST;
    sepa->rowsFlags[8] = CMR_SEPA_FIRST;

    sepa->columnsFlags[0] = CMR_SEPA_SECOND;
    sepa->columnsFlags[1] = CMR_SEPA_SECOND;
    sepa->columnsFlags[2] = CMR_SEPA_SECOND;
    sepa->columnsFlags[3] = CMR_SEPA_SECOND;
    sepa->columnsFlags[4] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[5] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[6] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[7] = CMR_SEPA_FIRST;
    sepa->columnsFlags[8] = CMR_SEPA_FIRST;

    CMR_CHRMAT* first = NULL;
    size_t firstSpecialRow;
    size_t firstSpecialColumn;
    ASSERT_CMR_CALL( CMRtwosumDecomposeFirst(cmr, matrix, sepa, &first, NULL, NULL, NULL, NULL, &firstSpecialRow,
      &firstSpecialColumn) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, "5 5 "
      "  1 1  1  1 -1 "
      "  0 0 -1  0  1 "
      "  1 0  0 -1  0 "
      "  0 1  0  0  1 "
      " -1 1  1  0  0 "
    ) );
    CMRchrmatPrintDense(cmr, first, stdout, '0', false);
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    size_t secondSpecialRow;
    size_t secondSpecialColumn;
    ASSERT_CMR_CALL( CMRtwosumDecomposeSecond(cmr, matrix, sepa, &second, NULL, NULL, NULL, NULL, &secondSpecialRow,
      &secondSpecialColumn) );
    CMR_CHRMAT* checkSecond = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkSecond, "5 5 "
      " 0  1  1  0  0 "
      " 1  1  0  1 -1 "
      "-1  0 -1  1  0 "
      " 0  0  0 -1  1 "
      " 1  0  1  1  0 "
    ) );
    CMRchrmatPrintDense(cmr, second, stdout, '0', false);
    ASSERT_TRUE( CMRchrmatCheckEqual(second, checkSecond) );

    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );

    /* Compose again. */
    CMR_CHRMAT* check = NULL;
    ASSERT_CMR_CALL( CMRtwosumCompose(cmr, second, first, &secondSpecialRow, &secondSpecialColumn, &firstSpecialRow,
      &firstSpecialColumn, 0, &check) );

    ASSERT_TRUE( CMRchrmatCheckEqual(matrix, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkFirst) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkSecond) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Separation, DeltasumComposition)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &first, "3 4 "
      " 1  1  0  0 "
      " 1  0  1  1 "
      " 0  1  0 -1 "
    ) );
    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &second, "3 4 "
      "-1  0  1  0 "
      " 1  1  0  1 "
      " 0  0  1  1 "
    ) );

    CMR_CHRMAT* check = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &check, " 4 4 "
      " 1  1  0  0 "
      " 1  0  1  0 "
      " 0  1  0  1 "
      " 0  0  1  1 "
    ) );

    CMR_CHRMAT* threesum = NULL;
    size_t specials[6] = { 2, 2, 3, 0, 0, 1 };
    ASSERT_CMR_CALL( CMRdeltasumCompose(cmr, first, second, &specials[0], &specials[1], &specials[3], &specials[4], 3,
      &threesum) );

    CMRchrmatPrintDense(cmr, threesum, stdout, '0', false);

    ASSERT_TRUE( CMRchrmatCheckEqual(threesum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &threesum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
  }

  {
    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &first, "3 4 "
      " 1  1  0  0 "
      " 1  0  1  1 "
      " 0  1  0  1 "
    ) );
    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &second, "3 4 "
      " 1  0  1  0 "
      " 1  1  0  1 "
      " 0  0  1  1 "
    ) );

    CMR_CHRMAT* check = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &check, " 4 4 "
      " 1  1  0  0 "
      " 1  0  1  0 "
      " 0  1  0  1 "
      " 0  0  1  1 "
    ) );

    CMR_CHRMAT* threesum = NULL;
    size_t specials[6] = { 2, 2, 3, 0, 0, 1 };
    ASSERT_CMR_CALL( CMRdeltasumCompose(cmr, first, second, &specials[0], &specials[1], &specials[3], &specials[4], 3,
      &threesum) );

    CMRchrmatPrintDense(cmr, threesum, stdout, '0', false);

    ASSERT_TRUE( CMRchrmatCheckEqual(threesum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &threesum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
  }

  /* Tests cases where special rows/columns are spread a bit. */
  {
    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &first, "5 6 "
      " 1  1  0  0  0  0 "
      " 1  0  0  1 -1  1 "
      " 0 -1  1  1  0 -1 "
      " 0  0 -1 -1  1  0 "
      " 0  1  1  1  0  1 "
    ) );
    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &second, "5 6 "
      " 1 -1 -1  1  0  0 "
      " 1  1  0  1  1 -1 "
      " 0  0  0 -1  0  1 "
      " 1  0  0  0 -1  0 "
      " 0  1  1  0  0  1 "
    ) );

    CMR_CHRMAT* check = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &check, " 8 8 "
      " 1  1  0  0  0  0  0  0 "
      " 0 -1  0 -1  1  1  1 -1 "
      " 0  0  1  0 -1 -1 -1  1 "
      " 0  1  0  1  1  1  1 -1 "
      "-1  0  1 -1  1  1  0  0 "
      " 0  0  0  0  0 -1  0  1 "
      " 0  0  0  0  1  0 -1  0 "
      " 1  0 -1  1  0  0  0  1 "
    ) );

    CMR_CHRMAT* threesum = NULL;
    size_t specials[6] = { 1, 2, 3, 1, 1, 2 };
    ASSERT_CMR_CALL( CMRdeltasumCompose(cmr, first, second, &specials[0], &specials[1], &specials[3], &specials[4], 3,
      &threesum) );

    CMRchrmatPrintDense(cmr, threesum, stdout, '0', false);
    CMRchrmatPrintDense(cmr, check, stdout, '0', false);

    ASSERT_TRUE( CMRchrmatCheckEqual(threesum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &threesum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
  }

  /* Tests cases where special rows/columns are spread a bit. */
  {
    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &first, "5 6 "
      " 1  1  0  0  0  0 "
      " 1  0  1 -1  1  0 "
      " 0 -1  1  0 -1  1 "
      " 0  0 -1 -1  0 -1 "
      " 0  1  1  0  1  1 "
    ) );
    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &second, "5 6 "
      " 1 -1 -1  0  0 -1 "
      " 1  1  1  1 -1  0 "
      " 0  0 -1  0  1  0 "
      " 1  0  0 -1  0  0 "
      " 0  1  0  0  1  1 "
    ) );

    CMR_CHRMAT* check = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &check, " 8 8 "
      " 1  1  0  0  0  0  0  0 "
      " 0 -1  0 -1  1  1  1 -1 "
      " 0  0 -1  0 -1 -1 -1  1 "
      " 0  1  0  1  1  1  1 -1 "
      "-1  0  1 -1  1 -1  0  0 "
      " 0  0  0  0  0 -1  0  1 "
      " 0  0  0  0  1  0 -1  0 "
      " 1  0 -1  1  0  0  0  1 "
    ) );

    CMR_CHRMAT* threesum = NULL;
    size_t specials[6] = { 1, 5, 2, 1, 1, 5 };
    ASSERT_CMR_CALL( CMRdeltasumCompose(cmr, first, second, &specials[0], &specials[1], &specials[3], &specials[4], 3,
      &threesum) );

    printf("Computed 3-sum:\n");
    CMRchrmatPrintDense(cmr, threesum, stdout, '0', false);
    printf("Expected result:\n");
    CMRchrmatPrintDense(cmr, check, stdout, '0', false);

    ASSERT_TRUE( CMRchrmatCheckEqual(threesum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &threesum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Separation, DeltasumDecomposition)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 4 4 "
      " 1  1  0  0 "
      " 1  0  1  0 "
      " 0  1  0  1 "
      " 0  0  1  1 "
    ) );
    CMR_CHRMAT* transpose = NULL;
    ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 4, 4, &sepa) );
    sepa->type = CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS;

    sepa->rowsFlags[0] = CMR_SEPA_FIRST;
    sepa->rowsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[3] = CMR_SEPA_SECOND;

    sepa->columnsFlags[0] = CMR_SEPA_FIRST;
    sepa->columnsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[3] = CMR_SEPA_SECOND;

    char epsilon;
    ASSERT_CMR_CALL( CMRdeltasumDecomposeEpsilon(cmr, matrix, transpose, sepa, &epsilon) );
    ASSERT_EQ(epsilon, -1);

    CMR_CHRMAT* first = NULL;
    size_t specials[3];
    ASSERT_CMR_CALL( CMRdeltasumDecomposeFirst(cmr, matrix, sepa, epsilon, &first, NULL, NULL, NULL, NULL, &specials[0],
      &specials[1]) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 3 4 "
      " 1  1  0  0 "
      " 1  0  1  1 "
      " 0  1  0 -1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( CMRdeltasumDecomposeSecond(cmr, matrix, sepa, epsilon, &second, NULL, NULL, NULL, NULL, NULL,
      NULL) );

    CMR_CHRMAT* checkSecond = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkSecond, " 3 4 "
      "-1  0  1  0 "
      " 1  1  0  1 "
      " 0  0  1  1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(second, checkSecond) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkSecond) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkFirst) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 4 4 "
      " 1  1  0  0 "
      " 1  0 -1  0 "
      " 0  1  0  1 "
      " 0  0  1  1 "
    ) );
    CMR_CHRMAT* transpose = NULL;
    ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 4, 4, &sepa) );
    sepa->type = CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS;

    sepa->rowsFlags[0] = CMR_SEPA_FIRST;
    sepa->rowsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[3] = CMR_SEPA_SECOND;

    sepa->columnsFlags[0] = CMR_SEPA_FIRST;
    sepa->columnsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[3] = CMR_SEPA_SECOND;

    char epsilon;
    ASSERT_CMR_CALL( CMRdeltasumDecomposeEpsilon(cmr, matrix, transpose, sepa, &epsilon) );
    ASSERT_EQ(epsilon, +1);

    CMR_CHRMAT* first = NULL;
    size_t specials[3];
    ASSERT_CMR_CALL( CMRdeltasumDecomposeFirst(cmr, matrix, sepa, epsilon, &first, NULL, NULL, NULL, NULL, &specials[0],
      &specials[1]) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 3 4 "
      " 1  1  0  0 "
      " 1  0 -1 -1 "
      " 0  1  0  1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( CMRdeltasumDecomposeSecond(cmr, matrix, sepa, epsilon, &second, NULL, NULL, NULL, NULL,
      &specials[0], &specials[1]) );

    CMR_CHRMAT* checkSecond = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkSecond, " 3 4 "
      " 1  0  1  0 "
      " 1  1  0  1 "
      " 0  0  1  1 "
    ) );

    ASSERT_TRUE( CMRchrmatCheckEqual(second, checkSecond) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkSecond) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkFirst) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 4 4 "
      " 1  1  0  0 "
      " 1  0  1  0 "
      " 0 -1  0  1 "
      " 0  0  1  1 "
    ) );
    CMR_CHRMAT* transpose = NULL;
    ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 4, 4, &sepa) );
    sepa->type = CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS;

    sepa->rowsFlags[0] = CMR_SEPA_FIRST;
    sepa->rowsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[3] = CMR_SEPA_SECOND;

    sepa->columnsFlags[0] = CMR_SEPA_FIRST;
    sepa->columnsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[3] = CMR_SEPA_SECOND;

    char epsilon;
    ASSERT_CMR_CALL( CMRdeltasumDecomposeEpsilon(cmr, matrix, transpose, sepa, &epsilon) );
    ASSERT_EQ(epsilon, +1);

    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( CMRdeltasumDecomposeFirst(cmr, matrix, sepa, epsilon, &first, NULL, NULL, NULL, NULL, NULL,
      NULL) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 3 4 "
      " 1  1  0  0 "
      " 1  0  1  1 "
      " 0 -1  0  1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( CMRdeltasumDecomposeSecond(cmr, matrix, sepa, epsilon, &second, NULL, NULL, NULL, NULL, NULL,
      NULL) );

    CMR_CHRMAT* checkSecond = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkSecond, " 3 4 "
      " 1  0  1  0 "
      " 1  1  0  1 "
      " 0  0  1  1 "
    ) );

    ASSERT_TRUE( CMRchrmatCheckEqual(second, checkSecond) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkSecond) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkFirst) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 4 4 "
      " 1  1  0  0 "
      " 1  0 -1  0 "
      " 0 -1  0  1 "
      " 0  0  1  1 "
    ) );
    CMR_CHRMAT* transpose = NULL;
    ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 4, 4, &sepa) );
    sepa->type = CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS;

    sepa->rowsFlags[0] = CMR_SEPA_FIRST;
    sepa->rowsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[3] = CMR_SEPA_SECOND;

    sepa->columnsFlags[0] = CMR_SEPA_FIRST;
    sepa->columnsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[3] = CMR_SEPA_SECOND;

    char epsilon;
    ASSERT_CMR_CALL( CMRdeltasumDecomposeEpsilon(cmr, matrix, transpose, sepa, &epsilon) );
    ASSERT_EQ(epsilon, -1);

    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( CMRdeltasumDecomposeFirst(cmr, matrix, sepa, epsilon, &first, NULL, NULL, NULL, NULL, NULL,
      NULL) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 3 4 "
      " 1  1  0  0 "
      " 1  0 -1 -1 "
      " 0 -1  0 -1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( CMRdeltasumDecomposeSecond(cmr, matrix, sepa, epsilon, &second, NULL, NULL, NULL, NULL, NULL,
      NULL) );

    CMR_CHRMAT* checkSecond = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkSecond, " 3 4 "
      "-1  0  1  0 "
      " 1  1  0  1 "
      " 0  0  1  1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(second, checkSecond) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkSecond) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkFirst) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Separation, ThreesumComposition)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &first, "4 5 "
      " 1  0  1  1 0 "
      " 0  1  1  1 0 "
      " 1  0  1  0 1 "
      " 0 -1  0 -1 1 "
    ) );
    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &second, "5 4 "
      " 1  1  0  0 "
      " 1  0  1  1 "
      " 0 -1  1  1 "
      " 1  0  1  0 "
      " 0 -1  0  1 "
    ) );

    CMR_CHRMAT* check = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &check, " 6 6 "
      " 1  0  1  1 0 0 "
      " 0  1  1  1 0 0 "
      " 1  0  1  0 1 1 "
      " 0 -1  0 -1 1 1 "
      " 1  0  1  0 1 0 "
      " 0 -1  0 -1 0 1 "
    ) );

    CMR_CHRMAT* threesum = NULL;
    size_t firstSpecialRows[2] = { 2, 3 };
    size_t firstSpecialColumns[3] = { 2, 3, 4 };
    size_t secondSpecialRows[3] = { 0, 1, 2 };
    size_t secondSpecialColumns[2] = { 0, 1 };
    ASSERT_CMR_CALL( CMRthreesumCompose(cmr, first, second, firstSpecialRows, firstSpecialColumns,
      secondSpecialRows, secondSpecialColumns, 3, &threesum) );

    CMRchrmatPrintDense(cmr, threesum, stdout, '0', false);

    ASSERT_TRUE( CMRchrmatCheckEqual(threesum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &threesum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
  }

  {
    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &first, "4 5 "
      " 1  0 -1  1 0 "
      " 0  1 -1  1 0 "
      " 1  0  1  0 1 "
      " 0 -1  1 -1 1 "
    ) );
    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &second, "5 4 "
      "-1  1  0  0 "
      " 1  0  1  1 "
      " 1 -1  1  1 "
      " 1  0  1  0 "
      " 0 -1  0  1 "
    ) );

    CMR_CHRMAT* check = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &check, " 6 6 "
      " 1  0 -1  1 0 0 "
      " 0  1 -1  1 0 0 "
      " 1  0  1  0 1 1 "
      " 0 -1  1 -1 1 1 "
      " 1  0  1  0 1 0 "
      "-1 -1  0 -1 0 1 "
    ) );

    CMR_CHRMAT* threesum = NULL;
    size_t firstSpecialRows[2] = { 2, 3 };
    size_t firstSpecialColumns[3] = { 2, 3, 4 };
    size_t secondSpecialRows[3] = { 0, 1, 2 };
    size_t secondSpecialColumns[2] = { 0, 1 };
    ASSERT_CMR_CALL( CMRthreesumCompose(cmr, first, second, firstSpecialRows, firstSpecialColumns,
      secondSpecialRows, secondSpecialColumns, 3, &threesum) );

    CMRchrmatPrintDense(cmr, threesum, stdout, '0', false);

    ASSERT_TRUE( CMRchrmatCheckEqual(threesum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &threesum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
  }

  {
    /* Inconsistent special matrix. */

    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &first, "4 5 "
      " 1  0  1  1 0 "
      " 0  1  1  1 0 "
      " 1  0  1  0 1 "
      " 0 -1  1 -1 1 "
    ) );
    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &second, "5 4 "
      " 1  1  0  0 "
      " 1  0  1  1 "
      " 0 -1  1  1 "
      " 1  0  1  0 "
      " 0 -1  0  1 "
    ) );

    CMR_CHRMAT* threesum = NULL;
    size_t firstSpecialRows[2] = { 2, 3 };
    size_t firstSpecialColumns[3] = { 2, 3, 4 };
    size_t secondSpecialRows[3] = { 0, 1, 2 };
    size_t secondSpecialColumns[2] = { 0, 1 };
    CMR_ERROR error = CMRthreesumCompose(cmr, first, second, firstSpecialRows, firstSpecialColumns,
      secondSpecialRows, secondSpecialColumns, 3, &threesum);

    ASSERT_EQ( error, CMR_ERROR_STRUCTURE );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


TEST(Separation, ThreesumDecomposition)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    /* Usual R_12 matrix. */

    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 6 6 "
      " 1  0  1  1 0 0 "
      " 0  1  1  1 0 0 "
      " 1  0  1  0 1 1 "
      " 0 -1  0 -1 1 1 "
      " 1  0  1  0 1 0 "
      " 0 -1  0 -1 0 1 "
    ) );
    CMR_CHRMAT* transpose = NULL;
    ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 6, 6, &sepa) );
    sepa->type = CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK;

    sepa->rowsFlags[0] = CMR_SEPA_FIRST;
    sepa->rowsFlags[1] = CMR_SEPA_FIRST;
    sepa->rowsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[3] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK2;
    sepa->rowsFlags[4] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[5] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK2;

    sepa->columnsFlags[0] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK2;
    sepa->columnsFlags[2] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[3] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK2;
    sepa->columnsFlags[4] = CMR_SEPA_SECOND;
    sepa->columnsFlags[5] = CMR_SEPA_SECOND;

    char gamma, beta;
    size_t specialRows[2];
    size_t specialColumns[2];
    ASSERT_CMR_CALL( CMRthreesumDecomposeSearchConnecting(cmr, matrix, transpose, sepa, specialRows, specialColumns,
      &gamma, &beta) );
    /* Rows are {2,4} and {3,5}. So in particular in [2,6), sum is odd. */
    ASSERT_GE( specialRows[0], 2 );
    ASSERT_LT( specialRows[0], 6 );
    ASSERT_GE( specialRows[1], 2 );
    ASSERT_LT( specialRows[1], 6 );
    ASSERT_EQ( (specialRows[0] + specialRows[1]) % 2, "one row odd, one row even" ? 1 : 0 );

    ASSERT_GE( specialColumns[0], 0 );
    ASSERT_LT( specialColumns[0], 4 );
    ASSERT_GE( specialColumns[1], 0 );
    ASSERT_LT( specialColumns[1], 4 );
    ASSERT_EQ( (specialColumns[0] + specialColumns[1]) % 2, "one column odd, one column even" ? 1 : 0 );

    /* Connecting matrix:
     *  1  1  0
     *  1  0  1
     *  0 -1  1
     */

    ASSERT_EQ( gamma, 1 );
    ASSERT_EQ( beta, 1 );

    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( CMRthreesumDecomposeFirst(cmr, matrix, sepa, specialRows, specialColumns, beta, &first, NULL,
      NULL, NULL, NULL, NULL, NULL) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 4 5 "
      " 1  0  1  1 0 "
      " 0  1  1  1 0 "
      " 1  0  1  0 1 "
      " 0 -1  0 -1 1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( CMRthreesumDecomposeSecond(cmr, matrix, sepa, specialRows, specialColumns, gamma, &second, NULL,
      NULL, NULL, NULL, NULL, NULL) );

    CMR_CHRMAT* checkSecond = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkSecond, " 5 4 "
      " 1  1 0 0 "
      " 1  0 1 1 "
      " 0 -1 1 1 "
      " 1  0 1 0 "
      " 0 -1 0 1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(second, checkSecond) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkSecond) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkFirst) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    /* R_12 matrix with flipped signs in (top-)left part. */

    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 6 6 "
      " -1  0  -1  1 0 0 "
      "  0 -1   1 -1 0 0 "
      " -1  0  -1  0 1 1 "
      "  0 -1   0 -1 1 1 "
      " -1  0  -1  0 1 0 "
      "  0 -1   0 -1 0 1 "
    ) );
    CMR_CHRMAT* transpose = NULL;
    ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 6, 6, &sepa) );
    sepa->type = CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK;

    sepa->rowsFlags[0] = CMR_SEPA_FIRST;
    sepa->rowsFlags[1] = CMR_SEPA_FIRST;
    sepa->rowsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[3] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK2;
    sepa->rowsFlags[4] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[5] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK2;

    sepa->columnsFlags[0] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK2;
    sepa->columnsFlags[2] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[3] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK2;
    sepa->columnsFlags[4] = CMR_SEPA_SECOND;
    sepa->columnsFlags[5] = CMR_SEPA_SECOND;

    char gamma, beta;
    size_t specialRows[2];
    size_t specialColumns[2];
    ASSERT_CMR_CALL( CMRthreesumDecomposeSearchConnecting(cmr, matrix, transpose, sepa, specialRows, specialColumns,
      &gamma, &beta) );
    /* Rows are {2,4} and {3,5}. So in particular in [2,6), sum is odd. */
    ASSERT_GE( specialRows[0], 2 );
    ASSERT_LT( specialRows[0], 6 );
    ASSERT_GE( specialRows[1], 2 );
    ASSERT_LT( specialRows[1], 6 );
    ASSERT_EQ( (specialRows[0] + specialRows[1]) % 2, "one row odd, one row even" ? 1 : 0 );

    ASSERT_GE( specialColumns[0], 0 );
    ASSERT_LT( specialColumns[0], 4 );
    ASSERT_GE( specialColumns[1], 0 );
    ASSERT_LT( specialColumns[1], 4 );
    ASSERT_EQ( (specialColumns[0] + specialColumns[1]) % 2, "one column odd, one column even" ? 1 : 0 );

    /* Connecting matrix:
     *  -1  1  0
     *  -1  0  1
     *   0 -1  1
     */

    ASSERT_EQ( gamma, -1 );
    ASSERT_EQ( beta, 1 );

    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( CMRthreesumDecomposeFirst(cmr, matrix, sepa, specialRows, specialColumns, beta, &first, NULL,
      NULL, NULL, NULL, NULL, NULL) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 4 5 "
      " -1  0  -1  1 0 "
      "  0 -1   1 -1 0 "
      " -1  0  -1  0 1 "
      "  0 -1   0 -1 1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( CMRthreesumDecomposeSecond(cmr, matrix, sepa, specialRows, specialColumns, gamma, &second, NULL,
      NULL, NULL, NULL, NULL, NULL) );

    CMR_CHRMAT* checkSecond = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkSecond, " 5 4 "
      " -1  1 0 0 "
      " -1  0 1 1 "
      "  0 -1 1 1 "
      " -1  0 1 0 "
      "  0 -1 0 1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(second, checkSecond) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkSecond) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkFirst) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    /* R_12 with some flipped signs in bottom(-right). */

    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 6 6 "
      " 1  0  1  1  0  0 "
      " 0  1  1  1  0  0 "
      " 1  0  1  0  1 -1 "
      " 0  1  0  1 -1  1 "
      " 1  0  1  0  1  0 "
      " 0  1  0  1 0   1 "
    ) );
    CMR_CHRMAT* transpose = NULL;
    ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 6, 6, &sepa) );
    sepa->type = CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK;

    sepa->rowsFlags[0] = CMR_SEPA_FIRST;
    sepa->rowsFlags[1] = CMR_SEPA_FIRST;
    sepa->rowsFlags[2] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[3] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK2;
    sepa->rowsFlags[4] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[5] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK2;

    sepa->columnsFlags[0] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK2;
    sepa->columnsFlags[2] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[3] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK2;
    sepa->columnsFlags[4] = CMR_SEPA_SECOND;
    sepa->columnsFlags[5] = CMR_SEPA_SECOND;

    char gamma, beta;
    size_t specialRows[2];
    size_t specialColumns[2];
    ASSERT_CMR_CALL( CMRthreesumDecomposeSearchConnecting(cmr, matrix, transpose, sepa, specialRows, specialColumns,
      &gamma, &beta) );
    /* Rows are {2,4} and {3,5}. So in particular in [2,6), sum is odd. */
    ASSERT_GE( specialRows[0], 2 );
    ASSERT_LT( specialRows[0], 6 );
    ASSERT_GE( specialRows[1], 2 );
    ASSERT_LT( specialRows[1], 6 );
    ASSERT_EQ( (specialRows[0] + specialRows[1]) % 2, "one row odd, one row even" ? 1 : 0 );

    ASSERT_GE( specialColumns[0], 0 );
    ASSERT_LT( specialColumns[0], 4 );
    ASSERT_GE( specialColumns[1], 0 );
    ASSERT_LT( specialColumns[1], 4 );
    ASSERT_EQ( (specialColumns[0] + specialColumns[1]) % 2, "one column odd, one column even" ? 1 : 0 );

    /* Connecting matrix:
     *  1  1  0
     *  1  0  1
     *  0  1 -1
     */

    ASSERT_EQ( gamma, 1 );
    ASSERT_EQ( beta, -1 );


    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( CMRthreesumDecomposeFirst(cmr, matrix, sepa, specialRows, specialColumns, beta, &first, NULL,
      NULL, NULL, NULL, NULL, NULL) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 4 5 "
      " 1  0  1  1  0 "
      " 0  1  1  1  0 "
      " 1  0  1  0  1 "
      " 0  1  0  1 -1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( CMRthreesumDecomposeSecond(cmr, matrix, sepa, specialRows, specialColumns, gamma, &second, NULL,
      NULL, NULL, NULL, NULL, NULL) );

    CMR_CHRMAT* checkSecond = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkSecond, " 5 4 "
      " 1  1 0  0 "
      " 1  0  1 -1 "
      " 0  1 -1  1 "
      " 1  0  1  0 "
      " 0  1  0  1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(second, checkSecond) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkSecond) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkFirst) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    /* A matrix with a non-identity special matrix (second child is not TU). */

    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 6 6 "
      " 0  1 -1  0 0 0 "
      " 0  0  1  1 0 0 "
      " 1  0  0  1 0 0 "
      " 0  0  0  0 1 1 "
      " 1 -1  0  0 1 0 "
      " 0 -1  0  0 0 1 "
    ) );
    CMR_CHRMAT* transpose = NULL;
    ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 6, 6, &sepa) );
    sepa->type = CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK;

    sepa->rowsFlags[0] = CMR_SEPA_FIRST;
    sepa->rowsFlags[1] = CMR_SEPA_FIRST;
    sepa->rowsFlags[2] = CMR_SEPA_FIRST;
    sepa->rowsFlags[3] = CMR_SEPA_SECOND;
    sepa->rowsFlags[4] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[5] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK2;

    sepa->columnsFlags[0] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK2;
    sepa->columnsFlags[2] = CMR_SEPA_FIRST;
    sepa->columnsFlags[3] = CMR_SEPA_FIRST;
    sepa->columnsFlags[4] = CMR_SEPA_SECOND;
    sepa->columnsFlags[5] = CMR_SEPA_SECOND;

    char gamma, beta;
    size_t specialRows[2];
    size_t specialColumns[2];
    ASSERT_CMR_CALL( CMRthreesumDecomposeSearchConnecting(cmr, matrix, transpose, sepa, specialRows, specialColumns,
      &gamma, &beta) );
    /* Rows are {2,4} and {3,5}. So in particular in [2,6), sum is odd. */
    ASSERT_GE( specialRows[0], 4 );
    ASSERT_LT( specialRows[0], 6 );
    ASSERT_GE( specialRows[1], 4 );
    ASSERT_LT( specialRows[1], 6 );
    ASSERT_EQ( (specialRows[0] + specialRows[1]) % 2, "one row odd, one row even" ? 1 : 0 );

    ASSERT_GE( specialColumns[0], 0 );
    ASSERT_LT( specialColumns[0], 2 );
    ASSERT_GE( specialColumns[1], 0 );
    ASSERT_LT( specialColumns[1], 2 );
    ASSERT_EQ( (specialColumns[0] + specialColumns[1]) % 2, "one column odd, one column even" ? 1 : 0 );

    /* Connecting matrix:
     * -1  1  0
     *  1 -1  1
     *  0 -1  1
     *
     * full det = 1 - 1 - 1   = -1 * 0 + ... (independent of gamma)
     */

    ASSERT_EQ( gamma, -1 );
    ASSERT_EQ( beta, 1 );

    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( CMRthreesumDecomposeFirst(cmr, matrix, sepa, specialRows, specialColumns, beta, &first, NULL,
      NULL, NULL, NULL, NULL, NULL) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 5 5 "
      " 0  1 -1  0 0 "
      " 0  0  1  1 0 "
      " 1  0  0  1 0 "
      " 1 -1  0  0 1 "
      " 0 -1  0  0 1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( CMRthreesumDecomposeSecond(cmr, matrix, sepa, specialRows, specialColumns, gamma, &second, NULL,
      NULL, NULL, NULL, NULL, NULL) );

    CMRchrmatPrintDense(cmr, second, stdout, '0', true);

    CMR_CHRMAT* checkSecond = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkSecond, " 4 4 "
      "-1  1  0 0 "
      " 0  0  1 1 "
      " 1 -1  1 0 "
      " 0 -1  0 1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(second, checkSecond) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkSecond) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkFirst) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    /* A matrix with a non-identity special matrix (second child is not TU). */

    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 6 6 "
      " 0  1 -1  0 0 0 "
      " 0  0  1  1 0 0 "
      " 1  0  0  1 0 0 "
      " 0  0  0  0 1 1 "
      " 1 -1  0  0 1 0 "
      " 0 -1  0  0 0 1 "
    ) );
    CMR_CHRMAT* transpose = NULL;
    ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 6, 6, &sepa) );
    sepa->type = CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK;

    sepa->rowsFlags[0] = CMR_SEPA_FIRST;
    sepa->rowsFlags[1] = CMR_SEPA_FIRST;
    sepa->rowsFlags[2] = CMR_SEPA_FIRST;
    sepa->rowsFlags[3] = CMR_SEPA_SECOND;
    sepa->rowsFlags[4] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[5] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK2;

    sepa->columnsFlags[0] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK2;
    sepa->columnsFlags[2] = CMR_SEPA_FIRST;
    sepa->columnsFlags[3] = CMR_SEPA_FIRST;
    sepa->columnsFlags[4] = CMR_SEPA_SECOND;
    sepa->columnsFlags[5] = CMR_SEPA_SECOND;

    char gamma, beta;
    size_t specialRows[2] = { 4, 5 };
    size_t specialColumns[2] = { 0, 1 };
    ASSERT_CMR_CALL( CMRthreesumDecomposeSignConnecting(cmr, matrix, transpose, sepa, specialRows, specialColumns,
      &gamma, &beta) );
    /* Rows are {2,4} and {3,5}. So in particular in [2,6), sum is odd. */
    ASSERT_EQ( specialRows[0], 4 );
    ASSERT_EQ( specialRows[1], 5 );
    ASSERT_EQ( specialColumns[0], 0 );
    ASSERT_EQ( specialColumns[1], 1 );

    /* Connecting matrix:
     * -1  1  0
     *  1 -1  1
     *  0 -1  1
     *
     * full det = 1 - 1 - 1   = -1 * 0 + ... (independent of gamma)
     */

    ASSERT_EQ( gamma, -1 );
    ASSERT_EQ( beta, 1 );

    CMR_CHRMAT* first = NULL;
    ASSERT_CMR_CALL( CMRthreesumDecomposeFirst(cmr, matrix, sepa, specialRows, specialColumns, beta, &first, NULL,
      NULL, NULL, NULL, NULL, NULL) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 5 5 "
      " 0  1 -1  0 0 "
      " 0  0  1  1 0 "
      " 1  0  0  1 0 "
      " 1 -1  0  0 1 "
      " 0 -1  0  0 1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    ASSERT_CMR_CALL( CMRthreesumDecomposeSecond(cmr, matrix, sepa, specialRows, specialColumns, gamma, &second, NULL,
      NULL, NULL, NULL, NULL, NULL) );

    CMRchrmatPrintDense(cmr, second, stdout, '0', true);

    CMR_CHRMAT* checkSecond = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkSecond, " 4 4 "
      "-1  1  0 0 "
      " 0  0  1 1 "
      " 1 -1  1 0 "
      " 0 -1  0 1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(second, checkSecond) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkSecond) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &checkFirst) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  {
    /* A matrix with an identity special matrix but non-working given special columns (second child is not TU). */

    CMR_CHRMAT* matrix = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, " 6 6 "
      " 0  1 -1  0 0 0 "
      " 0  0  1  1 0 0 "
      " 1  0  0  1 0 0 "
      " 0  0  0  0 1 1 "
      " 1  0  1  0 1 0 "
      " 0 -1  0  0 0 1 "
    ) );
    CMR_CHRMAT* transpose = NULL;
    ASSERT_CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );

    CMR_SEPA* sepa = NULL;
    ASSERT_CMR_CALL( CMRsepaCreate(cmr, 6, 6, &sepa) );
    sepa->type = CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK;

    sepa->rowsFlags[0] = CMR_SEPA_FIRST;
    sepa->rowsFlags[1] = CMR_SEPA_FIRST;
    sepa->rowsFlags[2] = CMR_SEPA_FIRST;
    sepa->rowsFlags[3] = CMR_SEPA_SECOND;
    sepa->rowsFlags[4] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK1;
    sepa->rowsFlags[5] = CMR_SEPA_SECOND | CMR_SEPA_FLAG_RANK2;

    sepa->columnsFlags[0] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[1] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK2;
    sepa->columnsFlags[2] = CMR_SEPA_FIRST | CMR_SEPA_FLAG_RANK1;
    sepa->columnsFlags[3] = CMR_SEPA_FIRST;
    sepa->columnsFlags[4] = CMR_SEPA_SECOND;
    sepa->columnsFlags[5] = CMR_SEPA_SECOND;

    char gamma, beta;
    size_t specialRows[2] = { 4, 5 };
    size_t specialColumns[2] = { 0, 1 };
    CMR_ERROR error = CMRthreesumDecomposeSignConnecting(cmr, matrix, transpose, sepa, specialRows, specialColumns,
      &gamma, &beta);
    ASSERT_EQ( error, CMR_ERROR_INPUT );

    ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &transpose) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
