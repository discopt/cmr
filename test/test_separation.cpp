#include <gtest/gtest.h>

#include "common.h"
#include <cmr/separation.h>

TEST(Separation, OneSum)
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
    ASSERT_CMR_CALL( CMRoneSumCompose(cmr, 3, matrices, &onesum) );

    ASSERT_TRUE( CMRchrmatCheckEqual(onesum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &onesum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &third) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );  
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Separation, TwoSumComposition)
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
    ASSERT_CMR_CALL( CMRtwoSumCompose(cmr, first, second, CMRrowToElement(1), CMRcolumnToElement(2), 3, &twosum) );

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
    ASSERT_CMR_CALL( CMRtwoSumCompose(cmr, first, second, CMRcolumnToElement(4), CMRrowToElement(0), 3, &twosum) );

    CMRchrmatPrintDense(cmr, twosum, stdout, '0', false);

    ASSERT_TRUE( CMRchrmatCheckEqual(twosum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &twosum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Separation, TwoSumDecomposition)
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
    CMR_ELEMENT firstMarker;
    ASSERT_CMR_CALL( CMRtwoSumDecomposeFirst(cmr, matrix, sepa, &first, NULL, NULL, NULL, NULL, &firstMarker) );

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
    CMR_ELEMENT secondMarker;
    ASSERT_CMR_CALL( CMRtwoSumDecomposeSecond(cmr, matrix, sepa, &second, NULL, NULL, NULL, NULL, &secondMarker) );
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
    ASSERT_CMR_CALL( CMRtwoSumCompose(cmr, first, second, firstMarker, secondMarker, 0, &check) );

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
    CMR_ELEMENT firstMarker;
    ASSERT_CMR_CALL( CMRtwoSumDecomposeFirst(cmr, matrix, sepa, &first, NULL, NULL, NULL, NULL, &firstMarker) );

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
    CMR_ELEMENT secondMarker;
    ASSERT_CMR_CALL( CMRtwoSumDecomposeSecond(cmr, matrix, sepa, &second, NULL, NULL, NULL, NULL, &secondMarker) );
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
    ASSERT_CMR_CALL( CMRtwoSumCompose(cmr, second, first, secondMarker, firstMarker, 0, &check) );

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

TEST(Separation, ThreeSumSeymourComposition)
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
    ASSERT_CMR_CALL( CMRthreeSumSeymourCompose(cmr, first, second, CMRrowToElement(2), CMRcolumnToElement(2),
      CMRcolumnToElement(3), CMRrowToElement(0), CMRcolumnToElement(0), CMRcolumnToElement(1), 3, &threesum) );

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
    ASSERT_CMR_CALL( CMRthreeSumSeymourCompose(cmr, first, second, CMRrowToElement(2), CMRcolumnToElement(2),
      CMRcolumnToElement(3), CMRrowToElement(0), CMRcolumnToElement(0), CMRcolumnToElement(1), 3, &threesum) );

    CMRchrmatPrintDense(cmr, threesum, stdout, '0', false);

    ASSERT_TRUE( CMRchrmatCheckEqual(threesum, check) );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &threesum) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &check) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &second) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &first) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Separation, ThreeSumSeymourDecomposition)
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
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeEpsilon(cmr, matrix, transpose, sepa, &epsilon) );
    ASSERT_EQ(epsilon, -1);

    CMR_CHRMAT* first = NULL;
    CMR_ELEMENT firstMarker1 = 0;
    CMR_ELEMENT firstMarker2 = 0;
    CMR_ELEMENT firstMarker3 = 0;
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeFirst(cmr, matrix, sepa, epsilon, &first, NULL, NULL, NULL, NULL,
      &firstMarker1, &firstMarker2, &firstMarker3) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 3 4 "
      " 1  1  0  0 "
      " 1  0  1  1 "
      " 0  1  0 -1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    CMR_ELEMENT secondMarker1 = 0;
    CMR_ELEMENT secondMarker2 = 0;
    CMR_ELEMENT secondMarker3 = 0;
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeSecond(cmr, matrix, sepa, epsilon, &second, NULL, NULL, NULL, NULL,
      &secondMarker1, &secondMarker2, &secondMarker3) );

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
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeEpsilon(cmr, matrix, transpose, sepa, &epsilon) );
    ASSERT_EQ(epsilon, +1);

    CMR_CHRMAT* first = NULL;
    CMR_ELEMENT firstMarker1 = 0;
    CMR_ELEMENT firstMarker2 = 0;
    CMR_ELEMENT firstMarker3 = 0;
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeFirst(cmr, matrix, sepa, epsilon, &first, NULL, NULL, NULL, NULL,
      &firstMarker1, &firstMarker2, &firstMarker3) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 3 4 "
      " 1  1  0  0 "
      " 1  0 -1 -1 "
      " 0  1  0  1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    CMR_ELEMENT secondMarker1 = 0;
    CMR_ELEMENT secondMarker2 = 0;
    CMR_ELEMENT secondMarker3 = 0;
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeSecond(cmr, matrix, sepa, epsilon, &second, NULL, NULL, NULL, NULL,
      &secondMarker1, &secondMarker2, &secondMarker3) );

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
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeEpsilon(cmr, matrix, transpose, sepa, &epsilon) );
    ASSERT_EQ(epsilon, +1);

    CMR_CHRMAT* first = NULL;
    CMR_ELEMENT firstMarker1 = 0;
    CMR_ELEMENT firstMarker2 = 0;
    CMR_ELEMENT firstMarker3 = 0;
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeFirst(cmr, matrix, sepa, epsilon, &first, NULL, NULL, NULL, NULL,
      &firstMarker1, &firstMarker2, &firstMarker3) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 3 4 "
      " 1  1  0  0 "
      " 1  0  1  1 "
      " 0 -1  0  1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    CMR_ELEMENT secondMarker1 = 0;
    CMR_ELEMENT secondMarker2 = 0;
    CMR_ELEMENT secondMarker3 = 0;
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeSecond(cmr, matrix, sepa, epsilon, &second, NULL, NULL, NULL, NULL,
      &secondMarker1, &secondMarker2, &secondMarker3) );

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
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeEpsilon(cmr, matrix, transpose, sepa, &epsilon) );
    ASSERT_EQ(epsilon, -1);

    CMR_CHRMAT* first = NULL;
    CMR_ELEMENT firstMarker1 = 0;
    CMR_ELEMENT firstMarker2 = 0;
    CMR_ELEMENT firstMarker3 = 0;
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeFirst(cmr, matrix, sepa, epsilon, &first, NULL, NULL, NULL, NULL,
      &firstMarker1, &firstMarker2, &firstMarker3) );

    CMR_CHRMAT* checkFirst = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &checkFirst, " 3 4 "
      " 1  1  0  0 "
      " 1  0 -1 -1 "
      " 0 -1  0 -1 "
    ) );
    ASSERT_TRUE( CMRchrmatCheckEqual(first, checkFirst) );

    CMR_CHRMAT* second = NULL;
    CMR_ELEMENT secondMarker1 = 0;
    CMR_ELEMENT secondMarker2 = 0;
    CMR_ELEMENT secondMarker3 = 0;
    ASSERT_CMR_CALL( CMRthreeSumSeymourDecomposeSecond(cmr, matrix, sepa, epsilon, &second, NULL, NULL, NULL, NULL,
      &secondMarker1, &secondMarker2, &secondMarker3) );

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
