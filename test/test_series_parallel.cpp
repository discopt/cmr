#include <gtest/gtest.h>

#include "common.h"
#include <cmr/series_parallel.h>

TEST(SeriesParallel, Empty)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_SP_REDUCTION reductions[2];
  size_t numReductions;

  {
    CMR_CHRMAT* mat0x0 = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &mat0x0, "0 0 "
    ) );

    ASSERT_CMR_CALL( CMRtestBinarySeriesParallel(cmr, mat0x0, true, NULL, reductions, &numReductions, NULL, NULL,
      NULL) );
    ASSERT_EQ( numReductions, 0 );

    ASSERT_CMR_CALL( CMRtestTernarySeriesParallel(cmr, mat0x0, true, NULL, reductions, &numReductions, NULL, NULL,
      NULL) );
    ASSERT_EQ( numReductions, 0 );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &mat0x0) );
  }

  {
    CMR_CHRMAT* mat2x0 = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &mat2x0, "2 0 "
    ) );

    ASSERT_CMR_CALL( CMRtestBinarySeriesParallel(cmr, mat2x0, true, NULL, reductions, &numReductions, NULL, NULL,
      NULL) );
    ASSERT_EQ( numReductions, 2 );
    ASSERT_EQ( reductions[0].element, -1 ); ASSERT_EQ( reductions[0].mate, 0 );
    ASSERT_EQ( reductions[1].element, -2 ); ASSERT_EQ( reductions[1].mate, 0 );

    ASSERT_CMR_CALL( CMRtestTernarySeriesParallel(cmr, mat2x0, true, NULL, reductions, &numReductions, NULL, NULL,
      NULL) );
    ASSERT_EQ( numReductions, 2 );
    ASSERT_EQ( reductions[0].element, -1 ); ASSERT_EQ( reductions[0].mate, 0 );
    ASSERT_EQ( reductions[1].element, -2 ); ASSERT_EQ( reductions[1].mate, 0 );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &mat2x0) );
  }

  {
    CMR_CHRMAT* mat0x2 = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &mat0x2, "0 2 "
    ) );

    ASSERT_CMR_CALL( CMRtestBinarySeriesParallel(cmr, mat0x2, true, NULL, reductions, &numReductions, NULL, NULL,
      NULL) );
    ASSERT_EQ( numReductions, 2 );
    ASSERT_EQ( reductions[0].element, 1 ); ASSERT_EQ( reductions[0].mate, 0 );
    ASSERT_EQ( reductions[1].element, 2 ); ASSERT_EQ( reductions[1].mate, 0 );

    ASSERT_CMR_CALL( CMRtestTernarySeriesParallel(cmr, mat0x2, true, NULL, reductions, &numReductions, NULL, NULL,
      NULL) );
    ASSERT_EQ( numReductions, 2 );
    ASSERT_EQ( reductions[0].element, 1 ); ASSERT_EQ( reductions[0].mate, 0 );
    ASSERT_EQ( reductions[1].element, 2 ); ASSERT_EQ( reductions[1].mate, 0 );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &mat0x2) );
  }

  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(SeriesParallel, BinaryReduction)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "20 20 "
    "0 1 0 1 1 1 1 0 1 1 1 0 1 0 0 0 1 1 1 1 "
    "1 1 0 0 0 1 0 0 0 1 1 0 1 0 0 0 1 1 0 0 "
    "1 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 1 "
    "1 1 0 0 0 1 0 1 0 1 1 0 1 0 0 0 1 0 0 0 "
    "1 1 0 1 1 1 0 0 1 1 1 1 1 0 0 1 1 1 0 1 "
    "1 1 0 0 0 1 0 0 0 1 1 0 1 0 0 0 1 1 0 0 "
    "0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 "
    "0 1 1 1 1 1 1 0 1 1 1 0 1 0 0 0 1 1 1 1 "
    "1 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 1 "
    "0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 "
    "0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 0 0 1 1 1 "
    "1 0 0 1 1 0 0 0 1 0 0 0 0 0 0 0 0 1 1 1 "
    "0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 1 1 0 "
    "0 1 0 0 0 1 1 1 0 1 1 0 1 0 0 0 1 1 0 1 "
    "1 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0 1 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "1 1 0 0 0 1 0 0 0 1 1 0 1 0 0 0 1 1 0 0 "
    "1 1 0 0 0 1 0 0 0 1 1 0 1 0 0 0 1 1 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 "
  ) );

  CMR_SP_REDUCTION operations[40];
  size_t numOperations;
  CMR_SUBMAT* submatrix = NULL;

  ASSERT_CMR_CALL( CMRtestBinarySeriesParallel(cmr, matrix, true, NULL, operations, &numOperations, &submatrix, NULL,
    NULL) );
  ASSERT_EQ( numOperations, 20);
  ASSERT_EQ( operations[0].element, -6);  ASSERT_EQ( operations[0].mate, -2);
  ASSERT_EQ( operations[1].element, -7);  ASSERT_EQ( operations[1].mate, 5);
  ASSERT_EQ( operations[2].element, -8);  ASSERT_EQ( operations[2].mate, 14);
  ASSERT_EQ( operations[3].element, -10); ASSERT_EQ( operations[3].mate, -3);
  ASSERT_EQ( operations[4].element, -16); ASSERT_EQ( operations[4].mate, -3);
  ASSERT_EQ( operations[5].element, -17); ASSERT_EQ( operations[5].mate, 0);
  ASSERT_EQ( operations[6].element, -18); ASSERT_EQ( operations[6].mate, -2);
  ASSERT_EQ( operations[7].element, -19); ASSERT_EQ( operations[7].mate, -2);
  ASSERT_EQ( operations[8].element, -20); ASSERT_EQ( operations[8].mate, 19);
  ASSERT_EQ( operations[9].element, 3);   ASSERT_EQ( operations[9].mate, -9);
  ASSERT_EQ( operations[10].element, 10); ASSERT_EQ( operations[10].mate, 6);
  ASSERT_EQ( operations[11].element, 11); ASSERT_EQ( operations[11].mate, 6);
  ASSERT_EQ( operations[12].element, 12); ASSERT_EQ( operations[12].mate, -5);
  ASSERT_EQ( operations[13].element, 13); ASSERT_EQ( operations[13].mate, 6);
  ASSERT_EQ( operations[14].element, 14); ASSERT_EQ( operations[14].mate, 0);
  ASSERT_EQ( operations[15].element, 16); ASSERT_EQ( operations[15].mate, -5);
  ASSERT_EQ( operations[16].element, 17); ASSERT_EQ( operations[16].mate, 6);
  ASSERT_EQ( operations[17].element, 2);  ASSERT_EQ( operations[17].mate, 6);
  ASSERT_EQ( operations[18].element, 5);  ASSERT_EQ( operations[18].mate, 4);
  ASSERT_EQ( operations[19].element, -1); ASSERT_EQ( operations[19].mate, -9);

  ASSERT_EQ( submatrix->numRows, 10 );
  ASSERT_EQ( submatrix->numColumns, 10 );

  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(SeriesParallel, BinaryShortWheel)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "20 20 "
    "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "1 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 0 0 "
    "0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 1 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 1 0 0 "
    "0 0 0 0 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 "
    "0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 1 1 1 0 1 0 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 "
  ) );

  CMR_SP_REDUCTION operations[40];
  size_t numOperations;
  CMR_SUBMAT* wheelSubmatrix = NULL;

  ASSERT_CMR_CALL( CMRtestBinarySeriesParallel(cmr, matrix, true, NULL, operations, &numOperations, NULL,
    &wheelSubmatrix, NULL) );
  ASSERT_EQ( numOperations, 8 );

  CMR_CHRMAT* wheelMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatFilterSubmat(cmr, matrix, wheelSubmatrix, &wheelMatrix) );

  CMRchrmatPrintDense(cmr, stdout, wheelMatrix, '0', true);
  
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &wheelMatrix) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &wheelSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(SeriesParallel, BinaryLongWheel)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "20 20 "
    "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "1 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 "
    "0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 "
    "0 0 0 0 1 1 1 1 0 0 0 0 1 0 0 0 0 1 0 0 "
    "0 0 0 0 1 1 1 1 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 "
    "0 0 0 0 1 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 1 1 1 0 0 0 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 "
  ) );

  CMR_SP_REDUCTION operations[40];
  size_t numOperations;
  CMR_SUBMAT* wheelSubmatrix = NULL;

  ASSERT_CMR_CALL( CMRtestBinarySeriesParallel(cmr, matrix, true, NULL, operations, &numOperations, NULL,
    &wheelSubmatrix, NULL) );
  ASSERT_EQ( numOperations, 8 );
  for (size_t o = 0; o < numOperations; ++o)
  {
    printf("%s\n", CMRspReductionString(operations[o], NULL));
  }
  
  CMR_CHRMAT* wheelMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatFilterSubmat(cmr, matrix, wheelSubmatrix, &wheelMatrix) );

  CMRchrmatPrintDense(cmr, stdout, wheelMatrix, '0', true);
  
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &wheelMatrix) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &wheelSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(SeriesParallel, BinarySpecialWheel)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "20 20 "
    "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "1 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 "
    "0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 "
    "0 0 0 0 1 1 1 1 0 0 0 0 1 0 0 0 0 1 0 0 "
    "0 0 0 0 1 1 1 1 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 "
    "0 0 0 0 1 1 1 0 1 0 0 0 1 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 1 1 1 0 0 0 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 "
  ) );

  CMR_SP_REDUCTION operations[40];
  size_t numOperations;
  CMR_SUBMAT* wheelSubmatrix = NULL;

  ASSERT_CMR_CALL( CMRtestBinarySeriesParallel(cmr, matrix, true, NULL, operations, &numOperations, NULL,
    &wheelSubmatrix, NULL) );
  ASSERT_EQ( numOperations, 8 );
  for (size_t o = 0; o < numOperations; ++o)
  {
    printf("%s\n", CMRspReductionString(operations[o], NULL));
  }

  CMR_CHRMAT* wheelMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatFilterSubmat(cmr, matrix, wheelSubmatrix, &wheelMatrix) );

  CMRchrmatPrintDense(cmr, stdout, wheelMatrix, '0', true);
  
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &wheelMatrix) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &wheelSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(SeriesParallel, BinarySeparation)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "20 20 "
    "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "1 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 "
    "0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 "
    "0 0 0 0 1 1 1 1 0 0 0 0 1 0 0 0 0 1 0 0 "
    "0 0 0 0 1 1 1 1 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 "
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 "
    "0 0 0 0 1 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 1 1 1 0 0 0 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 "
  ) );

  CMR_SP_REDUCTION operations[40];
  size_t numOperations;
  CMR_SUBMAT* wheelSubmatrix = NULL;
  CMR_ELEMENT separationRank1Elements[40];
  size_t numSeparationRank1Elements;

  ASSERT_CMR_CALL( CMRdecomposeBinarySeriesParallel(cmr, matrix, true, NULL, operations, &numOperations, NULL, &wheelSubmatrix,
    separationRank1Elements, &numSeparationRank1Elements, NULL) );
  ASSERT_EQ( numOperations, 8 );
  ASSERT_EQ( numSeparationRank1Elements, 17 );
  ASSERT_EQ( separationRank1Elements[0], CMRrowToElement(4) );
  ASSERT_EQ( separationRank1Elements[1], CMRrowToElement(5) );
  ASSERT_EQ( separationRank1Elements[2], CMRrowToElement(6) );
  ASSERT_EQ( separationRank1Elements[3], CMRrowToElement(7) );
  ASSERT_EQ( separationRank1Elements[4], CMRrowToElement(8) );
  ASSERT_EQ( separationRank1Elements[5], CMRrowToElement(9) );
  ASSERT_EQ( separationRank1Elements[6], CMRrowToElement(10) );
  ASSERT_EQ( separationRank1Elements[7], CMRrowToElement(11) );
  ASSERT_EQ( separationRank1Elements[8], CMRrowToElement(12) );
  ASSERT_EQ( separationRank1Elements[9], CMRcolumnToElement(4) );
  ASSERT_EQ( separationRank1Elements[10], CMRcolumnToElement(5) );
  ASSERT_EQ( separationRank1Elements[11], CMRcolumnToElement(6) );
  ASSERT_EQ( separationRank1Elements[12], CMRcolumnToElement(7) );
  ASSERT_EQ( separationRank1Elements[13], CMRcolumnToElement(8) );
  ASSERT_EQ( separationRank1Elements[14], CMRcolumnToElement(9) );
  ASSERT_EQ( separationRank1Elements[15], CMRcolumnToElement(10) );
  ASSERT_EQ( separationRank1Elements[16], CMRcolumnToElement(11) );
  
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &wheelSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(SeriesParallel, BinaryWheelAfterSeparation)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "20 20 "
    "1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 " // 1
    "1 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 " // 2
    "0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 " // 3
    "0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 " // 4
    "0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 " // 5
    "0 0 0 0 1 1 1 1 0 0 0 0 1 0 0 0 0 1 0 0 " // 6
    "0 0 0 0 1 1 1 1 0 0 0 0 1 1 0 0 0 0 0 0 " // 7
    "0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 " // 8
    "0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 " // 9
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 " // 10
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 " // 11
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 " // 12
    "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 " // 13
    "0 0 0 0 1 1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 1 1 1 0 0 0 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 "
  ) );

  CMR_SP_REDUCTION reductions[40];
  size_t numReductions;
  CMR_SUBMAT* wheelSubmatrix = NULL;

  ASSERT_CMR_CALL( CMRtestBinarySeriesParallel(cmr, matrix, true, NULL, reductions, &numReductions, NULL,
    &wheelSubmatrix, NULL) );
  ASSERT_EQ( numReductions, 8 );
  for (size_t o = 0; o < numReductions; ++o)
  {
    printf("%s\n", CMRspReductionString(reductions[o], NULL));
  }

  CMR_CHRMAT* wheelMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatFilterSubmat(cmr, matrix, wheelSubmatrix, &wheelMatrix) );

  CMRchrmatPrintDense(cmr, stdout, wheelMatrix, '0', true);
  
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &wheelMatrix) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &wheelSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


TEST(SeriesParallel, TernarySeriesParallel)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "10 10 "
    " 0 -1  0  0  1 1  0  1  0  1 "
    " 1  0  0  1  0 0  0  0  0  0 "
    "-1  0  0 -1  0 0  0  0  0  0 "
    "-1  1  0 -1 -1 0  1 -1 -1 -1 "
    " 0  0 -1  0  0 0  0  0  0  0 "
    " 0  0  0  0  0 0  0 -1  0  0 "
    " 0 -1  0  0  0 0  0  0  0  1 "
    " 0  1  0  0  0 0  0  0  0 -1 "
    " 1 -1  0  1  1 0 -1  1  1  1 "
    " 0  0 -1  0  0 0  0  0  0  0 "
  ) );

  CMR_SP_REDUCTION reductions[20];
  size_t numReductions;

  ASSERT_CMR_CALL( CMRtestTernarySeriesParallel(cmr, matrix, true, NULL, reductions, &numReductions, NULL, NULL, NULL) );
  ASSERT_EQ( numReductions, 20 );
  for (size_t o = 0; o < numReductions; ++o)
  {
    printf("%s\n", CMRspReductionString(reductions[o], NULL));
  }

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


TEST(SeriesParallel, TernaryNonbinary)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "10 10 "
    " 1 1  1 1  1 0  1 0 1 1 "
    " 0 1  0 0  0 0  0 0 0 0 " 
    " 1 1  1 1  1 0  1 1 0 1 "
    " 1 1 -1 1  1 0  1 1 0 1 "
    " 1 1  1 1 -1 0 -1 1 0 1 "
    " 1 0  1 0  0 1  0 0 0 0 "
    " 1 1  1 1  1 0  1 0 0 1 "
    " 0 1  0 0  0 0  0 0 0 0 "
    " 1 1 -1 1  1 0  1 1 0 1 "
    " 1 1  1 1  1 0  1 1 0 1 "
  ) );

  CMR_SP_REDUCTION reductions[10];
  size_t numReductions;
  CMR_SUBMAT* violatorSubmatrix = NULL;

  ASSERT_CMR_CALL( CMRtestTernarySeriesParallel(cmr, matrix, true, NULL, reductions, &numReductions, NULL,
    &violatorSubmatrix, NULL) );
  ASSERT_EQ( numReductions, 10 );
  for (size_t o = 0; o < numReductions; ++o)
  {
    printf("%s\n", CMRspReductionString(reductions[o], NULL));
  }

  CMR_CHRMAT* violatorMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatFilterSubmat(cmr, matrix, violatorSubmatrix, &violatorMatrix) );
  ASSERT_EQ( violatorMatrix->numNonzeros, 4 );
  int det = violatorMatrix->entryValues[0] * violatorMatrix->entryValues[3]
    - violatorMatrix->entryValues[1] * violatorMatrix->entryValues[2];
  ASSERT_EQ( abs(det), 2 );

  ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, stdout, violatorMatrix, '0', true) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &violatorMatrix) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(SeriesParallel, TernaryWheel)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "10 10 "
    "-1  0 -1  0  0  1  0  0  0  0 "
    " 0 -1  0  0  0  0  0  0  0  0 "
    "-1  0 -1  0 -1  1  0  0  0 -1 "
    "-1  0 -1  0 -1  1  0  0  0 -1 "
    " 0  1  0  0 -1  0  0  0  0  0 "
    " 0  1  0  0 -1  0  0  0  0  0 "
    " 1 -1  1 -1  0 -1 -1 +1  0  1 "
    " 0 -1  0  0  1  0  0  0 +1  0 "
    " 0  1  0  0 -1  0  0  0  0  0 "
    " 0 -1  0  0  1  0  0  0  0  0 "
  ) );

  CMR_SP_REDUCTION reductions[20];
  size_t numReductions;
  CMR_SUBMAT* violatorSubmatrix = NULL;

  ASSERT_CMR_CALL( CMRtestTernarySeriesParallel(cmr, matrix, true, NULL, reductions, &numReductions, NULL,
    &violatorSubmatrix, NULL) );
  ASSERT_EQ( numReductions, 14 );
  for (size_t o = 0; o < numReductions; ++o)
  {
    printf("%s\n", CMRspReductionString(reductions[o], NULL));
  }

  CMR_CHRMAT* violatorMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatFilterSubmat(cmr, matrix, violatorSubmatrix, &violatorMatrix) );
  ASSERT_EQ( violatorMatrix->numNonzeros, 6 );

  ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, stdout, violatorMatrix, '0', true) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &violatorMatrix) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
