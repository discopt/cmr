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

    ASSERT_CMR_CALL( CMRspTestBinary(cmr, mat0x0, NULL, reductions, &numReductions, NULL, NULL,
      NULL, DBL_MAX) );
    ASSERT_EQ( numReductions, 0UL );

    ASSERT_CMR_CALL( CMRspTestTernary(cmr, mat0x0, NULL, reductions, &numReductions, NULL, NULL,
      NULL, DBL_MAX) );
    ASSERT_EQ( numReductions, 0UL );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &mat0x0) );
  }

  {
    CMR_CHRMAT* mat2x0 = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &mat2x0, "2 0 "
    ) );

    ASSERT_CMR_CALL( CMRspTestBinary(cmr, mat2x0, NULL, reductions, &numReductions, NULL, NULL,
      NULL, DBL_MAX) );
    ASSERT_EQ( numReductions, 2UL );
    ASSERT_EQ( reductions[0].element, -1 ); ASSERT_EQ( reductions[0].mate, 0 );
    ASSERT_EQ( reductions[1].element, -2 ); ASSERT_EQ( reductions[1].mate, 0 );

    ASSERT_CMR_CALL( CMRspTestTernary(cmr, mat2x0, NULL, reductions, &numReductions, NULL, NULL,
      NULL, DBL_MAX) );
    ASSERT_EQ( numReductions, 2UL );
    ASSERT_EQ( reductions[0].element, -1 ); ASSERT_EQ( reductions[0].mate, 0 );
    ASSERT_EQ( reductions[1].element, -2 ); ASSERT_EQ( reductions[1].mate, 0 );

    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &mat2x0) );
  }

  {
    CMR_CHRMAT* mat0x2 = NULL;
    ASSERT_CMR_CALL( stringToCharMatrix(cmr, &mat0x2, "0 2 "
    ) );

    ASSERT_CMR_CALL( CMRspTestBinary(cmr, mat0x2, NULL, reductions, &numReductions, NULL, NULL,
      NULL, DBL_MAX) );
    ASSERT_EQ( numReductions, 2UL );
    ASSERT_EQ( reductions[0].element, 1 ); ASSERT_EQ( reductions[0].mate, 0 );
    ASSERT_EQ( reductions[1].element, 2 ); ASSERT_EQ( reductions[1].mate, 0 );

    ASSERT_CMR_CALL( CMRspTestTernary(cmr, mat0x2, NULL, reductions, &numReductions, NULL, NULL,
      NULL, DBL_MAX) );
    ASSERT_EQ( numReductions, 2UL );
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

  CMR_SP_REDUCTION reductions[40];
  size_t numReductions;
  CMR_SUBMAT* submatrix = NULL;

  ASSERT_CMR_CALL( CMRspTestBinary(cmr, matrix, NULL, reductions, &numReductions, &submatrix, NULL, NULL,
    DBL_MAX) );
  ASSERT_EQ( numReductions, 20UL);
  ASSERT_EQ( reductions[0].element, -6);  ASSERT_EQ( reductions[0].mate, -2);
  ASSERT_EQ( reductions[1].element, -7);  ASSERT_EQ( reductions[1].mate, 5);
  ASSERT_EQ( reductions[2].element, -8);  ASSERT_EQ( reductions[2].mate, 14);
  ASSERT_EQ( reductions[3].element, -10); ASSERT_EQ( reductions[3].mate, -3);
  ASSERT_EQ( reductions[4].element, -16); ASSERT_EQ( reductions[4].mate, -3);
  ASSERT_EQ( reductions[5].element, -17); ASSERT_EQ( reductions[5].mate, 0);
  ASSERT_EQ( reductions[6].element, -18); ASSERT_EQ( reductions[6].mate, -2);
  ASSERT_EQ( reductions[7].element, -19); ASSERT_EQ( reductions[7].mate, -2);
  ASSERT_EQ( reductions[8].element, -20); ASSERT_EQ( reductions[8].mate, 19);
  ASSERT_EQ( reductions[9].element, 3);   ASSERT_EQ( reductions[9].mate, -9);
  ASSERT_EQ( reductions[10].element, 10); ASSERT_EQ( reductions[10].mate, 6);
  ASSERT_EQ( reductions[11].element, 11); ASSERT_EQ( reductions[11].mate, 6);
  ASSERT_EQ( reductions[12].element, 12); ASSERT_EQ( reductions[12].mate, -5);
  ASSERT_EQ( reductions[13].element, 13); ASSERT_EQ( reductions[13].mate, 6);
  ASSERT_EQ( reductions[14].element, 14); ASSERT_EQ( reductions[14].mate, 0);
  ASSERT_EQ( reductions[15].element, 16); ASSERT_EQ( reductions[15].mate, -5);
  ASSERT_EQ( reductions[16].element, 17); ASSERT_EQ( reductions[16].mate, 6);
  ASSERT_EQ( reductions[17].element, 2);  ASSERT_EQ( reductions[17].mate, 6);
  ASSERT_EQ( reductions[18].element, 5);  ASSERT_EQ( reductions[18].mate, 4);
  ASSERT_EQ( reductions[19].element, -1); ASSERT_EQ( reductions[19].mate, -9);

  ASSERT_EQ( submatrix->numRows, 10UL );
  ASSERT_EQ( submatrix->numColumns, 10UL );

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

  ASSERT_CMR_CALL( CMRspTestBinary(cmr, matrix, NULL, operations, &numOperations, NULL,
    &wheelSubmatrix, NULL, DBL_MAX) );
  ASSERT_EQ( numOperations, 8UL );

  CMR_CHRMAT* wheelMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatSlice(cmr, matrix, wheelSubmatrix, &wheelMatrix) );

  CMRchrmatPrintDense(cmr, wheelMatrix, stdout, '0', true);
  
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

  CMR_SP_REDUCTION reductions[40];
  size_t numReductions;
  CMR_SUBMAT* wheelSubmatrix = NULL;

  ASSERT_CMR_CALL( CMRspTestBinary(cmr, matrix, NULL, reductions, &numReductions, NULL, &wheelSubmatrix,
    NULL, DBL_MAX) );
  ASSERT_EQ( numReductions, 8UL );
  for (size_t o = 0; o < numReductions; ++o)
  {
    printf("%s\n", CMRspReductionString(reductions[o], NULL));
  }
  
  CMR_CHRMAT* wheelMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatSlice(cmr, matrix, wheelSubmatrix, &wheelMatrix) );

  CMRchrmatPrintDense(cmr, wheelMatrix, stdout, '0', true);
  
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

  CMR_SP_REDUCTION reductions[40];
  size_t numReductions;
  CMR_SUBMAT* wheelSubmatrix = NULL;

  ASSERT_CMR_CALL( CMRspTestBinary(cmr, matrix, NULL, reductions, &numReductions, NULL,
    &wheelSubmatrix, NULL, DBL_MAX) );
  ASSERT_EQ( numReductions, 8UL );
  for (size_t o = 0; o < numReductions; ++o)
  {
    printf("%s\n", CMRspReductionString(reductions[o], NULL));
  }

  CMR_CHRMAT* wheelMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatSlice(cmr, matrix, wheelSubmatrix, &wheelMatrix) );

  CMRchrmatPrintDense(cmr, wheelMatrix, stdout, '0', true);
  
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &wheelMatrix) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &wheelSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(SeriesParallel, BinarySeparationFirstSearch)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "11 11 "
    " 1 1 1 0 0 0 0 0 0 0 0 "
    " 0 0 1 1 0 0 0 0 0 0 0 "
    " 0 1 0 1 0 0 0 0 0 0 0 "
    " 1 0 0 0 1 1 0 0 0 0 0 "
    " 1 0 0 0 0 0 0 0 0 0 0 "
    " 0 0 0 0 1 1 0 0 0 0 1 "
    " 0 0 0 0 0 1 1 0 0 0 0 "
    " 0 0 0 0 1 0 1 0 1 0 0 "
    " 0 0 0 0 0 0 1 1 0 0 0 "
    " 0 0 0 0 0 0 0 0 1 1 0 "
    " 0 0 0 0 0 0 0 0 0 1 0 "
  ) );

  CMR_SP_REDUCTION reductions[22];
  size_t numReductions;
  CMR_SUBMAT* wheelSubmatrix = NULL;
  CMR_SEPA* sepa = NULL;


  CMR_SUBMAT* reducedSubmatrix = NULL;
  CMR_CHRMAT* reducedMatrix = NULL;

  ASSERT_CMR_CALL( CMRspDecomposeBinary(cmr, matrix, NULL, reductions, SIZE_MAX, &numReductions,
    &reducedSubmatrix, &wheelSubmatrix, &sepa, NULL, DBL_MAX) );
  ASSERT_EQ( numReductions, 8UL );


  ASSERT_CMR_CALL( CMRchrmatSlice(cmr, matrix, reducedSubmatrix, &reducedMatrix) );

  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);

  CMRchrmatPrintDense(cmr, reducedMatrix, stdout, '0', true);

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &reducedMatrix) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &reducedSubmatrix) );

  // TODO: CHECK flags.

  ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &wheelSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(SeriesParallel, BinarySeparationSecondSearch)
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

  CMR_SP_REDUCTION reductions[40];
  size_t numReductions;
  CMR_SUBMAT* wheelSubmatrix = NULL;
  CMR_SEPA* sepa = NULL;

  ASSERT_CMR_CALL( CMRspDecomposeBinary(cmr, matrix, NULL, reductions, SIZE_MAX, &numReductions, NULL,
    &wheelSubmatrix, &sepa, NULL, DBL_MAX) );
  ASSERT_EQ( numReductions, 8UL );

  // TODO: CHECK flags.


  ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
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

  ASSERT_CMR_CALL( CMRspTestBinary(cmr, matrix, NULL, reductions, &numReductions, NULL,
    &wheelSubmatrix, NULL, DBL_MAX) );
  ASSERT_EQ( numReductions, 8UL );
  for (size_t o = 0; o < numReductions; ++o)
  {
    printf("%s\n", CMRspReductionString(reductions[o], NULL));
  }

  CMR_CHRMAT* wheelMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatSlice(cmr, matrix, wheelSubmatrix, &wheelMatrix) );

  CMRchrmatPrintDense(cmr, wheelMatrix, stdout, '0', true);
  
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

  ASSERT_CMR_CALL( CMRspTestTernary(cmr, matrix, NULL, reductions, &numReductions, NULL, NULL, NULL,
    DBL_MAX) );
  ASSERT_EQ( numReductions, 20U );
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

  ASSERT_CMR_CALL( CMRspTestTernary(cmr, matrix, NULL, reductions, &numReductions, NULL,
    &violatorSubmatrix, NULL, DBL_MAX) );
  ASSERT_EQ( numReductions, 10UL );
  for (size_t o = 0; o < numReductions; ++o)
  {
    printf("%s\n", CMRspReductionString(reductions[o], NULL));
  }

  CMR_CHRMAT* violatorMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatSlice(cmr, matrix, violatorSubmatrix, &violatorMatrix) );
  ASSERT_EQ( violatorMatrix->numNonzeros, 4UL );
  int det = violatorMatrix->entryValues[0] * violatorMatrix->entryValues[3]
    - violatorMatrix->entryValues[1] * violatorMatrix->entryValues[2];
  ASSERT_EQ( abs(det), 2 );

  ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, violatorMatrix, stdout, '0', true) );

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

  ASSERT_CMR_CALL( CMRspTestTernary(cmr, matrix, NULL, reductions, &numReductions, NULL,
    &violatorSubmatrix, NULL, DBL_MAX) );
  ASSERT_EQ( numReductions, 14UL );
  for (size_t o = 0; o < numReductions; ++o)
  {
    printf("%s\n", CMRspReductionString(reductions[o], NULL));
  }

  CMR_CHRMAT* violatorMatrix = NULL;
  ASSERT_CMR_CALL( CMRchrmatSlice(cmr, matrix, violatorSubmatrix, &violatorMatrix) );
  ASSERT_EQ( violatorMatrix->numNonzeros, 6UL );

  ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, violatorMatrix, stdout, '0', true) );

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &violatorMatrix) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(SeriesParallel, TernarySeparation)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "20 20 "
    "1 0 0 0  0  0  0  0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "1 0 0 0  1 -1  1 -1 0 0 0 0 0 0 0 0 1 0 0 0 "
    "0 0 1 1  0  0  0  0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 1  1  0  0  0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0  1 -1  1 -1 0 0 0 0 0 0 0 0 1 0 0 0 "
    "0 0 0 0  1 -1  1 -1 0 0 0 0 1 0 0 0 0 1 0 0 "
    "0 0 0 0 -1  1 -1  1 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0  0  0  0  0 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0  0  0  0  0 0 0 0 0 0 1 1 0 0 0 0 0 "
    "0 0 0 0  0  0  0  0 0 0 0 0 0 0 1 1 0 0 0 0 "
    "0 0 0 0  0  0  0  0 0 0 0 0 0 0 0 1 1 0 0 1 "
    "0 0 0 0  0  0  0  0 0 0 0 0 0 0 0 0 0 0 1 1 "
    "0 0 0 0  0  0  0  0 0 0 0 0 0 0 0 0 1 1 1 0 "
    "0 0 0 0  1  1  1  0 1 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0  0  0  0  0 1 1 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0  0  0  0  0 0 0 1 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0  0  0  0  0 0 1 1 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0  1  0  1  0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0  0  1  1  1 0 0 0 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0  0  1  1  0 0 0 0 0 0 0 0 0 0 0 0 0 "
  ) );

  CMR_SP_REDUCTION reductions[20];
  size_t numReductions;
  CMR_SUBMAT* violatorSubmatrix = NULL;
  CMR_SEPA* sepa = NULL;

  ASSERT_CMR_CALL( CMRspDecomposeTernary(cmr, matrix, NULL, reductions, SIZE_MAX, &numReductions, NULL,
    &violatorSubmatrix, &sepa, NULL, DBL_MAX) );

  ASSERT_FALSE( violatorSubmatrix );

  // TODO: CHECK flags.

  ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}


TEST(SeriesParallel, TernaryBadSeparation)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "20 20 "
    "1 0 0 0 0  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "1 0 0 0 1  1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 "
    "0 0 1 1 0  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 1 1  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 1  1 1 1 0 0 0 0 0 0 0 0 1 0 0 0 "
    "0 0 0 0 1 -1 1 1 0 0 0 0 1 0 0 0 0 1 0 0 "
    "0 0 0 0 1  1 1 1 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0 0  0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 "
    "0 0 0 0 0  0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 "
    "0 0 0 0 0  0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 "
    "0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 "
    "0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 "
    "0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 "
    "0 0 0 0 1  1 1 0 1 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0  0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0  0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0  0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 1  0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0  1 1 1 0 0 0 1 0 0 0 0 0 0 0 0 "
    "0 0 0 0 0  1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 "
  ) );

  CMR_SP_REDUCTION reductions[20];
  size_t numReductions;
  CMR_SUBMAT* reducedSubmatrix = NULL;
  CMR_SUBMAT* violatorSubmatrix = NULL;
  CMR_SEPA* sepa = NULL;

  ASSERT_CMR_CALL( CMRspDecomposeTernary(cmr, matrix, NULL, reductions, SIZE_MAX, &numReductions,
    &reducedSubmatrix, &violatorSubmatrix, &sepa, NULL, DBL_MAX) );

  ASSERT_TRUE( violatorSubmatrix );
  ASSERT_FALSE( sepa );
  CMR_CHRMAT* reducedMatrix = NULL;
  CMR_CHRMAT* violatorMatrix = NULL;

  ASSERT_CMR_CALL( CMRchrmatSlice(cmr, matrix, reducedSubmatrix, &reducedMatrix) );

  ASSERT_CMR_CALL( CMRchrmatSlice(cmr, reducedMatrix, violatorSubmatrix, &violatorMatrix) );

  CMRchrmatPrintDense(cmr, violatorMatrix, stdout, '0', true);

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &violatorMatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &reducedMatrix) );
  ASSERT_CMR_CALL( CMRsepaFree(cmr, &sepa) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
  ASSERT_CMR_CALL( CMRsubmatFree(cmr, &reducedSubmatrix) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
