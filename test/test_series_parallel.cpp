#include <gtest/gtest.h>

#include "common.h"
#include <tu/series_parallel.h>

TEST(SeriesParallel, Empty)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  TU_SP operations[2];
  size_t numOperations;

  {
    TU_CHRMAT* mat0x0 = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &mat0x0, "0 0 "
    ) );

    ASSERT_TU_CALL( TUfindSeriesParallel(tu, mat0x0, operations, &numOperations, NULL, NULL, NULL, NULL, true) );
    ASSERT_EQ( numOperations, 0 );

    ASSERT_TU_CALL( TUchrmatFree(tu, &mat0x0) );
  }

  {
    TU_CHRMAT* mat2x0 = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &mat2x0, "2 0 "
    ) );

    ASSERT_TU_CALL( TUfindSeriesParallel(tu, mat2x0, operations, &numOperations, NULL, NULL, NULL, NULL, true) );
    ASSERT_EQ( numOperations, 2 );
    ASSERT_EQ( operations[0].element, -1 ); ASSERT_EQ( operations[0].mate, 0 );
    ASSERT_EQ( operations[1].element, -2 ); ASSERT_EQ( operations[1].mate, 0 );

    ASSERT_TU_CALL( TUchrmatFree(tu, &mat2x0) );
  }

  {
    TU_CHRMAT* mat0x2 = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &mat0x2, "0 2 "
    ) );

    ASSERT_TU_CALL( TUfindSeriesParallel(tu, mat0x2, operations, &numOperations, NULL, NULL, NULL, NULL, true) );
    ASSERT_EQ( numOperations, 2 );
    ASSERT_EQ( operations[0].element, 1 ); ASSERT_EQ( operations[0].mate, 0 );
    ASSERT_EQ( operations[1].element, 2 ); ASSERT_EQ( operations[1].mate, 0 );

    ASSERT_TU_CALL( TUchrmatFree(tu, &mat0x2) );
  }

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(SeriesParallel, Reduction)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "20 20 "
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

  TU_SP operations[40];
  size_t numOperations;
  TU_SUBMAT* submatrix = NULL;

  ASSERT_TU_CALL( TUfindSeriesParallel(tu, matrix, operations, &numOperations, &submatrix, NULL, NULL, NULL, true) );
  ASSERT_EQ( numOperations, 20 );
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

  ASSERT_TU_CALL( TUsubmatFree(tu, &submatrix) );
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(SeriesParallel, FirstAttemptShortWheel)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "20 20 "
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

  TU_SP operations[40];
  size_t numOperations;
  TU_SUBMAT* wheelSubmatrix = NULL;

  ASSERT_TU_CALL( TUfindSeriesParallel(tu, matrix, operations, &numOperations, NULL, &wheelSubmatrix, NULL, NULL, true) );
  ASSERT_EQ( numOperations, 8 );

  TU_CHRMAT* wheelMatrix = NULL;
  ASSERT_TU_CALL( TUchrmatFilterSubmat(tu, matrix, wheelSubmatrix, &wheelMatrix) );

  TUchrmatPrintDense(tu, stdout, wheelMatrix, '0', true);
  
  ASSERT_TU_CALL( TUchrmatFree(tu, &wheelMatrix) );
  ASSERT_TU_CALL( TUsubmatFree(tu, &wheelSubmatrix) );
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(SeriesParallel, SecondAttemptLongWheel)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "20 20 "
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

  TU_SP operations[40];
  size_t numOperations;
  TU_SUBMAT* wheelSubmatrix = NULL;

  ASSERT_TU_CALL( TUfindSeriesParallel(tu, matrix, operations, &numOperations, NULL, &wheelSubmatrix, NULL, NULL, true) );
  ASSERT_EQ( numOperations, 8 );
  for (size_t o = 0; o < numOperations; ++o)
  {
    printf("%s\n", TUspString(operations[o], NULL));
  }
  
  TU_CHRMAT* wheelMatrix = NULL;
  ASSERT_TU_CALL( TUchrmatFilterSubmat(tu, matrix, wheelSubmatrix, &wheelMatrix) );

  TUchrmatPrintDense(tu, stdout, wheelMatrix, '0', true);
  
  ASSERT_TU_CALL( TUchrmatFree(tu, &wheelMatrix) );
  ASSERT_TU_CALL( TUsubmatFree(tu, &wheelSubmatrix) );
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(SeriesParallel, SecondAttemptShortWheel)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "20 20 "
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

  TU_SP operations[40];
  size_t numOperations;
  TU_SUBMAT* wheelSubmatrix = NULL;

  ASSERT_TU_CALL( TUfindSeriesParallel(tu, matrix, operations, &numOperations, NULL, &wheelSubmatrix, NULL, NULL, true) );
  ASSERT_EQ( numOperations, 8 );
  for (size_t o = 0; o < numOperations; ++o)
  {
    printf("%s\n", TUspString(operations[o], NULL));
  }

  TU_CHRMAT* wheelMatrix = NULL;
  ASSERT_TU_CALL( TUchrmatFilterSubmat(tu, matrix, wheelSubmatrix, &wheelMatrix) );

  TUchrmatPrintDense(tu, stdout, wheelMatrix, '0', true);
  
  ASSERT_TU_CALL( TUchrmatFree(tu, &wheelMatrix) );
  ASSERT_TU_CALL( TUsubmatFree(tu, &wheelSubmatrix) );
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(SeriesParallel, Separation)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "20 20 "
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

  TU_SP operations[40];
  size_t numOperations;
  TU_SUBMAT* wheelSubmatrix = NULL;
  TU_ELEMENT separationRank1Elements[40];
  size_t numSeparationRank1Elements;

  ASSERT_TU_CALL( TUfindSeriesParallel(tu, matrix, operations, &numOperations, NULL, &wheelSubmatrix,
    separationRank1Elements, &numSeparationRank1Elements, true) );
  ASSERT_EQ( numOperations, 8 );
  ASSERT_EQ( numSeparationRank1Elements, 17 );
  ASSERT_EQ( separationRank1Elements[0], TUrowToElement(4) );
  ASSERT_EQ( separationRank1Elements[1], TUrowToElement(5) );
  ASSERT_EQ( separationRank1Elements[2], TUrowToElement(6) );
  ASSERT_EQ( separationRank1Elements[3], TUrowToElement(7) );
  ASSERT_EQ( separationRank1Elements[4], TUrowToElement(8) );
  ASSERT_EQ( separationRank1Elements[5], TUrowToElement(9) );
  ASSERT_EQ( separationRank1Elements[6], TUrowToElement(10) );
  ASSERT_EQ( separationRank1Elements[7], TUrowToElement(11) );
  ASSERT_EQ( separationRank1Elements[8], TUrowToElement(12) );
  ASSERT_EQ( separationRank1Elements[9], TUcolumnToElement(4) );
  ASSERT_EQ( separationRank1Elements[10], TUcolumnToElement(5) );
  ASSERT_EQ( separationRank1Elements[11], TUcolumnToElement(6) );
  ASSERT_EQ( separationRank1Elements[12], TUcolumnToElement(7) );
  ASSERT_EQ( separationRank1Elements[13], TUcolumnToElement(8) );
  ASSERT_EQ( separationRank1Elements[14], TUcolumnToElement(9) );
  ASSERT_EQ( separationRank1Elements[15], TUcolumnToElement(10) );
  ASSERT_EQ( separationRank1Elements[16], TUcolumnToElement(11) );
  
  ASSERT_TU_CALL( TUsubmatFree(tu, &wheelSubmatrix) );
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(SeriesParallel, ThirdAttemptAfterSeparation)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  TU_CHRMAT* matrix = NULL;
  ASSERT_TU_CALL( stringToCharMatrix(tu, &matrix, "20 20 "
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

  TU_SP operations[40];
  size_t numOperations;
  TU_SUBMAT* wheelSubmatrix = NULL;

  ASSERT_TU_CALL( TUfindSeriesParallel(tu, matrix, operations, &numOperations, NULL, &wheelSubmatrix, NULL, NULL, true) );
  ASSERT_EQ( numOperations, 8 );
  for (size_t o = 0; o < numOperations; ++o)
  {
    printf("%s\n", TUspString(operations[o], NULL));
  }

  TU_CHRMAT* wheelMatrix = NULL;
  ASSERT_TU_CALL( TUchrmatFilterSubmat(tu, matrix, wheelSubmatrix, &wheelMatrix) );

  TUchrmatPrintDense(tu, stdout, wheelMatrix, '0', true);
  
  ASSERT_TU_CALL( TUchrmatFree(tu, &wheelMatrix) );
  ASSERT_TU_CALL( TUsubmatFree(tu, &wheelSubmatrix) );
  ASSERT_TU_CALL( TUchrmatFree(tu, &matrix) );
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}
