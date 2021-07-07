#include <gtest/gtest.h>

#include "common.h"
#include <tu/series_parallel.h>

TEST(SeriesParallel, Binary)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);

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

  ASSERT_TU_CALL( TUfindSeriesParallel(tu, matrix, operations, &numOperations, NULL, /*NULL, */true) );
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

  TUchrmatFree(tu, &matrix);
  TUfreeEnvironment(&tu);
}
