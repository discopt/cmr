#include <gtest/gtest.h>

#include "common.h"
#include <tu/series_parallel.h>

TEST(SeriesParallel, Char)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);

  TU_CHRMAT* matrix = NULL;
  stringToCharMatrix(tu, &matrix, "10 10 "
    "1 0 0 0 0 0 0 0 0 0 "
    "0 1 0 0 0 0 0 0 0 0 "
    "0 0 1 0 0 0 0 0 0 0 "
    "0 0 0 1 0 0 0 0 0 0 "
    "0 0 0 0 1 1 1 0 1 0 "
    "0 0 0 0 0 1 1 0 1 0 "
    "0 0 0 0 0 0 1 0 0 0 "
    "0 0 0 0 0 0 1 1 0 0 "
    "0 0 0 0 0 0 0 0 1 -1 "
    "0 0 0 0 0 0 1 1 1 -1 "
  );

  TU_SERIES_PARALLEL operations[20];
  size_t numOperations;

  ASSERT_TU_CALL( TUfindSeriesParallel(tu, matrix, operations, &numOperations, NULL, NULL, true) );

  TUchrmatFree(tu, &matrix);
  TUfreeEnvironment(&tu);
}
