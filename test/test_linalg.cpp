#include <gtest/gtest.h>

#include "common.h"
#include "../src/cmr/linalg.h"
#include <vector>
#include <stdio.h>

TEST(Linalg, ComputeUpperDiagonal)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  {
    // Create a random dense matrix.
    size_t numRows = 10;
    size_t numColumns = 16;
    double nonzeroProbability = 0.5;
    int maxEntry = 9;
    std::vector<std::vector<int>> denseMatrix(numRows);
    size_t numNonzeros = 0;
    for (size_t row = 0; row < numRows; ++row)
    {
      denseMatrix[row].resize(numColumns, 0);
      for (size_t column = 0; column < numColumns; ++column)
      {
        if ((rand() * 1.0 / RAND_MAX) < nonzeroProbability)
        {
          int x = (int)((rand() * 1.0 / RAND_MAX) * 2.0 * maxEntry) - maxEntry;
          if (x >= 0)
            ++x;
          denseMatrix[row][column] = x;
          ++numNonzeros;
        }
      }
    }

    // Create sparse representation.
    CMR_INTMAT* mat = NULL;
    ASSERT_CMR_CALL( CMRintmatCreate(cmr, &mat, numRows, numColumns, numNonzeros) );
    mat->rowSlice[0] = 0;
    numNonzeros = 0;
    for (size_t row = 0; row < mat->numRows; ++row)
    {
      for (size_t column = 0; column < mat->numColumns; ++column)
      {
        int x = denseMatrix[row][column];
        if (x == 0)
          continue;
        mat->entryColumns[numNonzeros] = column;
        mat->entryValues[numNonzeros] = x;
        ++numNonzeros;
      }
      mat->rowSlice[row + 1] = numNonzeros;
    }

    CMR_INTMAT* transformed = NULL;
    size_t rank;
    CMR_SUBMAT* permutations = NULL;
    ASSERT_CMR_CALL( CMRintmatComputeUpperDiagonal(cmr, mat, true, &rank, &permutations, &transformed, NULL) );

//     ASSERT_EQ(rank, 4);

    ASSERT_CMR_CALL( CMRsubmatWriteToStream(cmr, permutations, mat->numRows, mat->numColumns, stdout) );
    std::cout << std::endl;
    ASSERT_CMR_CALL( CMRintmatPrintDense(cmr, transformed, stdout, '0', false) );

    CMR_INTMAT* zoomed = NULL;
    ASSERT_CMR_CALL( CMRintmatZoomSubmat(cmr, transformed, permutations, &zoomed) );
    CMRintmatSortNonzeros(cmr, zoomed);
    ASSERT_CMR_CALL( CMRintmatPrintDense(cmr, zoomed, stdout, '0', false) );
    ASSERT_CMR_CALL( CMRintmatFree(cmr, &zoomed) );

    ASSERT_CMR_CALL( CMRsubmatFree(cmr, &permutations) );
    ASSERT_CMR_CALL( CMRintmatFree(cmr, &transformed) );
    ASSERT_CMR_CALL( CMRintmatFree(cmr, &mat) );
  }

  CMRfreeEnvironment(&cmr);
}
