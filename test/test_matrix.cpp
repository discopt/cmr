#include <gtest/gtest.h>

#include <tu/matrix.hpp>

#include "common.hpp"

TEST(DenseMatrix, Basic)
{
  std::size_t i;
  auto matrix = stringToDenseMatrix<int>("2 3 "
    "4 0 5 "
    "0 -6 0 ");

  ASSERT_EQ(matrix.numRows(), 2);
  ASSERT_EQ(matrix.numColumns(), 3);

  // Test entry-wise.
  
  ASSERT_EQ(matrix.get(0, 0), 4);
  ASSERT_EQ(matrix.get(0, 1), 0);
  ASSERT_EQ(matrix.get(0, 2), 5);
  ASSERT_EQ(matrix.get(1, 0), 0);
  ASSERT_EQ(matrix.get(1, 1), -6);
  ASSERT_EQ(matrix.get(1, 2), 0);

  // Test row/column count.

  ASSERT_EQ(matrix.countRowNonzeros(0), 2);
  ASSERT_EQ(matrix.countRowNonzeros(1), 1);
  ASSERT_EQ(matrix.countColumnNonzeros(0), 1);
  ASSERT_EQ(matrix.countColumnNonzeros(1), 1);
  ASSERT_EQ(matrix.countColumnNonzeros(2), 1);

  // Test row iteration.

  tu::Matrix<int>::Nonzero row0[2] = { tu::Matrix<int>::Nonzero(0, 0, 4), tu::Matrix<int>::Nonzero(0, 2, 5) };
  i = 0;
  for (auto nz : matrix.iterateRowNonzeros(0))
  {
    ASSERT_EQ(nz.row, 0);
    ASSERT_EQ(nz.column, row0[i].column);
    ASSERT_EQ(nz.value, row0[i].value);
    ++i;
  }

  tu::Matrix<int>::Nonzero row1[1] = { tu::Matrix<int>::Nonzero(1, 1, -6) };
  i = 0;
  for (auto nz : matrix.iterateRowNonzeros(1))
  {
    ASSERT_EQ(nz.row, 1);
    ASSERT_EQ(nz.column, row1[i].column);
    ASSERT_EQ(nz.value, row1[i].value);
    ++i;
  }

  // Test column iteration.

  tu::Matrix<int>::Nonzero column0[1] = { tu::Matrix<int>::Nonzero(0, 0, 4) };
  i = 0;
  for (auto nz : matrix.iterateColumnNonzeros(0))
  {
    ASSERT_EQ(nz.row, column0[i].row);
    ASSERT_EQ(nz.column, 0);
    ASSERT_EQ(nz.value, column0[i].value);
    ++i;
  }

  tu::Matrix<int>::Nonzero column1[1] = { tu::Matrix<int>::Nonzero(1, 1, -6) };
  i = 0;
  for (auto nz : matrix.iterateColumnNonzeros(1))
  {
    ASSERT_EQ(nz.row, column1[i].row);
    ASSERT_EQ(nz.column, 1);
    ASSERT_EQ(nz.value, column1[i].value);
    ++i;
  }
  
  tu::Matrix<int>::Nonzero column2[1] = { tu::Matrix<int>::Nonzero(0, 2, 5) };
  i = 0;
  for (auto nz : matrix.iterateColumnNonzeros(2))
  {
    ASSERT_EQ(nz.row, column2[i].row);
    ASSERT_EQ(nz.column, 2);
    ASSERT_EQ(nz.value, column2[i].value);
    ++i;
  }
}
