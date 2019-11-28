#include <gtest/gtest.h>

#include <tu/matrix.hpp>

#include "common.hpp"

template<typename M>
void testMatrixBasic(const M& matrix)
{
  std::size_t i;

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

  struct NZ
  {
    std::size_t row;
    std::size_t column;
    int value;

    NZ(std::size_t r, std::size_t c, int v)
      : row(r), column(c), value(v)
    {

    }
  };

  NZ row0[2] = { NZ(0, 0, 4), NZ(0, 2, 5) };
  i = 0;
  for (auto nz : matrix.iterateRowNonzeros(0))
  {
    ASSERT_EQ(nz.row, 0);
    ASSERT_EQ(nz.column, row0[i].column);
    ASSERT_EQ(nz.value, row0[i].value);
    ++i;
  }

  NZ row1[1] = { NZ(1, 1, -6) };
  i = 0;
  for (auto nz : matrix.iterateRowNonzeros(1))
  {
    ASSERT_EQ(nz.row, 1);
    ASSERT_EQ(nz.column, row1[i].column);
    ASSERT_EQ(nz.value, row1[i].value);
    ++i;
  }

  // Test column iteration.

  NZ column0[1] = { NZ(0, 0, 4) };
  i = 0;
  for (auto nz : matrix.iterateColumnNonzeros(0))
  {
    ASSERT_EQ(nz.row, column0[i].row);
    ASSERT_EQ(nz.column, 0);
    ASSERT_EQ(nz.value, column0[i].value);
    ++i;
  }

  NZ column1[1] = { NZ(1, 1, -6) };
  i = 0;
  for (auto nz : matrix.iterateColumnNonzeros(1))
  {
    ASSERT_EQ(nz.row, column1[i].row);
    ASSERT_EQ(nz.column, 1);
    ASSERT_EQ(nz.value, column1[i].value);
    ++i;
  }

  NZ column2[1] = { NZ(0, 2, 5) };
  i = 0;
  for (auto nz : matrix.iterateColumnNonzeros(2))
  {
    ASSERT_EQ(nz.row, column2[i].row);
    ASSERT_EQ(nz.column, 2);
    ASSERT_EQ(nz.value, column2[i].value);
    ++i;
  }
}

TEST(DenseMatrix, Basic)
{
  testMatrixBasic(stringToDenseMatrix<int>("2 3 "
    "4 0 5 "
    "0 -6 0 "));
}


TEST(SparseMatrix, Basic)
{
  testMatrixBasic(stringToSparseMatrix<int>("2 3 "
    "4 0 5 "
    "0 -6 0 "));
}
