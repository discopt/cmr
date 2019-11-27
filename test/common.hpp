#pragma once

#include <string>
#include <tu/matrix.hpp>

template <typename V>
tu::DenseMatrix<V> stringToDenseMatrix(const std::string& str)
{
  std::istringstream stream(str);
  std::size_t numRows, numColumns;
  stream >> numRows >> numColumns;

  tu::DenseMatrix<V> result(numRows, numColumns);
  for (std::size_t row = 0; row < numRows; ++row)
  {
    for (std::size_t column = 0; column < numColumns; ++column)
    {
      V entry;
      stream >> entry;
      result.set(row, column, entry);
    }
  }

  return result;
}
