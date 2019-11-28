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

template <typename V>
tu::SparseMatrix<V> stringToSparseMatrix(const std::string& str)
{
  std::istringstream stream(str);
  std::size_t numRows, numColumns;
  stream >> numRows >> numColumns;

  std::vector<std::size_t> begin;
  begin.reserve(numRows);
  std::vector<std::size_t> columns;
  std::vector<V> values;
  for (std::size_t row = 0; row < numRows; ++row)
  {
    begin.push_back(columns.size());
    for (std::size_t column = 0; column < numColumns; ++column)
    {
      V value;
      stream >> value;
      if (value != 0)
      {
        columns.push_back(column);
        values.push_back(value);
      }
    }
  }

  return tu::SparseMatrix<V>(true, numRows, numColumns, values.size(), &begin[0], nullptr,
    &columns[0], &values[0], true, true);
}
