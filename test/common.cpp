#include "common.hpp"

#include <sstream>

/*
tu::Matrix stringToMatrix(const std::string& str)
{
  tu::Matrix result;
  std::stringstream stream(str);
  std::size_t numRows, numColumns;
  stream >> numRows >> numColumns;

  std::vector<tu::Matrix::Entry> begin;
  begin.reserve(numRows);
  std::vector<tu::Matrix::Index> columns;
  std::vector<tu::Matrix::Value> values;
  for (std::size_t row = 0; row < numRows; ++row)
  {
    begin.push_back(columns.size());
    for (std::size_t column = 0; column < numColumns; ++column)
    {
        int entry;
        stream >> entry;
        if (entry != 0)
        {
          columns.push_back(column);
          values.push_back(entry);
        }
    }
  }

  result.set(true, numRows, numColumns, values.size(), &begin[0], nullptr, &columns[0], &values[0],
    true, true);

  return result;
}
*/
