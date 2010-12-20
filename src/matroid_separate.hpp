/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef MATROID_SEPARATE_HPP_
#define MATROID_SEPARATE_HPP_

namespace tu
{

  template <typename MatroidType, typename MatrixType>
  void matroid_separate(MatroidType& matroid, MatrixType& matrix, const int type, const std::pair <size_t, size_t>& split,
      integer_matroid& upper_left_matroid, integer_matrix& upper_left_matrix, integer_matroid& lower_right_matroid,
      integer_matrix& lower_right_matrix)
  {
    if (type == 1)
    {
      // Upper left
      upper_left_matroid.resize(split.first, split.second);
      upper_left_matrix.resize(split.first, split.second, false);
      for (size_t row = 0; row < split.first; ++row)
        upper_left_matroid.name1(row) = matroid.name1(row);
      for (size_t column = 0; column < split.second; ++column)
        upper_left_matroid.name2(column) = matroid.name2(column);
      for (size_t row = 0; row < split.first; ++row)
      {
        for (size_t column = 0; column < split.second; ++column)
          upper_left_matrix(row, column) = matrix(row, column);
      }

      // Lower right
      lower_right_matroid.resize(matroid.size1() - split.first, matroid.size2() - split.second);
      lower_right_matrix.resize(matrix.size1() - split.first, matrix.size2() - split.second);
      for (size_t row = 0; row < lower_right_matroid.size1(); ++row)
        lower_right_matroid.name1(row) = matroid.name1(split.first + row);
      for (size_t column = 0; column < lower_right_matroid.size2(); ++column)
        lower_right_matroid.name2(column) = matroid.name2(split.second + column);
      for (size_t row = 0; row < lower_right_matrix.size1(); ++row)
      {
        for (size_t column = 0; column < lower_right_matrix.size2(); ++column)
          lower_right_matrix(row, column) = matrix(split.first + row, split.second + column);
      }
    }
    else
    {
      throw std::logic_error("Not yet implemented!");
    }
  }

}

#endif /* MATROID_SEPARATE_HPP_ */
