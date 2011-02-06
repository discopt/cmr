/*
 * linear_algebra.hpp
 *
 *  Created on: Feb 5, 2011
 *      Author: xammy
 */

#ifndef LINEAR_ALGEBRA_HPP_
#define LINEAR_ALGEBRA_HPP_

#include "gcd.hpp"

namespace unimod
{
  template <typename Matrix>
  bool matrix_row_gcd(Matrix& matrix, size_t row1, size_t row2, size_t column)
  {
    if (matrix(row2, column) == 0)
      return false;

    int s, t;
    int g = gcd(matrix(row1, column), matrix(row2, column), s, t);
    int a = matrix(row2, column) / g;
    int b = -matrix(row1, column) / g;
    for (size_t c = 0; c < matrix.size2(); ++c)
    {
      int x = matrix(row1, c);
      int y = matrix(row2, c);
      matrix(row1, c) = s * x + t * y;
      matrix(row2, c) = a * x + b * y;
    }

    return true;
  }

  template <typename Matrix>
  bool matrix_column_gcd(Matrix& matrix, size_t column1, size_t column2, size_t row)
  {
    matrix_transposed <Matrix> transposed = make_transposed_matrix(matrix);
    return matrix_row_gcd(transposed, column1, column2, row);
  }
}

#endif /* LINEAR_ALGEBRA_HPP_ */
