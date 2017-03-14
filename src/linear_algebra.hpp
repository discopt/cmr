/*
 * linear_algebra.hpp
 *
 *  Created on: Feb 5, 2011
 *      Author: xammy
 */

#ifndef LINEAR_ALGEBRA_HPP_
#define LINEAR_ALGEBRA_HPP_

#include "gcd.hpp"
#include "matrix_transposed.hpp"
#include "matrix_permuted.hpp"
#include "matrix.hpp"

namespace unimod
{
  template <typename Matrix>
  void matrix_row_combine(Matrix& matrix, size_t row1, size_t row2, long long ul, long long ur, long long ll, long long lr)
  {
    assert(abs(ul * lr - ll * ur) == 1);
    for (size_t c = 0; c < matrix.size2(); ++c)
    {
      long long x = matrix(row1, c);
      long long y = matrix(row2, c);
      if ((abs(x) > (1L << 31)) || (abs(y) > (1L << 31))
       || (abs(ul) > (1L << 31)) || (abs(ur) > (1L << 31)))
      {
        throw std::runtime_error("Numbers too large.");
      }
      matrix(row1, c) = ul * x + ur * y;
      matrix(row2, c) = ll * x + lr * y;
    }
  }

  template <typename Matrix>
  bool matrix_row_gcd(Matrix& matrix, size_t row1, size_t row2, size_t column)
  {
    if (matrix(row2, column) == 0)
      return false;

    long long x = matrix(row1, column);
    long long y = matrix(row2, column);
    if (x > (1L << 31) || x < -1 * (1L << 32) || y > (1L << 31) || y < -1 * (1L << 32))
    {
      throw std::runtime_error("Numbers too large.");
    }
    long long s, t;
    long long g = gcd(x, y, s, t);
    if (s == 0)
    {
      /// Special case: If both are |.| = g then s should not be zero.
      matrix_row_combine(matrix, row1, row2, 1, 0, -x / y, 1);
    }
    else
      matrix_row_combine(matrix, row1, row2, s, t, y / g, -x / g);

    return true;
  }

  template <typename Matrix>
  bool matrix_column_gcd(Matrix& matrix, size_t column1, size_t column2, size_t row)
  {
    matrix_transposed <Matrix> transposed = make_transposed_matrix(matrix);
    return matrix_row_gcd(transposed, column1, column2, row);
  }

  template <typename InputMatrix, typename OutputMatrix>
  size_t matrix_find_column_basis_and_transform_integral(const InputMatrix& input_matrix, OutputMatrix& output_matrix,
      std::vector <size_t>& column_basis)
  {
    output_matrix = input_matrix;
    column_basis.clear();
    matrix_permuted <OutputMatrix> permuted_matrix(output_matrix);
    size_t rank = 0;

    for (rank = 0; rank < permuted_matrix.size1(); ++rank)
    {
      /// Swap a non-zero column to the pivot column.
      bool found = false;
      for (size_t c = rank; c < permuted_matrix.size2(); ++c)
      {
        if (matrix_column_zero(permuted_matrix, c, rank, permuted_matrix.size1()))
          continue;

        found = true;
        matrix_permute2(permuted_matrix, rank, c);
        break;
      }

      if (!found)
        break;

      /// Swap the row to the pivot row (in the output matrix. This is necessary later anyway.
      for (size_t r = rank; r < permuted_matrix.size1(); ++r)
      {
        if (permuted_matrix(r, rank) != 0)
        {
          matrix_permute1(output_matrix, rank, r);
          break;
        }
      }

      /// Reduce the entries below the pivot which is at (rank, rank).
      for (size_t r = rank + 1; r < permuted_matrix.size1(); ++r)
      {
        matrix_row_gcd(permuted_matrix, rank, r, rank);
      }
    }

    for (size_t p = rank; p > 0; --p)
    {
      long long entry = permuted_matrix(p - 1, p - 1);
      assert(entry != 0);
      long long g = entry;
      for (size_t c = p; c < permuted_matrix.size2(); ++c)
        g = gcd(g, permuted_matrix(p - 1, c));
      if (entry * g < 0)
        g = -g;
      assert(g != 0);
      entry /= g;

      /// Divide row by gcd of all entries.
      if (g != 1)
      {
        for (size_t c = p - 1; c < permuted_matrix.size2(); ++c)
        {
          assert(permuted_matrix(p - 1, c) % g == 0);
          permuted_matrix(p - 1, c) /= g;
        }
      }

      /// Reduce column above entry
      for (size_t r = 0; r < p - 1; ++r)
      {
        long long factor = permuted_matrix(r, p - 1) / entry;
        for (size_t c = p - 1; c < permuted_matrix.size2(); ++c)
        {
          permuted_matrix(r, c) -= factor * permuted_matrix(p - 1, c);
        }
      }
    }

    output_matrix.resize(rank, output_matrix.size2());

    for (size_t r = 0; r < rank; ++r)
      column_basis.push_back(permuted_matrix.perm2()(r));

    return rank;
  }
}

#endif /* LINEAR_ALGEBRA_HPP_ */
