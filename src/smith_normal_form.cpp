/*
 * smith_normal_form.cpp
 *
 *  Created on: Feb 3, 2011
 *      Author: xammy
 */

#include "smith_normal_form.hpp"

#include "matrix.hpp"
#include "matrix_permuted.hpp"
#include "matrix_transposed.hpp"
#include "linear_algebra.hpp"

#include <assert.h>
#include <set>

#include <boost/numeric/ublas/io.hpp>

namespace unimod
{
  template <typename Matrix>
  bool find_smallest_nonzero_matrix_entry(const Matrix& matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond,
      size_t& row, size_t& column)
  {
    bool result = false;
    int current_value = 0;
    for (size_t r = row_first; r != row_beyond; ++r)
    {
      for (size_t c = column_first; c != column_beyond; ++c)
      {
        int value = matrix(r, c);
        if (value == 0)
          continue;

        value = value >= 0 ? value : -value;

        if (!result || value < current_value)
        {
          result = true;
          row = r;
          column = c;
          current_value = value;
        }
      }
    }
    return result;
  }

  void smith_normal_form_diagonal(const integer_matrix& input_matrix, std::vector <int>& diagonal)
  {
    integer_matrix matrix = input_matrix;
    matrix_permuted <integer_matrix> permuted_matrix(matrix);
    size_t handled = 0;
    size_t row = 0;
    size_t column = 0;
    diagonal.resize(matrix.size1() < matrix.size2() ? matrix.size1() : matrix.size2(), 0);

    while (find_smallest_nonzero_matrix_entry(permuted_matrix, handled, permuted_matrix.size1(), handled, permuted_matrix.size2(), row, column))
    {
      matrix_permute1(permuted_matrix, handled, row);
      matrix_permute2(permuted_matrix, handled, column);

      /// Ensure it is positive
      if (permuted_matrix(handled, handled) < 0)
      {
        for (size_t r = 0; r < permuted_matrix.size1(); ++r)
          permuted_matrix(r, handled) = -permuted_matrix(r, handled);
      }
      diagonal[handled] = permuted_matrix(handled, handled);

      enum tests_t
      {
        COLUMN_ZERO = 1, ROW_ZERO = 2, GCD = 4,
      };

      int tests = COLUMN_ZERO | ROW_ZERO;
      while (true)
      {
        /// Reduce entries below pivot
        if ((tests & COLUMN_ZERO) && !matrix_column_zero(permuted_matrix, handled, handled + 1, permuted_matrix.size1()))
        {
          bool changed = false;
          for (size_t r = handled + 1; r < matrix.size1(); ++r)
          {
            changed = matrix_row_gcd(permuted_matrix, handled, r, handled) || changed;
          }
          if (changed)
          {
            tests = ROW_ZERO;
            continue;
          }
        }
        /// Reduce entries right of pivot
        if ((tests & ROW_ZERO) && !matrix_row_zero(permuted_matrix, handled, handled + 1, permuted_matrix.size2()))
        {
          bool changed = false;
          for (size_t c = handled + 1; c < permuted_matrix.size2(); ++c)
          {
            changed = matrix_column_gcd(permuted_matrix, handled, c, handled) || changed;
          }
          if (changed)
          {
            tests = COLUMN_ZERO;
            continue;
          }
        }
        /// Test remaining for gcd and move to column
        bool changed = false;
        for (size_t r = handled + 1; r < permuted_matrix.size1(); ++r)
        {
          for (size_t c = handled + 1; c < permuted_matrix.size2(); ++c)
          {
            if (gcd(permuted_matrix(r, c), permuted_matrix(handled, handled)) != permuted_matrix(handled, handled))
            {
              for (size_t r = handled; r < permuted_matrix.size1(); ++r)
              {
                permuted_matrix(r, handled) += permuted_matrix(r, c);
              }
              changed = true;
              break;
            }
          }
          if (changed)
            break;
        }
        if (changed)
          tests = COLUMN_ZERO;
        else
          break;
      }

      handled++;
    }

  }
}
