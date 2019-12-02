#pragma once

#include "common.hpp"
#include "gcd.hpp"
#include <tu/linear_algebra.hpp>
#include <tu/matrix.hpp>
#include <tu/matrix_permuted.hpp>
#include <tu/matrix_transposed.hpp>

namespace tu
{
  template <typename Matrix>
  void smith_normal_form_diagonal(Matrix& input_matrix, std::vector <int>& diagonal)
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

} /* namespace tu */
