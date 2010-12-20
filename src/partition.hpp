/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef PARTITION_HPP_
#define PARTITION_HPP_

#include <vector>

#include "matrix_transposed.hpp"
#include "matrix_permuted.hpp"
#include "binary_linear_space.hpp"

namespace tu
{

  typedef std::pair <size_t, size_t> size_pair_t;

  namespace detail
  {

    /**
     * Copies a partial row of a matrix into a vector.
     *
     * @param matrix A given matrix
     * @param vector Vector to be filled
     * @param row Index of the row
     * @param first_column Index of the first column
     * @param beyond_column Index of the column beyond the partial row
     */

    template <typename MatrixType, typename VectorType>
    inline void copy_partial_row(const MatrixType& matrix, VectorType& vector, size_t row, size_t first_column, size_t beyond_column)
    {
      size_t j = 0;
      for (size_t i = first_column; i < beyond_column; ++i, ++j)
      {
        vector[j] = matrix(row, i);
      }
    }

    /**
     * Copies a partial column of a matrix into a vector.
     *
     * @param matrix A given matrix
     * @param vector Vector to be filled
     * @param column Index of the column
     * @param first_row Index of the first row
     * @param beyond_row Index of the row beyond the partial column
     */

    template <typename MatrixType, typename VectorType>
    inline void copy_partial_column(const MatrixType& matrix, VectorType& vector, size_t column, size_t first_row, size_t beyond_row)
    {
      copy_partial_row(matrix_transposed <const MatrixType> (matrix), vector, column, first_row, beyond_row);
    }
  }

  /**
   * Types of distributions of the rank 2
   */

  enum rank_distribution
  {
    RANK_TOO_HIGH = 0, RANK_BL_1_TR_1 = 1, RANK_BL_2_TR_0 = 2, RANK_BL_0_TR_2 = 3
  };

  /**
   * Partitioning routine. Takes a 3-separation of a minor and tries to enlarge it
   * to a 3-separation of the whole matroid.
   *
   * @param matrix Representation matrix of the matroid
   * @param top_left_height Number of rows in the top-left separation part
   * @param top_left_width Number of columns in the top-left separation part
   * @param bottom_right_height Number of rows in the bottom-right separation part
   * @param bottom_right_width Number of columns in the bottom-right separation part
   * @return A rank distribution, if enlarging was possible or RANK_TOO_HIGH if not
   */

  template <typename MatrixType>
  inline rank_distribution partition(matrix_permuted <const MatrixType>& matrix, size_t& top_left_height, size_t& top_left_width,
      size_t& bottom_right_height, size_t& bottom_right_width)
  {
    const size_t height = matrix.size1();
    size_t free_rows_first = top_left_height;
    size_t free_rows_beyond = height - bottom_right_height;

    const size_t width = matrix.size2();
    size_t free_columns_first = top_left_width;
    size_t free_columns_beyond = width - bottom_right_width;

    rank_distribution result = RANK_TOO_HIGH;

    /// Repeat until no rows/columns can be shifted.
    bool changed = true;
    while (changed)
    {
      /// Setup bottom left rows
      binary_linear_space bottom_left_row_space(top_left_width);
      std::vector <bool> bottom_left_row(top_left_width);

      for (size_t row = free_rows_beyond; row < height; ++row)
      {
        detail::copy_partial_row(matrix, bottom_left_row, row, 0, top_left_width);
        bottom_left_row_space.insert_checked(bottom_left_row);
      }

      /// Setup top right rows
      binary_linear_space top_right_row_space(bottom_right_width);
      std::vector <bool> top_right_row(bottom_right_width);

      for (size_t row = 0; row < top_left_height; ++row)
      {
        detail::copy_partial_row(matrix, top_right_row, row, free_columns_beyond, width);
        top_right_row_space.insert_checked(top_right_row);
      }

      /// Check we have rank sum of 2
      if (bottom_left_row_space.dimension() + top_right_row_space.dimension() != 2)
        return RANK_TOO_HIGH;

      /// Check all other rows
      changed = false;
      for (size_t row = free_rows_first; row < free_rows_beyond; ++row)
      {
        detail::copy_partial_row(matrix, top_right_row, row, free_columns_beyond, width);
        if (!top_right_row_space.is_spanned(top_right_row))
        {
          //          std::cout << "Row " << row << " is not spanned by top-right." << std::endl;
          detail::copy_partial_row(matrix, bottom_left_row, row, 0, top_left_width);
          if (!bottom_left_row_space.is_spanned(bottom_left_row))
            return RANK_TOO_HIGH;

          ++bottom_right_height;
          --free_rows_beyond;
          matrix_permute1(matrix, row, free_rows_beyond);
          --row;
          changed = true;
        }
      }

      /// Setup top right columns
      binary_linear_space top_right_column_space(top_left_height);
      std::vector <bool> top_right_column(top_left_height);

      for (size_t column = free_columns_beyond; column < width; ++column)
      {
        detail::copy_partial_column(matrix, top_right_column, column, 0, top_left_height);
        top_right_column_space.insert_checked(top_right_column);
      }

      /// Setup bottom left columns
      binary_linear_space bottom_left_column_space(bottom_right_height);
      std::vector <bool> bottom_left_column(bottom_right_height);

      for (size_t column = 0; column < top_left_width; ++column)
      {
        detail::copy_partial_column(matrix, bottom_left_column, column, free_rows_beyond, height);
        bottom_left_column_space.insert_checked(bottom_left_column);
      }

      /// Assert we have rank sum of 2
      assert (bottom_left_row_space.dimension () + top_right_row_space.dimension () == 2);

      /// Check all other columns
      for (size_t column = free_columns_first; column < free_columns_beyond; ++column)
      {
        detail::copy_partial_column(matrix, bottom_left_column, column, free_rows_beyond, height);
        if (!bottom_left_column_space.is_spanned(bottom_left_column))
        {
          detail::copy_partial_column(matrix, top_right_column, column, 0, top_left_height);
          if (!top_right_column_space.is_spanned(top_right_column))
            return RANK_TOO_HIGH;

          ++bottom_right_width;
          --free_columns_beyond;
          matrix_permute2(matrix, column, free_columns_beyond);
          --column;
          changed = true;
        }
      }

      switch (bottom_left_column_space.dimension())
      {
      case 2:
        result = RANK_BL_2_TR_0;
      break;
      case 1:
        result = RANK_BL_1_TR_1;
      break;
      case 0:
        result = RANK_BL_0_TR_2;
      break;
      default:
        assert (false);
      }
    }

    top_left_height = height - bottom_right_height;
    top_left_width = width - bottom_right_width;

    return result;
  }

}

#endif /* PARTITION_HPP_ */
