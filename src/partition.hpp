/**
 *
 * Copyright (c) 2010 Matthias Walter (xammy@xammy.homelinux.net)
 *
 * Authors: Matthias Walter
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#ifndef PARTITION_HPP_
#define PARTITION_HPP_

#include "binary_linear_space.hpp"

namespace tu {

  typedef std::pair <size_t, size_t> size_pair_t;

  namespace detail {

    template <typename MatrixType, typename VectorType>
    inline void copy_partial_row (const MatrixType& matrix, VectorType& vector, size_t row, size_t first_column, size_t beyond_column)
    {
      size_t j = 0;
      for (size_t i = first_column; i < beyond_column; ++i, ++j)
      {
        vector[j] = matrix (row, i);
      }
    }

    template <typename MatrixType, typename VectorType>
    inline void copy_partial_column (const MatrixType& matrix, VectorType& vector, size_t column, size_t first_row, size_t beyond_row)
    {
      copy_partial_row (matrix_transposed <const MatrixType> (matrix), vector, column, first_row, beyond_row);
    }
  }

  enum rank_distribution
  {
    RANK_TOO_HIGH = 0, RANK_BL_1_TR_1 = 1, RANK_BL_2_TR_0 = 2, RANK_BL_0_TR_2 = 3
  };

  template <typename MatrixType>
  inline rank_distribution partition (matrix_permuted <const MatrixType>& matrix, size_t& top_left_height, size_t& top_left_width,
      size_t& bottom_right_height, size_t& bottom_right_width)
  {
    const size_t height = matrix.size1 ();
    size_t free_rows_first = top_left_height;
    size_t free_rows_beyond = height - bottom_right_height;

    const size_t width = matrix.size2 ();
    size_t free_columns_first = top_left_width;
    size_t free_columns_beyond = width - bottom_right_width;

    //    std::cout << "Partition check on matrix:\n";
    //    matrix_print (matrix);
    //    std::cout << "Blocks are " << top_left_height << " x " << top_left_width << " and " << bottom_right_height << " x "
    //        << bottom_right_width << std::endl;

    rank_distribution result;

    bool changed = true;
    while (changed)
    {
      /// Setup bottom left rows
      binary_linear_space bottom_left_row_space (top_left_width);
      std::vector <bool> bottom_left_row (top_left_width);

      for (size_t row = free_rows_beyond; row < height; ++row)
      {
        detail::copy_partial_row (matrix, bottom_left_row, row, 0, top_left_width);
        bottom_left_row_space.insert_checked (bottom_left_row);
      }

      //      std::cout << "bottom left row rank = " << bottom_left_row_space.dimension () << std::endl;

      /// Setup top right rows
      binary_linear_space top_right_row_space (bottom_right_width);
      std::vector <bool> top_right_row (bottom_right_width);

      for (size_t row = 0; row < top_left_height; ++row)
      {
        detail::copy_partial_row (matrix, top_right_row, row, free_columns_beyond, width);
        top_right_row_space.insert_checked (top_right_row);
      }

      //      std::cout << "top right row rank = " << top_right_row_space.dimension () << std::endl;

      /// Check we have rank sum of 2
      if (bottom_left_row_space.dimension () + top_right_row_space.dimension () != 2)
        return RANK_TOO_HIGH;

      /// Check all other rows
      changed = false;
      for (size_t row = free_rows_first; row < free_rows_beyond; ++row)
      {
        detail::copy_partial_row (matrix, top_right_row, row, free_columns_beyond, width);
        if (!top_right_row_space.is_spanned (top_right_row))
        {
          //          std::cout << "Row " << row << " is not spanned by top-right." << std::endl;
          detail::copy_partial_row (matrix, bottom_left_row, row, 0, top_left_width);
          if (!bottom_left_row_space.is_spanned (bottom_left_row))
            return RANK_TOO_HIGH;

          //          std::cout << "Row " << row << " is spanned by bottom-left." << std::endl;

          ++bottom_right_height;
          --free_rows_beyond;
          matrix_permute1 (matrix, row, free_rows_beyond);
          --row;
          changed = true;

          //          std::cout << "Row " << (row + 1) << " was swapped with " << free_rows_beyond << std::endl;
          //          matrix_print (matrix);
        }
        //        else
        //        {
        //          std::cout << "Row " << row << " is spanned by top-right." << std::endl;
        //        }
      }

      /// Setup top right columns
      binary_linear_space top_right_column_space (top_left_height);
      std::vector <bool> top_right_column (top_left_height);

      for (size_t column = free_columns_beyond; column < width; ++column)
      {
        detail::copy_partial_column (matrix, top_right_column, column, 0, top_left_height);
        top_right_column_space.insert_checked (top_right_column);
      }

      //      std::cout << "top right column rank = " << top_right_column_space.dimension () << std::endl;

      /// Setup bottom left columns 
      binary_linear_space bottom_left_column_space (bottom_right_height);
      std::vector <bool> bottom_left_column (bottom_right_height);

      for (size_t column = 0; column < top_left_width; ++column)
      {
        detail::copy_partial_column (matrix, bottom_left_column, column, free_rows_beyond, height);
        bottom_left_column_space.insert_checked (bottom_left_column);
      }

      //      std::cout << "bottom left column space:\n" << bottom_left_column_space << std::endl;

      //      std::cout << "bottom left column rank = " << bottom_left_column_space.dimension () << std::endl;

      /// Assert we have rank sum of 2
      assert (bottom_left_row_space.dimension () + top_right_row_space.dimension () == 2);

      /// Check all other columns
      //      changed = false;
      for (size_t column = free_columns_first; column < free_columns_beyond; ++column)
      {
        //        std::cout << "column = " << column << ", free_columns_beyond = " << free_columns_beyond << ", height = "
        //            << height << std::endl;

        detail::copy_partial_column (matrix, bottom_left_column, column, free_rows_beyond, height);
        if (!bottom_left_column_space.is_spanned (bottom_left_column))
        {
          //          std::cout << "Column " << column << " is not spanned by bottom-left." << std::endl;
          detail::copy_partial_column (matrix, top_right_column, column, 0, top_left_height);
          if (!top_right_column_space.is_spanned (top_right_column))
            return RANK_TOO_HIGH;

          //          std::cout << "Column " << column << " is spanned by top-right." << std::endl;

          ++bottom_right_width;
          --free_columns_beyond;
          //          std::cout << "swapping columns " << column << " with " << free_columns_beyond << std::endl;
          matrix_permute2 (matrix, column, free_columns_beyond);
          --column;
          changed = true;

          //          std::cout << "Column " << (column + 1) << " was swapped with " << free_columns_beyond << std::endl;
          //          matrix_print (matrix);
        }
        //        else
        //        {
        //          std::cout << "Column " << column << " is spanned by bottom-left." << std::endl;
        //          std::cout << "space:\n" << bottom_left_column_space << std::endl;
        //          std::copy (bottom_left_column.begin (), bottom_left_column.end (), std::ostream_iterator <bool> (std::cout,
        //              " "));
        //        }
      }

      switch (bottom_left_column_space.dimension ())
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

      //      std::cout << "Repeating? changed = " << changed << std::endl;
    }

    //    std::cout << "Pre-final blocks are " << top_left_height << " x " << top_left_width << " and "
    //        << bottom_right_height << " x " << bottom_right_width << std::endl;

    top_left_height = height - bottom_right_height;
    top_left_width = width - bottom_right_width;

    //    std::cout << "After partition check on matrix:\n";
    //    matrix_print (matrix);
    //    std::cout << "Blocks are " << top_left_height << " x " << top_left_width << " and " << bottom_right_height << " x "
    //        << bottom_right_width << std::endl;

    return result;
  }

}

#endif /* PARTITION_HPP_ */
