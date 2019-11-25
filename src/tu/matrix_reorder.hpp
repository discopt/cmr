/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef MATRIX_REORDER_HPP_
#define MATRIX_REORDER_HPP_

#include <tu/permutations.hpp>
#include "matroid_transposed.hpp"
#include <tu/matrix.hpp>

namespace tu
{

  /**
   * Applies a row permutation to a given matrix, i.e. performs the
   * necessary row swaps at the matrix directly.
   *
   * @param matrix The matrix to be changed.
   * @param perm The given permutation
   */

  template <typename MatrixType>
  inline void matrix_apply_row_permutation(MatrixType& matrix, const permutation& perm)
  {
    permutation p = perm;
    for (size_t row = 0; row < matrix.size1(); ++row)
    {
      size_t target_row = p(row);
      matrix_permute1(matrix, row, target_row);
      p.rswap(row, target_row);
    }
  }

  /**
   * Applies a columnn permutation to a given matrix, i.e. performs the
   * necessary column swaps at the matrix directly.
   *
   * @param matrix The matrix to be changed.
   * @param perm The given permutation
   */

  template <typename MatrixType>
  inline void matrix_apply_column_permutation(MatrixType& matrix, const permutation& perm)
  {
    matroid_transposed <MatrixType> transposed(matrix);
    matrix_apply_row_permutation(transposed, perm);
  }

  /**
   * A functor to compare rows of a matrix lexicographically.
   */

  template <typename MatrixType, typename Less>
  struct matrix_reorder_row_less
  {
    /**
     * Creates the functor.
     *
     * @param matrix The given matrix
     * @param column_first First column to compare at
     * @param column_beyond Beyond column to compare at
     * @param less Used functor to compare elements
     */

    matrix_reorder_row_less(const MatrixType& matrix, size_t column_first, size_t column_beyond, Less less) :
      matrix_(matrix), column_first_(column_first), column_beyond_(column_beyond), less_(less)
    {

    }

    /**
     * Applies the functor.
     *
     * @param a First row index
     * @param b Second row index
     * @return false if and only if the first row is lexicographically greater than the second.
     */

    bool operator()(size_t a, size_t b)
    {
      for (size_t column = column_first_; column < column_beyond_; ++column)
      {
        typename MatrixType::value_type va = matrix_(a, column);
        typename MatrixType::value_type vb = matrix_(b, column);

        if (less_(va, vb))
          return true;
        else if (less_(vb, va))
          return false;
      }

      return false;
    }

    const MatrixType& matrix_;
    size_t column_first_;
    size_t column_beyond_;
    Less less_;
  };

  /**
   * Reorders the rows of a matrix by setting up a permutation to view the ordering.
   *
   * @param matrix The given matrix, which is unchanged.
   * @param row_first First row to reorder
   * @param row_beyond Beyond row to reorder
   * @param column_first First column for comparison
   * @param column_beyond Beyond column for comparison
   * @param element_less Compare functor
   * @param result_permutation Resulting permutation
   */

  template <typename MatrixType, typename ElementLess>
  inline void matrix_reorder_rows(const MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond,
      ElementLess element_less, permutation& result_permutation)
  {
    result_permutation.reset(matrix.size1());
    matrix_reorder_row_less <MatrixType, ElementLess> less(matrix, column_first, column_beyond, element_less);
    sort(result_permutation, row_first, row_beyond, less);
  }

  /**
   * Reorders the columns of a matrix by setting up a permutation to view the ordering.
   *
   * @param matrix The given matrix, which is unchanged.
   * @param row_first First row for comparison
   * @param row_beyond Beyond row for comparison
   * @param column_first First column to reorder
   * @param column_beyond Beyond column to reorder
   * @param element_less Compare functor
   * @param result_permutation Resulting permutation
   */

  template <typename MatrixType, typename ElementLess>
  inline void matrix_reorder_columns(const MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond,
      ElementLess element_less, permutation& result_permutation)
  {
    const matrix_transposed <MatrixType> transposed(matrix);
    matrix_reorder_rows(transposed, column_first, column_beyond, row_first, row_beyond, element_less, result_permutation);
  }

  /**
   * Reorders the rows of a matrix.
   *
   * @param matrix The given matrix, which is unchanged.
   * @param row_first First row for comparison
   * @param row_beyond Beyond row for comparison
   * @param column_first First column to reorder
   * @param column_beyond Beyond column to reorder
   * @param element_less Compare functor
   */

  template <typename MatrixType, typename ElementLess>
  inline void matrix_reorder_rows(MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond,
      ElementLess element_less)
  {
    permutation perm;
    matrix_reorder_rows(matrix, row_first, row_beyond, column_first, column_beyond, element_less, perm);
    matrix_apply_row_permutation(matrix, perm);
  }

  /**
   * Reorders the columns of a matrix.
   *
   * @param matrix The given matrix, which is unchanged.
   * @param row_first First row for comparison
   * @param row_beyond Beyond row for comparison
   * @param column_first First column to reorder
   * @param column_beyond Beyond column to reorder
   * @param element_less Compare functor
   */

  template <typename MatrixType, typename ElementLess>
  inline void matrix_reorder_columns(MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond,
      ElementLess element_less)
  {
    matrix_transposed <MatrixType> transposed(matrix);
    matrix_reorder_rows(matrix, column_first, column_beyond, row_first, row_beyond, element_less);
  }

}

#endif /* MATRIX_REORDER_HPP_ */
