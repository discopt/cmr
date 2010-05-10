/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef MATROID_REORDER_HPP_
#define MATROID_REORDER_HPP_

#include "matrix_reorder.hpp"

namespace tu {

  /**
   * Applies a row permutation to a given matroid and its representation matrix, i.e. performs the
   * necessary row swaps at the matroid / matrix directly.
   *
   * @param matroid The given matroid to be changed
   * @param matrix Representation matrix of the given matroid
   * @param perm The given permutation
   */

  template <typename MatroidType, typename MatrixType>
  inline void matroid_apply_row_permutation (MatroidType& matroid, MatrixType& matrix, const permutation& perm)
  {
    permutation p = perm;
    for (size_t row = 0; row < matroid.size1(); ++row)
    {
      size_t target_row = p(row);

      matroid_permute1(matroid, matrix, row, target_row);
      p.rswap(row, target_row);
    }
  }

  /**
   * Applies a column permutation to a given matroid and its representation matrix, i.e. performs the
   * necessary columns swaps at the matroid / matrix directly.
   *
   * @param matroid The given matroid to be changed
   * @param matrix Representation matrix of the given matroid
   * @param perm The given permutation
   */

  template <typename MatroidType, typename MatrixType>
  inline void matroid_apply_column_permutation (MatroidType& matroid, MatrixType& matrix, const permutation& perm)
  {
    matroid_transposed <MatroidType> transposed_matroid(matroid);
    matrix_transposed <MatrixType> transposed_matrix(matrix);
    matroid_apply_row_permutation(transposed_matroid, transposed_matrix, perm);
  }

  /**
   * Reorders the rows of a matroid and its representation matrix by setting up a permutation to view the ordering.
   *
   * @param matroid The given matrix, which is unchanged.
   * @param matrix Representation matrix of the given matroid
   * @param row_first First row to reorder
   * @param row_beyond Beyond row to reorder
   * @param column_first First column for comparison
   * @param column_beyond Beyond column for comparison
   * @param element_less Compare functor
   */

  template <typename MatroidType, typename MatrixType, typename ElementLess>
  inline void matroid_reorder_rows (MatroidType& matroid, MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first,
      size_t column_beyond, ElementLess element_less)
  {
    permutation perm;
    matrix_reorder_rows(matrix, row_first, row_beyond, column_first, column_beyond, element_less, perm);
    matroid_apply_row_permutation(matroid, matrix, perm);
  }

  /**
   * Reorders the columns of a matroid and its representation matrix by setting up a permutation to view the ordering.
   *
   * @param matroid The given matrix, which is unchanged.
   * @param matrix Representation matrix of the given matroid
   * @param row_first First row for comparison
   * @param row_beyond Beyond row for comparison
   * @param column_first First column to reorder
   * @param column_beyond Beyond column to reorder
   * @param element_less Compare functor
   */

  template <typename MatroidType, typename MatrixType, typename ElementLess>
  inline void matroid_reorder_columns (MatroidType& matroid, MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first,
      size_t column_beyond, ElementLess element_less)
  {
    matroid_transposed <MatroidType> transposed_matroid(matroid);
    matrix_transposed <MatrixType> transposed_matrix(matrix);
    matroid_reorder_rows(transposed_matroid, transposed_matrix, column_first, column_beyond, row_first, row_beyond, element_less);
  }

}

#endif /* MATROID_REORDER_HPP_ */
