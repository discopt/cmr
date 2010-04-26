
//          Copyright Matthias Walter 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef MATRIX_REORDER_HPP_
#define MATRIX_REORDER_HPP_

#include "../config.h"
#include "permutations.hpp"
#include "matroid_transposed.hpp"

namespace tu {

  template <typename MatrixType>
  inline void matrix_apply_row_permutation (MatrixType& matrix, const permutation& perm)
  {
    permutation p = perm;
    for (size_t row = 0; row < matrix.size1 (); ++row)
    {
      size_t target_row = p (row);
      matrix_permute1 (matrix, row, target_row);
      p.rswap (row, target_row);
    }
  }

  template <typename MatrixType>
  inline void matrix_apply_column_permutation (MatrixType& matrix, const permutation& perm)
  {
    matroid_transposed <MatrixType> transposed (matrix);
    matrix_apply_row_permutation (transposed, perm);
  }

  template <typename MatrixType, typename Less>
  struct matrix_reorder_row_less
  {
    matrix_reorder_row_less (const MatrixType& matrix, size_t column_first, size_t column_beyond, Less less) :
      matrix_ (matrix), column_first_ (column_first), column_beyond_ (column_beyond), less_ (less)
    {

    }

    bool operator() (size_t a, size_t b)
    {
      //        std::cout << "comparing row " << a << " with " << b << std::endl;
      for (size_t column = column_first_; column < column_beyond_; ++column)
      {
        typename MatrixType::value_type va = matrix_ (a, column);
        typename MatrixType::value_type vb = matrix_ (b, column);

        if (less_ (va, vb))
        {
          //                std::cout << "column " << column << " said a < b" << std::endl;
          return true;
        }
        else if (less_ (vb, va))
        {
          //                std::cout << "column " << column << " said a > b" << std::endl;
          return false;
        }
      }

      //        std::cout << "a = b" << std::endl;
      return false;
    }

    const MatrixType& matrix_;
    size_t column_first_;
    size_t column_beyond_;
    Less less_;
  };

  template <typename MatrixType, typename ElementLess>
  inline void matrix_reorder_rows (const MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond,
      ElementLess element_less, permutation& result_permutation)
  {
    //    std::cout << "matrix_reorder_rows (with permutation)" << std::endl;
    result_permutation.reset (matrix.size1 ());
    //    std::cout << "constructing comparator" << std::endl;
    matrix_reorder_row_less <MatrixType, ElementLess> less (matrix, column_first, column_beyond, element_less);
    //    result_permutation.write (std::cout);
    //    std::cout << std::endl;
    //    std::cout << "sorting in [" << row_first << ", " << row_beyond << ")" << std::endl;
    sort (result_permutation, row_first, row_beyond, less);
    //    std::cout << "sorted" << std::endl;
  }

  template <typename MatrixType, typename ElementLess>
  inline void matrix_reorder_columns (const MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond,
      ElementLess element_less, permutation& result_permutation)
  {
    const matrix_transposed <MatrixType> transposed (matrix);
    matrix_reorder_rows (transposed, column_first, column_beyond, row_first, row_beyond, element_less, result_permutation);
  }

  template <typename MatrixType, typename ElementLess>
  inline void matrix_reorder_rows (MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond,
      ElementLess element_less)
  {
    permutation perm;
    //    std::cout << "matrix_reorder_rows" << std::endl;
    matrix_reorder_rows (matrix, row_first, row_beyond, column_first, column_beyond, element_less, perm);
    //    std::cout << "applying permutation: " << std::endl;
    //    perm.write(std::cout);
    //    std::cout << std::endl;
    //    perm.revert();
    matrix_apply_row_permutation (matrix, perm);
    //    std::cout << "applied" << std::endl;
  }

  template <typename MatrixType, typename ElementLess>
  inline void matrix_reorder_columns (MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond,
      ElementLess element_less)
  {
    matrix_transposed <MatrixType> transposed (matrix);
    matrix_reorder_rows (matrix, column_first, column_beyond, row_first, row_beyond, element_less);
  }

}

#endif /* MATRIX_REORDER_HPP_ */
