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

#ifndef VIOLATOR_SEARCH_HPP_
#define VIOLATOR_SEARCH_HPP_

#include "algorithm.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp> 

namespace tu {
  namespace detail {

    template <typename MatrixType, typename IndirectArrayType>
    void shrink_rows (const MatrixType& input_matrix, IndirectArrayType& row_indices, IndirectArrayType& column_indices)
    {
      std::cout << "Shrinking rows of input matrix:\n";
      matrix_print (input_matrix);

      std::cout << "\nCurrent matrix:\n";
      boost::numeric::ublas::matrix_indirect <const MatrixType, IndirectArrayType> indirect_matrix (input_matrix, row_indices, column_indices);
      matrix_print (indirect_matrix);

      // row_indices.
    }

    /**
     * Identifies a violator in a given matrix, which is known to be not totally unimodular. 
     * 
     * @param input_matrix The given matrix 
     * @param violator Output parameter to put the violating submatrix into.
     */

    inline void search_violator (const integer_matrix& input_matrix, submatrix_indices& violator)
    {
      
      
//      violator.rows = submatrix_indices::indirect_array_type (input_matrix.size1 ());
//      violator.columns = submatrix_indices::indirect_array_type (input_matrix.size2 ());
//
//      for (size_t row = 0; row < violator.rows.size (); ++row)
//      {
//        violator.rows[row] = row;
//      }
//
//      for (size_t column = 0; column < violator.columns.size (); ++column)
//      {
//        violator .columns[column] = column;
//      }
//
//      shrink_rows (input_matrix, violator.rows, violator.columns);
//
//      shrink_rows (view_matrix_transposed (input_matrix), violator.columns, violator.rows);

      throw std::runtime_error ("Search for violator is not yet implemented.");
    }

  }
}

#endif /* VIOLATOR_SEARCH_HPP_ */
