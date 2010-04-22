/*
 * total_unimodularity.cpp
 *
 *  Created on: Dec 20, 2009
 *      Author: xammy
 */

#include "../config.h"
#include "algorithm.hpp"
#include "total_unimodularity.hpp"
#include "matroid.hpp"
#include "violator_search.hpp"

#include <boost/numeric/ublas/io.hpp>

namespace tu {

  /// Returns true, iff the given matrix is totally unimodular.

  bool is_totally_unimodular (const integer_matrix& input_matrix)
  {
    std::cout << "Checking for a -1/0/1 matrix" << std::endl;

    if (!is_zero_plus_minus_one_matrix (input_matrix))
    {
      return false;
    }

    std::cout << "Checking for signedness" << std::endl;

    if (!is_signed_matrix (input_matrix))
    {
      //        submatrix_indices submatrix;
      //        is_signed_matrix (input_matrix, submatrix);
      //
      //        for (size_t col = 0; col < submatrix.columns.size (); col++)
      //        {
      //            std::cout << "column " << col << " = " << submatrix.columns[col] << std::endl;
      //        }
      //        for (size_t row = 0; row < submatrix.rows.size (); row++)
      //        {
      //            std::cout << "row " << row << " = " << submatrix.rows[row] << std::endl;
      //        }
      //
      //        const boost::numeric::ublas::matrix_indirect<const boost::numeric::ublas::matrix<int>,
      //                submatrix_indices::indirect_array_type> indirect_matrix (input_matrix, submatrix.rows,
      //                submatrix.columns);
      //
      //        boost::numeric::ublas::matrix<float> matrix (indirect_matrix);
      //        
      //        std::cout << matrix << std::endl;

      return false;
    }

    std::cout << "Copying matrix and creating matroid" << std::endl;

    integer_matrix matrix (input_matrix);
    integer_matroid matroid (matrix.size1 (), matrix.size2 ());

    support_matrix (matrix);

    std::pair <bool, decomposed_matroid*> result = decompose_binary_matroid (matroid, matrix, false);

    assert (result.second == NULL);

    return result.first;
  }

  /// Returns true, iff the given matrix is totally unimodular.
  /// decomposition points to a k-sum-decomposition (k=1,2,3) in graphic, cographic, R10 and maybe irregular components.

  bool is_totally_unimodular (const integer_matrix& input_matrix, decomposed_matroid*& decomposition)
  {
    if (!is_zero_plus_minus_one_matrix (input_matrix) || !is_signed_matrix (input_matrix))
    {
      decomposition = NULL;
      return false;
    }

    integer_matrix matrix (input_matrix);
    integer_matroid matroid (matrix.size1 (), matrix.size2 ());

    support_matrix (matrix);

    std::pair <bool, decomposed_matroid*> result = decompose_binary_matroid (matroid, matrix, true);
    decomposition = result.second;

    return result.first;
  }

  /// Returns true, iff the given matrix is totally unimodular.
  /// decomposition points to a k-sum-decomposition (k=1,2,3) in graphic, cographic, R10 and maybe irregular components.

  bool is_totally_unimodular (const integer_matrix& input_matrix, decomposed_matroid*& decomposition,
      submatrix_indices& violator)
  {
    /// Check each entry 
    std::pair <unsigned int, unsigned int> entry;
    if (!is_zero_plus_minus_one_matrix (input_matrix, entry))
    {
      violator.rows = submatrix_indices::indirect_array_type (1);
      violator.rows[0] = entry.first;
      violator.columns = submatrix_indices::indirect_array_type (1);
      violator.columns[0] = entry.second;
      decomposition = NULL;
      return false;
    }

    /// Signing test
    if (!is_signed_matrix (input_matrix, violator))
    {
      return false;
    }

    integer_matroid matroid (input_matrix.size1 (), input_matrix.size2 ());
    integer_matrix matrix (input_matrix);

    /// Remove sign from matrix
    support_matrix (matrix);

    /// Matroid decomposition
    std::pair <bool, decomposed_matroid*> result = decompose_binary_matroid (matroid, matrix, true);
    decomposition = result.second;
    if (result.first)
    {
      return true;
    }

    detail::search_violator (input_matrix, violator);
    return false;
  }

  /// Returns true, iff the given matrix is totally unimodular.
  /// If this is not the case, violator describes a violating submatrix.  

  bool is_totally_unimodular (const integer_matrix& matrix, submatrix_indices& violator)
  {
    if (is_totally_unimodular (matrix))
      return true;

    detail::search_violator (matrix, violator);
    return false;
  }

}
