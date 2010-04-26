
//          Copyright Matthias Walter 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef TOTAL_UNIMODULARITY_HPP_
#define TOTAL_UNIMODULARITY_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace tu {

  /// Represents subsets of row/column index sets

  struct submatrix_indices
  {
    typedef boost::numeric::ublas::vector <size_t> vector_type;
    typedef boost::numeric::ublas::indirect_array <vector_type> indirect_array_type;

    indirect_array_type rows;
    indirect_array_type columns;
  };

  typedef boost::numeric::ublas::matrix <int> integer_matrix;
  typedef boost::numeric::ublas::matrix_indirect <const integer_matrix, submatrix_indices::indirect_array_type> integer_submatrix;

  class decomposed_matroid;

  /// Returns a decomposition of a given binary matroid into a k-sum-decomposition (k=1,2,3)
  /// in graphic, cographic, R10 and maybe irregular components. 

  decomposed_matroid* decompose_binary_matroid (const boost::numeric::ublas::matrix <int>& matrix);

  /// Returns true, iff the given matrix is totally unimodular.

  bool is_totally_unimodular (const integer_matrix& matrix);

  /// Returns true, iff the given matrix is totally unimodular.
  /// decomposition points to a k-sum-decomposition (k=1,2,3) in graphic, cographic, R10 and maybe irregular components.

  bool is_totally_unimodular (const integer_matrix& matrix, decomposed_matroid*& decomposition);

  /// Returns true, iff the given matrix is totally unimodular.
  /// If this is not the case, violator describes a violating submatrix.  

  bool is_totally_unimodular (const integer_matrix& matrix, submatrix_indices& violator);

  /// Returns true, iff the given matrix is totally unimodular.
  /// decomposition points to a k-sum-decomposition (k=1,2,3) in graphic, cographic, R10 and maybe irregular components.

  bool is_totally_unimodular (const integer_matrix& matrix, decomposed_matroid*& decomposition, submatrix_indices& violator);

  /// Returns true, iff the given matrix contains only values {-1, 0, 1}.

  bool is_zero_plus_minus_one_matrix (const integer_matrix& matrix);

  /// Returns true, iff the given matrix contains only values {-1, 0, 1}.
  /// If not, position will contain the position of the wrong entry.
  /// Running time: O(height * width) 

  bool is_zero_plus_minus_one_matrix (const integer_matrix& matrix, std::pair <size_t, size_t>& position);

  /// Returns true, iff the given matrix contains only values {0, 1}.

  bool is_zero_one_matrix (const integer_matrix& matrix);

  /// Returns true, iff the given matrix contains only values {0, 1}.
  /// If not, position will contain the position of the wrong entry.
  /// Running time: O(height * width) 

  bool is_zero_one_matrix (const integer_matrix& matrix, std::pair <size_t, size_t>& position);

  /// Returns true, iff the given matrix is already a signed version of its support matrix.
  /// Running time: O(height * width * min (height, width) )

  bool is_signed_matrix (const integer_matrix& matrix);

  /// Returns true, iff the given matrix is already a signed version of its support matrix.
  /// If not, violator describes a violating submatrix.
  /// Running time: O(height * width * min (height, width) )

  bool is_signed_matrix (const integer_matrix& matrix, submatrix_indices& violator);

  /// Returns true, iff the given matrix is already a signed version of its support matrix.
  /// If not, the given matrix is changed to be such a signed version.
  /// Running time: O(height * width * min (height, width) )

  bool sign_matrix (integer_matrix& matrix);

  /// Drops all signs in the given matrix.

  void support_matrix (integer_matrix& matrix);

  /// Calculates a subdeterminant of the given matrix.

  int determinant_submatrix (const integer_matrix& matrix, const submatrix_indices& submatrix);

  /// Returns true, iff the given matrix is totally unimodular by checking all subdeterminants.

  bool determinant_is_totally_unimodular (const integer_matrix& matrix);

  /// Returns true, iff the given matrix is totally unimodular by checking all subdeterminants.
  /// If this is not the case, violator describes a violating submatrix.

  bool determinant_is_totally_unimodular (const integer_matrix& matrix, submatrix_indices& violator);

  /// Returns true, iff the given matrix is totally unimodular via ghouila-houri's criterion by enumeration of row subsets.

  bool ghouila_houri_is_totally_unimodular_enum_rows (const integer_matrix& matrix);

  /// Returns true, iff the given matrix is totally unimodular via ghouila-houri's criterion by enumeration of column subsets.

  bool ghouila_houri_is_totally_unimodular_enum_columns (const integer_matrix& matrix);

  /// Returns true, iff the given matrix is totally unimodular via ghouila-houri's criterion.
  /// Calls either the row or column version, depending on the dimensions of the matrix.

  bool ghouila_houri_is_totally_unimodular (const integer_matrix& matrix);

}

#endif /* TOTAL_UNIMODULARITY_HPP_ */

