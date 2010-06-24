/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef TOTAL_UNIMODULARITY_HPP_
#define TOTAL_UNIMODULARITY_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace tu {

  /**
   * Represents subsets of row and column index sets
   */

  struct submatrix_indices
  {
    typedef boost::numeric::ublas::vector <size_t> vector_type;
    typedef boost::numeric::ublas::indirect_array <vector_type> indirect_array_type;

    indirect_array_type rows;
    indirect_array_type columns;
  };

  /**
   * Integer matrix
   */

  typedef boost::numeric::ublas::matrix <int> integer_matrix;

  /**
   * Indirect integer matrix
   */

  typedef boost::numeric::ublas::matrix_indirect <const integer_matrix, submatrix_indices::indirect_array_type> integer_submatrix;

  /// Node of a decomposition tree

  class decomposed_matroid;

  enum log_level
  {
    LOG_QUIET, LOG_UPDATING, LOG_VERBOSE
  };

  /**
   * Returns a decomposition of a given binary matroid into a k-sum-decomposition (k=1,2,3)
   * in graphic, cographic, R10 and maybe irregular components.
   *
   * @param matrix Representation matrix of a binary matroid
   * @param level Log level
   * @return Root of decomposition tree
   */

  decomposed_matroid* decompose_binary_matroid (const integer_matrix& matrix, log_level level = LOG_QUIET);

  /**
   * Tests for total unimodularity without certificates.
   *
   * @param matrix The matrix to be tested
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  bool is_totally_unimodular (const integer_matrix& matrix, log_level level = LOG_QUIET);

  /**
   * Tests for total unimodularity, returning a decomposition of a given binary matroid
   * into a k-sum-decomposition (k=1,2,3) in graphic, cographic, R10 and maybe irregular
   * components if the matrix is totally unimodular.
   *
   * @param matrix The matrix to be tested
   * @param decomposition Returns root of decomposition tree
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  bool is_totally_unimodular (const integer_matrix& matrix, decomposed_matroid*& decomposition, log_level level = LOG_QUIET);

  /**
   * Tests for total unimodularity, returning a decomposition of a given binary matroid
   * into a k-sum-decomposition (k=1,2,3) in graphic, cographic, R10 and maybe irregular
   * components if the matrix is totally unimodular and the indices of a violating submatrix
   * otherwise.
   *
   * @param matrix The matrix to be tested
   * @param decomposition Returns root of decomposition tree
   * @param violator Returns violator indices
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  bool is_totally_unimodular (const integer_matrix& matrix, submatrix_indices& violator, log_level level = LOG_QUIET);

  /**
   * Tests for total unimodularity, returning the indices of a violating submatrix
   * if the matrix is not totally unimodular.
   *
   * @param matrix The matrix to be tested
   * @param violator Returns violator indices
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  bool is_totally_unimodular (const integer_matrix& matrix, decomposed_matroid*& decomposition, submatrix_indices& violator, log_level level =
      LOG_QUIET);

  /**
   * Tests if the given matrix contains only -1,0,+1 entries.
   *
   * @param matrix The given matrix
   * @return true if and only if it is a -1,0,+1 matrix
   */

  bool is_zero_plus_minus_one_matrix (const integer_matrix& matrix);

  /**
   * Tests if the given matrix contains only -1,0,+1 entries,
   * returning the violating position if this is not the case.
   *
   * @param matrix The given matrix
   * @param position Returns a violating entry
   * @return true if and only if it is a -1,0,+1 matrix
   */

  bool is_zero_plus_minus_one_matrix (const integer_matrix& matrix, std::pair <size_t, size_t>& position);

  /**
   * Tests if the given matrix contains only 0 or 1 entries.
   *
   * @param matrix The given matrix
   * @return true if and only if it is a 0-1 matrix
   */

  bool is_zero_one_matrix (const integer_matrix& matrix);

  /**
   * Tests if the given matrix contains only 0 or 1 entries,
   * returning the violating position if this is not the case.
   *
   * @param matrix The given matrix
   * @param position Returns a violating entry
   * @return true if and only if it is a 0-1 matrix
   */

  bool is_zero_one_matrix (const integer_matrix& matrix, std::pair <size_t, size_t>& position);

  /**
   * Tests if a given matrix is a signed version of its support matrix already.
   * Running time: O(height * width * min(height, width))
   *
   * @param matrix The matrix to be tested
   * @return true if and only if the support matrix can be signed to the orignal
   */

  bool is_signed_matrix (const integer_matrix& matrix);

  /**
   * Tests if a given matrix is a signed version of its support matrix already,
   * returning the indices of a violating submatrix if this is not the case.
   * Running time: O(height * width * min(height, width))
   *
   * @param matrix The matrix to be tested
   * @param violator Returns violator indices
   * @return true if and only if the support matrix can be signed to the original
   */

  bool is_signed_matrix (const integer_matrix& matrix, submatrix_indices& violator);

  /**
   * Signes a given matrix to be a signed version of its support matrix.
   * Running time: O(height * width * min(height, width))
   *
   * @param matrix The matrix to be signed
   * @return true if and only if any change was necessary
   */

  bool sign_matrix (integer_matrix& matrix);

  /**
   * Makes the matrix its own support matrix.
   *
   * @param matrix The given matrix
   */

  void support_matrix (integer_matrix& matrix);

  /**
   * Calculates a subdeterminant of the given matrix.
   *
   * @param matrix A given integer matrix
   * @param submatrix Matrix-indices describing a submatrix
   * @return The submatrix' determinant
   */

  int determinant_submatrix (const integer_matrix& matrix, const submatrix_indices& submatrix);

  /**
   * Checks all subdeterminants to test a given matrix for total unimodularity.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimdular
   */

  bool determinant_is_totally_unimodular (const integer_matrix& matrix);

  /**
   * Checks all subdeterminants to test a given matrix for total unimodularity.
   * If this is not the case, violator describes a violating submatrix.
   *
   * @param matrix The given matrix
   * @param violator The violating submatrix, if the result is false
   * @return true if and only if the this matrix is totally unimodular
   */

  bool determinant_is_totally_unimodular (const integer_matrix& matrix, submatrix_indices& violator);

  /**
   * Tests a given matrix to be totally unimodular using ghouila-houri's criterion by enumeration of row subsets.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimodular
   */

  bool ghouila_houri_is_totally_unimodular_enum_rows (const integer_matrix& matrix);

  /**
   * Tests a given matrix to be totally unimodular using ghouila-houri's criterion by enumeration of column subsets.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimodular
   */

  bool ghouila_houri_is_totally_unimodular_enum_columns (const integer_matrix& matrix);

  /**
   * Tests a given matrix to be totally unimodular using ghouila-houri's criterion by enumeration of either
   * row or column subsets by choosing the one which induces fewer enumerations.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimodular
   */

  bool ghouila_houri_is_totally_unimodular (const integer_matrix& matrix);

}

#endif /* TOTAL_UNIMODULARITY_HPP_ */

