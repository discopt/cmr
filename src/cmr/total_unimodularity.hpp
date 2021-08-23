#pragma once

#include <cmr/config.h>
#include <cmr/export.h>

#include "common.hpp"

namespace tu
{
  /// Node of a decomposition tree

  class decomposed_matroid;

  /**
   * Returns a decomposition of a given binary matroid into a k-sum-decomposition (k=1,2,3)
   * in graphic, cographic, R10 and maybe irregular components.
   *
   * @param matrix Representation matrix of a binary matroid
   * @param level Log level
   * @return Root of decomposition tree
   */

  CMR_EXPORT
  decomposed_matroid* decompose_binary_matroid(const integer_matrix& matrix, log_level level = LOG_QUIET);

  /**
   * Tests for total unimodularity without certificates.
   * A matrix is totally unimodular if and only if every square submatrix
   * has a determinant of -1, 0 or +1.
   *
   * @param matrix The matrix to be tested
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  CMR_EXPORT
  bool is_totally_unimodular(const integer_matrix& matrix, log_level level = LOG_QUIET);

  /**
   * Tests for total unimodularity with a positive certificate.
   * A matrix is totally unimodular if and only if every square submatrix
   * has a determinant of -1, 0 or +1.
   * If the matrix has this property, the routine returns a decomposition of
   * the underlying binary matroid into a k-sum-decomposition (k=1,2,3)
   * in graphic, cographic, R10 and maybe irregular components.
   *
   * @param matrix The matrix to be tested
   * @param decomposition Returns the root of the decomposition tree
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  CMR_EXPORT
  bool is_totally_unimodular(const integer_matrix& matrix, decomposed_matroid*& decomposition, log_level level = LOG_QUIET);

  /**
   * Tests for total unimodularity with a negative certificate.
   * A matrix is totally unimodular if and only if every square submatrix
   * has a determinant of -1, 0 or +1.
   * If the matrix does not have this property,the routine returns the indices
   * of a violating submatrix.
   *
   * @param matrix The matrix to be tested
   * @param violator Returns violator indices
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  CMR_EXPORT
  bool is_totally_unimodular(const integer_matrix& matrix, submatrix_indices& violator, log_level level = LOG_QUIET);

  /**
   * Tests for total unimodularity with certificates.
   * A matrix is totally unimodular if and only if every square submatrix
   * has a determinant of -1, 0 or +1.
   * If the matrix has this property, the routine returns a decomposition of
   * the underlying binary matroid into a k-sum-decomposition (k=1,2,3)
   * in graphic, cographic, R10 and maybe irregular components.
   * If the matrix does not have this property,the routine returns the indices
   * of a violating submatrix.
   *
   * @param matrix The matrix to be tested
   * @param decomposition Returns the root of the decomposition tree.
   * @param violator Returns violator indices
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  CMR_EXPORT
  bool is_totally_unimodular(const integer_matrix& matrix, decomposed_matroid*& decomposition, submatrix_indices& violator, log_level level =
      LOG_QUIET);

  /**
   * Tests if the given matrix contains only -1,0,+1 entries.
   *
   * @param matrix The given matrix
   * @return true if and only if it is a -1,0,+1 matrix
   */

  CMR_EXPORT
  bool is_zero_plus_minus_one_matrix(const integer_matrix& matrix);

  /**
   * Tests if the given matrix contains only -1,0,+1 entries,
   * returning the violating position if this is not the case.
   *
   * @param matrix The given matrix
   * @param position Returns a violating entry
   * @return true if and only if it is a -1,0,+1 matrix
   */

  CMR_EXPORT
  bool is_zero_plus_minus_one_matrix(const integer_matrix& matrix, std::pair <integer_matrix::size_type, integer_matrix::size_type>& position);

  /**
   * Tests if the given matrix contains only 0 or 1 entries.
   *
   * @param matrix The given matrix
   * @return true if and only if it is a 0-1 matrix
   */

  CMR_EXPORT
  bool is_zero_one_matrix(const integer_matrix& matrix);

  /**
   * Tests if the given matrix contains only 0 or 1 entries,
   * returning the violating position if this is not the case.
   *
   * @param matrix The given matrix
   * @param position Returns a violating entry
   * @return true if and only if it is a 0-1 matrix
   */

  CMR_EXPORT
  bool is_zero_one_matrix(const integer_matrix& matrix, std::pair <integer_matrix::size_type, integer_matrix::size_type>& position);

  /**
   * Tests if a given matrix is a signed version of its support matrix already.
   * Running time: O(height * width * min(height, width))
   *
   * @param matrix The matrix to be tested
   * @return true if and only if the support matrix can be signed to the orignal
   */

  CMR_EXPORT
  bool is_signed_matrix(const integer_matrix& matrix);

  /**
   * Tests if a given matrix is a signed version of its support matrix already,
   * returning the indices of a violating submatrix if this is not the case.
   * Running time: O(height * width * min(height, width))
   *
   * @param matrix The matrix to be tested
   * @param violator Returns violator indices
   * @return true if and only if the support matrix can be signed to the original
   */

  CMR_EXPORT
  bool is_signed_matrix(const integer_matrix& matrix, submatrix_indices& violator);

  /**
   * Signes a given matrix to be a signed version of its support matrix.
   * Running time: O(height * width * min(height, width))
   *
   * @param matrix The matrix to be signed
   * @return true if and only if any change was necessary
   */

  CMR_EXPORT
  bool sign_matrix(integer_matrix& matrix);

  /**
   * Makes the matrix its own support matrix.
   *
   * @param matrix The given matrix
   */

  CMR_EXPORT
  void support_matrix(integer_matrix& matrix);

  /**
   * Calculates a subdeterminant of the given matrix.
   *
   * @param matrix A given integer matrix
   * @param submatrix Matrix-indices describing a submatrix
   * @return The submatrix' determinant
   */

  CMR_EXPORT
  int submatrix_determinant(const integer_matrix& matrix, const submatrix_indices& submatrix);

  /**
   * Checks all subdeterminants to test a given matrix for total unimodularity.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimdular
   */

  CMR_EXPORT
  bool determinant_is_totally_unimodular(const integer_matrix& matrix);

  /**
   * Checks all subdeterminants to test a given matrix for total unimodularity.
   * If this is not the case, violator describes a violating submatrix.
   *
   * @param matrix The given matrix
   * @param violator The violating submatrix, if the result is false
   * @return true if and only if the this matrix is totally unimodular
   */

  CMR_EXPORT
  bool determinant_is_totally_unimodular(const integer_matrix& matrix, submatrix_indices& violator);

  /**
   * Tests a given matrix to be totally unimodular using ghouila-houri's criterion by enumeration of row subsets.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimodular
   */

  CMR_EXPORT
  bool ghouila_houri_is_totally_unimodular_enum_rows(const integer_matrix& matrix);

  /**
   * Tests a given matrix to be totally unimodular using ghouila-houri's criterion by enumeration of column subsets.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimodular
   */

  CMR_EXPORT
  bool ghouila_houri_is_totally_unimodular_enum_columns(const integer_matrix& matrix);

  /**
   * Tests a given matrix to be totally unimodular using ghouila-houri's criterion by enumeration of either
   * row or column subsets by choosing the one which induces fewer enumerations.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimodular
   */

  CMR_EXPORT
  bool ghouila_houri_is_totally_unimodular(const integer_matrix& matrix);

} /* namespace tu */
