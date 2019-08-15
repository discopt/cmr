/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef UNIMODULARITY_HPP_
#define UNIMODULARITY_HPP_

#include "common.hpp"

namespace unimod
{
  /**
   * Tests for unimodularity without certificates.
   * A matrix of rank r is unimodular if and only if for every
   * column basis B the gcd (greatest common divisor) of the
   * determinants of regular r x r submatrices of B is 1.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank k
   * @param level Log level
   * @return true if and only if the matrix is unimodular
   */

  bool is_unimodular(const integer_matrix& matrix, size_t& rank, log_level level = LOG_QUIET);

  /**
   * Tests for strong unimodularity without certificates.
   * A matrix is strongly unimodular if and only if
   * it and its transpose are unimodular.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank r of the matrix
   * @param level Log level
   * @return true if and only if the matrix is strongly unimodular
   */

  bool is_strongly_unimodular(const integer_matrix& matrix, size_t& rank, log_level level = LOG_QUIET);

  /**
   * Tests for k-modularity without certificates.
   * A matrix of rank r is k-modular if and only if for every
   * column basis B the gcd (greatest common divisor) of the
   * determinants of regular r x r submatrices of B has the same
   * absolute value.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank r
   * @param level Log level
   * @return true if and only if the matrix is k-modular
   */

  bool is_k_modular(const integer_matrix& matrix, size_t& rank, log_level level = LOG_QUIET);

  /**
   * Tests for k-modularity without certificates.
   * A matrix of rank r is k-modular if and only if for every
   * column basis B the gcd (greatest common divisor) of the
   * determinants of regular r x r submatrices of B has the same
   * absolute value k.
   * It also computes k.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank r
   * @param k The common absolute value
   * @param level Log level
   * @return true if and only if the matrix is k-modular
   */

  bool is_k_modular(const integer_matrix& matrix, size_t& rank, unsigned int& k, log_level level = LOG_QUIET);

  /**
   * Tests for strong k-modularity without certificates.
   * A matrix of rank r is k-modular if and only if
   * it and its transpose are k-modular.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank r
   * @param level Log level
   * @return true if and only if the matrix is k-modular
   */

  bool is_strongly_k_modular(const integer_matrix& matrix, size_t& rank, log_level level = LOG_QUIET);

  /**
   * Tests for strong k-modularity without certificates.
   * A matrix of rank r is k-modular if and only if
   * it and its transpose are k-modular.
   * It also computes k.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank r
   * @param k The common absolute value
   * @param level Log level
   * @return true if and only if the matrix is k-modular
   */

  bool is_strongly_k_modular(const integer_matrix& matrix, size_t& rank, unsigned int& k, log_level level = LOG_QUIET);

  /**
   * In case a matrix A is k-modular, it may lead to q-integrality of the
   * polyhedron A*x = b, x >= 0 if B*x = b*q for some basis B of A. This method
   * finds a minimal integer q, assuming that A is k-modular.
   * If q = 1, the polyhedron is integral, if q = 2 it is half-integral, etc.
   *
   * @param matrix The matrix A, being k-modular.
   * @param rhs The rhs, given as a m x 1 matrix.
   * @return Minimal q as defined above
   */

  unsigned int get_k_modular_integrality(const integer_matrix& matrix, const integer_matrix& rhs);

  /**
   * In case a matrix A is k-modular, it may lead to integrality of the
   * polyhedron A*x = b, x >= 0 if B*x = b for some basis B of A.
   * This method tests of that property.
   *
   * @param matrix The matrix A, having the Dantzig property.
   * @param rhs The rhs, given as a m x 1 matrix.
   * @return true if and only if the polyhedron is integral
   */

  bool is_k_modular_integral(const integer_matrix& matrix, const integer_matrix& rhs);

  /**
   * Tests if a matrix A is complement totally unimodular.
   *
   * @param matrix The matrix A, having the Dantzig property.
   * @param level log_level
   * @return true if and only if the matrix is ctu.
   */

  bool is_complement_total_unimodular(const integer_matrix& matrix, log_level level = LOG_QUIET);
}

#endif /* UNIMODULARITY_HPP_ */
