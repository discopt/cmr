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
   * A matrix of rank k is unimodular if and only if for every
   * column basis B the gcd (greatest common divisor) of the
   * determinants of regular k x k submatrices of B is 1.
   *
   * @param matrix The matrix to be tested
   * @param level Log level
   * @return true if and only if the matrix is unimodular
   */

  bool is_unimodular(const integer_matrix& matrix, log_level level = LOG_QUIET);

}

#endif /* UNIMODULARITY_HPP_ */
