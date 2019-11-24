/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace unimod
{

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

  typedef boost::numeric::ublas::matrix <long long> integer_matrix;

  /**
   * Indirect integer matrix
   */

  typedef boost::numeric::ublas::matrix_indirect <const integer_matrix, submatrix_indices::indirect_array_type> integer_submatrix;

  enum log_level
  {
    LOG_QUIET, LOG_PROGRESSIVE, LOG_VERBOSE
  };

}
#endif /* COMMON_HPP_ */
