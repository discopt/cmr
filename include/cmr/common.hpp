#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace tu
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

} /* namespace tu */
