/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#include "../config.h"
#include <map>
#include <vector>

#include <boost/numeric/ublas/lu.hpp>
#include <boost/dynamic_bitset.hpp>

#include "total_unimodularity.hpp"

namespace unimod
{

  /**
   * Calculates a subdeterminant of the given matrix.
   *
   * @param matrix A given integer matrix
   * @param submatrix Matrix-indices describing a submatrix
   * @return The submatrix' determinant
   */

  int submatrix_determinant(const integer_matrix& matrix, const submatrix_indices& submatrix)
  {
    typedef boost::numeric::ublas::matrix_indirect <const integer_matrix, submatrix_indices::indirect_array_type> indirect_matrix_t;

    assert (submatrix.rows.size() == submatrix.columns.size());

    const indirect_matrix_t indirect_matrix(matrix, submatrix.rows, submatrix.columns);

    boost::numeric::ublas::matrix <float> float_matrix(indirect_matrix);
    boost::numeric::ublas::permutation_matrix <size_t> permutation_matrix(float_matrix.size1());

    /// LU-factorize
    int result = boost::numeric::ublas::lu_factorize(float_matrix, permutation_matrix);
    if (result != 0)
      return 0;

    /// Calculate the determinant as a product of the diagonal entries, taking the permutation matrix into account.
    float det = 1.0f;
    for (size_t i = 0; i < float_matrix.size1(); ++i)
    {
      det *= float_matrix(i, i);
      if (i != permutation_matrix(i))
        det = -det;
    }

    return (int) det;
  }

  /**
   * Checks all subdeterminants to test a given matrix for total unimodularity.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimdular
   */

  bool determinant_is_totally_unimodular(const integer_matrix& matrix)
  {
    submatrix_indices indices;

    return determinant_is_totally_unimodular(matrix, indices);
  }

  /**
   * Checks all subdeterminants to test a given matrix for total unimodularity.
   * If this is not the case, violator describes a violating submatrix.
   *
   * @param matrix The given matrix
   * @param violator The violating submatrix, if the result is false
   * @return true if and only if the this matrix is totally unimodular
   */

  bool determinant_is_totally_unimodular(const integer_matrix& matrix, submatrix_indices& violator)
  {
    typedef unsigned long long int bitset_type;

    if (matrix.size1() > 8 * sizeof(bitset_type) || matrix.size2() > 8 * sizeof(bitset_type))
    {
      throw std::runtime_error("Cannot test such a large matrix for total unimodularity via determinants!");
    }

    bitset_type row_max = ((bitset_type) 1) << (bitset_type) matrix.size1();
    bitset_type column_max = ((bitset_type) 1) << (bitset_type) matrix.size2();

    /// Enumerate rows
    for (bitset_type row_choice = 1; row_choice < row_max; ++row_choice)
    {
      size_t row_cardinality = 0;
      for (size_t i = 0; i < matrix.size1(); ++i)
      {
        if ((row_choice & (((bitset_type) 1) << i)) != 0)
          row_cardinality++;
      }

      for (bitset_type column_choice = 1; column_choice < column_max; ++column_choice)
      {
        size_t column_cardinality = 0;
        for (size_t i = 0; i < matrix.size2(); ++i)
        {
          if ((column_choice & (((bitset_type) 1) << i)) != 0)
            column_cardinality++;
        }

        if (column_cardinality != row_cardinality)
          continue;

        submatrix_indices sub;
        submatrix_indices::vector_type indirect_array(row_cardinality);
        size_t current = 0;
        for (size_t i = 0; i < matrix.size1(); ++i)
        {
          if ((row_choice & (((bitset_type) 1) << i)) != 0)
            indirect_array[current++] = i;
        }
        sub.rows = submatrix_indices::indirect_array_type(row_cardinality, indirect_array);

        current = 0;
        for (size_t i = 0; i < matrix.size2(); ++i)
        {
          if ((column_choice & (((bitset_type) 1) << i)) != 0)
            indirect_array[current++] = i;
        }

        sub.columns = submatrix_indices::indirect_array_type(column_cardinality, indirect_array);

        /// Examine the determinant
        int det = submatrix_determinant(matrix, sub);
        if (det < -1 || det > 1)
        {
          violator = sub;
          return false;
        }
      }
    }

    return true;
  }

}
