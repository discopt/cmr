/*
 * determinant.cpp
 *
 *  Created on: Oct 14, 2009
 *      Author: xammy
 */

#include "total_unimodularity.hpp"
#include "../config.h"

#include <boost/numeric/ublas/lu.hpp>
#include <boost/dynamic_bitset.hpp>

#include <map>
#include <vector>

namespace tu {

  /// Calculates a subdeterminant of the given matrix.

  int determinant_submatrix (const boost::numeric::ublas::matrix <int>& input_matrix, const submatrix_indices& submatrix)
  {
    const boost::numeric::ublas::matrix_indirect <const boost::numeric::ublas::matrix <int>, submatrix_indices::indirect_array_type> indirect_matrix (
        input_matrix, submatrix.rows, submatrix.columns);

    boost::numeric::ublas::matrix <float> matrix (indirect_matrix);
    boost::numeric::ublas::permutation_matrix <size_t> permutation_matrix (matrix.size1 ());

    int result = boost::numeric::ublas::lu_factorize (matrix, permutation_matrix);
    if (result != 0)
    {
      return 0;
    }

    float det = 1.0f;
    for (size_t i = 0; i < matrix.size1 (); ++i)
    {
      det *= matrix (i, i);
      if (i != permutation_matrix (i))
        det = -det;
    }

    return det;
  }

  /// Returns true, iff the given matrix is totally unimodular by checking all subdeterminants.

  bool determinant_is_totally_unimodular (const boost::numeric::ublas::matrix <int>& matrix)
  {
    submatrix_indices indices;

    return determinant_is_totally_unimodular (matrix, indices);
  }

  /// Returns true, iff the given matrix is totally unimodular by checking all subdeterminants.
  /// If this is not the case, violator describes a violating submatrix.

  bool determinant_is_totally_unimodular (const boost::numeric::ublas::matrix <int>& matrix, submatrix_indices& violator)
  {
    typedef unsigned long long int bitset_type;

    bitset_type row_max = 1 << matrix.size1 ();
    bitset_type column_max = 1 << matrix.size2 ();

    // collect all choices for column range with the cardinality as key
    std::map <size_t, std::vector <bitset_type> > column_bitsets;
    for (bitset_type choice = 1; choice < column_max; ++choice)
    {
      size_t cardinality = 0;
      for (size_t i = 0; i < matrix.size2 (); ++i)
      {
        if ((choice & (1 << i)) != 0)
          cardinality++;
      }
      column_bitsets[cardinality].push_back (choice);
    }

    for (bitset_type row_choice = 1; row_choice < row_max; ++row_choice)
    {
      size_t cardinality = 0;
      for (size_t i = 0; i < matrix.size2 (); ++i)
      {
        if ((row_choice & (1 << i)) != 0)
          cardinality++;
      }
      const std::vector <bitset_type>& column_choices = column_bitsets[cardinality];
      for (std::vector <bitset_type>::const_iterator iter = column_choices.begin (); iter != column_choices.end (); ++iter)
      {
        const bitset_type column_choice = *iter;

        submatrix_indices sub;
        submatrix_indices::vector_type indirect_array (cardinality);
        size_t current = 0;
        for (size_t i = 0; i < matrix.size1 (); ++i)
        {
          if ((row_choice & (1 << i)) != 0)
            indirect_array[current++] = i;
        }
        sub.rows = submatrix_indices::indirect_array_type (cardinality, indirect_array);
        current = 0;
        for (size_t i = 0; i < matrix.size2 (); ++i)
        {
          if ((column_choice & (1 << i)) != 0)
            indirect_array[current++] = i;
        }
        sub.columns = submatrix_indices::indirect_array_type (cardinality, indirect_array);

        int det = determinant_submatrix (matrix, sub);
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
