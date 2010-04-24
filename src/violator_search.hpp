/**
 *
 * Copyright (c) 2010 Matthias Walter (xammy@xammy.homelinux.net)
 *
 * Authors: Matthias Walter
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 */

#ifndef VIOLATOR_SEARCH_HPP_
#define VIOLATOR_SEARCH_HPP_

#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "algorithm.hpp"
#include "matroid.hpp"

namespace tu {
  namespace detail {

    inline matroid_element_set find_smallest_irregular_minor (const decomposed_matroid* decomposition, bool collect_extra_elements = true)
    {
      if (decomposition->is_leaf ())
      {
        const decomposed_matroid_leaf* leaf = (decomposed_matroid_leaf*) decomposition;

        if (leaf->is_regular ())
          return matroid_element_set ();

        matroid_element_set result;
        std::copy (leaf->elements ().begin (), leaf->elements ().end (), std::inserter (result, result.end ()));
        if (collect_extra_elements)
          std::copy (leaf->extra_elements ().begin (), leaf->extra_elements ().end (), std::inserter (result, result.end ()));
        return result;
      }
      else
      {
        const decomposed_matroid_separator* separator = (decomposed_matroid_separator*) decomposition;

        matroid_element_set first_elements = find_smallest_irregular_minor (separator->first (), collect_extra_elements);
        matroid_element_set second_elements = find_smallest_irregular_minor (separator->second (), collect_extra_elements);
        if (first_elements.empty ())
          return second_elements;
        else if (second_elements.empty ())
          return first_elements;
        else
          return (first_elements.size () < second_elements.size ()) ? first_elements : second_elements;
      }
    }

    template <typename InputIterator, typename OutputIterator1, typename OutputIterator2>
    std::pair <OutputIterator1, OutputIterator2> split_elements (InputIterator first, InputIterator beyond, OutputIterator1 rows,
        OutputIterator2 columns)
    {
      for (; first != beyond; ++first)
      {
        if (*first > 0)
          *columns++ = *first;
        else
          *rows++ = *first;
      }
      return std::make_pair (rows, columns);
    }

    template <typename MatrixType>
    void create_indirect_matroid (const MatrixType& input_matrix, const matroid_element_set& row_elements,
        const matroid_element_set& column_elements, integer_matroid& sub_matroid, submatrix_indices& sub_indices)
    {
      sub_matroid.resize (row_elements.size (), column_elements.size ());
      submatrix_indices::vector_type row_vector (row_elements.size ());
      submatrix_indices::vector_type column_vector (column_elements.size ());

      size_t index = row_elements.size () - 1;
      for (matroid_element_set::const_iterator iter = row_elements.begin (); iter != row_elements.end (); ++iter)
      {
        sub_matroid.name1 (index) = *iter;
        row_vector[index] = -1 - *iter;
        --index;
      }

      index = 0;
      for (matroid_element_set::const_iterator iter = column_elements.begin (); iter != column_elements.end (); ++iter)
      {
        sub_matroid.name2 (index) = *iter;
        column_vector[index] = -1 + *iter;
        ++index;
      }

      sub_indices.rows = submatrix_indices::indirect_array_type (row_vector.size (), row_vector);
      sub_indices.columns = submatrix_indices::indirect_array_type (column_vector.size (), column_vector);
    }

    /**
     * Identifies a violator in a given matrix, which is known to be not totally unimodular. 
     * 
     * @param input_matrix The given matrix 
     * @param violator Output parameter to put the violating submatrix into.
     */

    inline void search_violator (const integer_matrix& input_matrix, const matroid_element_set& row_elements,
        const matroid_element_set& column_elements, submatrix_indices& violator)
    {
      integer_matroid sub_matroid;
      submatrix_indices sub_indices;

      create_indirect_matroid (input_matrix, row_elements, column_elements, sub_matroid, sub_indices);

      boost::numeric::ublas::matrix_indirect <const integer_matrix, submatrix_indices::indirect_array_type> sub_matrix (input_matrix,
          sub_indices.rows, sub_indices.columns);

      matroid_print (sub_matroid, sub_matrix);

      //      violator.rows = submatrix_indices::indirect_array_type (input_matrix.size1 ());
      //      violator.columns = submatrix_indices::indirect_array_type (input_matrix.size2 ());
      //
      //      for (size_t row = 0; row < violator.rows.size (); ++row)
      //      {
      //        violator.rows[row] = row;
      //      }
      //
      //      for (size_t column = 0; column < violator.columns.size (); ++column)
      //      {
      //        violator .columns[column] = column;
      //      }
      //
      //      shrink_rows (input_matrix, violator.rows, violator.columns);
      //
      //      shrink_rows (view_matrix_transposed (input_matrix), violator.columns, violator.rows);

      throw std::runtime_error ("Search for violator is not yet implemented.");
    }

  }
}

#endif /* VIOLATOR_SEARCH_HPP_ */
