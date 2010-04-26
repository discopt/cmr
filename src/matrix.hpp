
//          Copyright Matthias Walter 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include "../config.h"
#include <boost/numeric/ublas/matrix.hpp>
#include "matrix_transposed.hpp"

namespace tu {

  class matrix_binary_pivot_exception: public std::exception
  {
  public:
    const char* what () const throw ()
    {
      return "Cannot pivot on a zero entry!";
    }

  };

  template <typename MatrixType>
  inline void matrix_set_value (MatrixType& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {
    //    std::cout << "Changing matrix value at " << row << "," << column << std::endl;

    matrix (row, column) = value;
  }

  template <typename MatrixType>
  inline void matrix_set_value (const MatrixType& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {
    // This routine should not be called, but must exist to compile.
    assert (false);
  }

  template <typename MatrixType>
  inline void matrix_permute1 (MatrixType& matrix, size_t index1, size_t index2)
  {
    for (size_t index = 0; index < matrix.size2 (); ++index)
    {
      std::swap (matrix (index1, index), matrix (index2, index));
    }
  }

  template <typename MatrixType>
  inline void matrix_permute2 (MatrixType& matrix, size_t index1, size_t index2)
  {
    for (size_t index = 0; index < matrix.size1 (); ++index)
    {
      std::swap (matrix (index, index1), matrix (index, index2));
    }
  }

  template <typename MatrixType>
  void matrix_binary_pivot (MatrixType& matrix, size_t i, size_t j)
  {
    typedef typename MatrixType::value_type value_type;
    const value_type& base_value = matrix (i, j);

    if (base_value == 0)
    {
      throw matrix_binary_pivot_exception ();
    }

    for (size_t row = 0; row < matrix.size1 (); ++row)
    {
      if (row == i)
      {
        continue;
      }
      const value_type& first = matrix (row, j);
      if (first == 0)
      {
        continue;
      }

      for (size_t column = 0; column < matrix.size2 (); ++column)
      {
        if (column == j)
        {
          continue;
        }
        const value_type& second = matrix (i, column);
        if (second == 0)
        {
          continue;
        }
        matrix (row, column) = 1 - matrix (row, column);
      }
    }
  }

  template <typename MatrixType, typename PropertyCheck>
  inline size_t matrix_count_property_row_series (const MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first,
      size_t column_beyond, PropertyCheck check)
  {
    for (size_t row = row_first; row < row_beyond; ++row)
    {
      for (size_t column = column_first; column < column_beyond; ++column)
        check (matrix (row, column));
      if (!check ())
        return row - row_first;
    }
    return row_beyond - row_first;
  }

  template <typename MatrixType, typename PropertyCheck>
  inline size_t matrix_count_property_column_series (const MatrixType& matrix, size_t row_first, size_t row_beyond, size_t column_first,
      size_t column_beyond, PropertyCheck check)
  {
    matrix_transposed <const MatrixType> transposed (matrix);
    return matrix_count_property_row_series (transposed, column_first, column_beyond, row_first, row_beyond, check);
  }

  template <typename MatrixType>
  inline void matrix_print (const MatrixType& matrix)
  {
    for (size_t row = 0; row < matrix.size1 (); ++row)
    {
      for (size_t column = 0; column < matrix.size2 (); ++column)
      {
        std::cout << " " << matrix (row, column);
      }
      std::cout << "\n";
    }
    std::cout << std::flush;
  }

  template <typename MatrixType1, typename MatrixType2>
  inline bool matrix_equals (const MatrixType1& matrix1, const MatrixType2& matrix2)
  {
    if (matrix1.size1 () != matrix2.size1 ())
      return false;
    if (matrix1.size2 () != matrix2.size2 ())
      return false;

    for (size_t row = 0; row < matrix1.size1 (); ++row)
    {
      for (size_t column = 0; column < matrix1.size2 (); ++column)
      {
        if (matrix1 (row, column) != matrix2 (row, column))
          return false;
      }
    }
    return true;
  }

}

#endif /* MATRIX_HPP_ */
