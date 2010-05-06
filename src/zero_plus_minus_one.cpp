/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#include <utility>
#include <boost/numeric/ublas/matrix.hpp>

#include "../config.h"

namespace tu {

  /// Returns true, iff the given matrix contains only values {-1, 0, 1}.
  /// If not, position will contain the position of the wrong entry. 

  bool is_zero_plus_minus_one_matrix (const boost::numeric::ublas::matrix <int>& matrix, std::pair <size_t, size_t>& position)
  {
    for (size_t row = 0; row < matrix.size1(); ++row)
    {
      for (size_t column = 0; column < matrix.size2(); ++column)
      {
        const int value = matrix(row, column);
        if (value < -1 || value > 1)
        {
          position.first = row;
          position.second = column;
          return false;
        }
      }
    }

    return true;
  }

  /// Returns true, iff the given matrix contains only values {-1, 0, 1}. 

  bool is_zero_plus_minus_one_matrix (const boost::numeric::ublas::matrix <int>& matrix)
  {
    std::pair <size_t, size_t> result;

    return is_zero_plus_minus_one_matrix(matrix, result);
  }

  /// Returns true, iff the given matrix contains only values {0, 1}.
  /// If not, position will contain the position of the wrong entry. 

  bool is_zero_one_matrix (const boost::numeric::ublas::matrix <int>& matrix, std::pair <size_t, size_t>& position)
  {
    for (size_t i = 0; i < matrix.size1(); ++i)
    {
      for (size_t j = 0; j < matrix.size2(); ++j)
      {
        const int value = matrix(i, j);
        if (value < 0 || value > 1)
        {
          position .first = i;
          position.second = j;
          return false;
        }
      }
    }

    return true;
  }

  /// Returns true, iff the given matrix contains only values {0, 1}. 

  bool is_zero_one_matrix (const boost::numeric::ublas::matrix <int>& matrix)
  {
    std::pair <size_t, size_t> result;

    return is_zero_one_matrix(matrix, result);
  }

}

