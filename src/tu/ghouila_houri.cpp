/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#include <vector>

#include <tu/total_unimodularity.hpp>

namespace tu
{

  /**
   * Class to enumerate subset-2-partitions for a ghouila-houri test
   */

  class ghouila_houri_enumerator
  {
  public:

    typedef signed char choice_type;
    typedef std::vector <choice_type> choice_vector_type;

    /**
     * Constructs an enumerator.
     *
     * @param matrix The matrix to be tested
     */

    ghouila_houri_enumerator(const integer_matrix& matrix) :
      _matrix(matrix)
    {
      _choice.resize(_matrix.size1());
      for (choice_vector_type::iterator iter = _choice.begin(); iter != _choice.end(); ++iter)
      {
        *iter = 0;
      }
    }

    /**
     * Tests a sum induced by a 2-partition of the selected rows
     *
     * @return true if and only if the sum along each column is valid
     */

    bool check_sum()
    {
      for (size_t column = 0; column < _matrix.size2(); ++column)
      {
        int sum = 0;
        for (size_t row = 0; row < _matrix.size1(); ++row)
        {
          sum += _choice[row] * _matrix(row, column);
        }
        if (sum < -1 || sum > +1)
          return false;
      }
      return true;
    }

    /**
     * Recursively enumerates both possible choices for a given row.
     *
     * @param row The given row
     * @return true if and only if one of the enumerations succeeded
     */

    bool choose_partition(size_t row = 0)
    {
      while (row < _choice.size())
      {
        if (_choice[row])
        {
          _choice[row] = -1;
          if (choose_partition(row + 1))
            return true;
          _choice[row] = 1;
          if (choose_partition(row + 1))
            return true;
          return false;
        }
        else
          ++row;
      }
      return check_sum();
    }

    /**
     * Recursively enumerates all subset-choices for a given row.
     * Finally it starts the enumeration of all the partitions of this subset.
     *
     * @param row The given row
     * @return true if and only if the further enumerations succeeded all.
     */

    bool choose_subset(size_t row = 0)
    {
      if (row < _choice.size())
      {
        _choice[row] = 0;
        if (!choose_subset(row + 1))
          return false;
        _choice[row] = 1;
        if (!choose_subset(row + 1))
          return false;
        return true;
      }

      return choose_partition();
    }

    /**
     * Tests the given matrix for total unimodularity using rows.
     *
     * @return true if and only if it is totally unimodular
     */

    inline bool check()
    {
      return choose_subset();
    }

  private:

    choice_vector_type _choice;
    const integer_matrix& _matrix;
  };

  /**
   * Tests a given matrix to be totally unimodular using ghouila-houri's criterion by enumeration of row subsets.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimodular
   */

  bool ghouila_houri_is_totally_unimodular_enum_rows(const integer_matrix& matrix)
  {
    ghouila_houri_enumerator enumerator(matrix);
    return enumerator.check();
  }

  /**
   * Tests a given matrix to be totally unimodular using ghouila-houri's criterion by enumeration of column subsets.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimodular
   */

  bool ghouila_houri_is_totally_unimodular_enum_columns(const integer_matrix& matrix)
  {
    const boost::numeric::ublas::matrix <int> transposed = boost::numeric::ublas::trans(matrix);

    return ghouila_houri_is_totally_unimodular_enum_rows(transposed);
  }

  /**
   * Tests a given matrix to be totally unimodular using ghouila-houri's criterion by enumeration of either
   * row or column subsets by choosing the one which induces fewer enumerations.
   *
   * @param matrix The given matrix
   * @return true if and only if this matrix is totally unimodular
   */

  bool ghouila_houri_is_totally_unimodular(const integer_matrix& matrix)
  {
    if (matrix.size1() > matrix.size2())
    {
      return ghouila_houri_is_totally_unimodular_enum_columns(matrix);
    }
    else
    {
      return ghouila_houri_is_totally_unimodular_enum_rows(matrix);
    }
  }

}
