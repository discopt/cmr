/*
 * ghouila_houri.cpp
 *
 *  Created on: Oct 18, 2009
 *      Author: xammy
 */

#include <boost/numeric/ublas/matrix.hpp>

#include "../config.h"

namespace tu {

  class ghouila_houri_enumerator
  {
  private:
    typedef signed char choice_type;
    typedef std::vector <choice_type> choice_vector_type;
    typedef boost::numeric::ublas::matrix <int> matrix_type;

    choice_vector_type _choice;
    const matrix_type& _matrix;
  public:
    ghouila_houri_enumerator (const matrix_type& matrix) :
      _matrix (matrix)
    {
      _choice.resize (_matrix.size1 ());
      for (choice_vector_type::iterator iter = _choice.begin (); iter != _choice.end (); ++iter)
      {
        *iter = 0;
      }
    }

    bool check_sum ()
    {
      for (size_t column = 0; column < _matrix.size2 (); ++column)
      {
        int sum = 0;
        for (size_t row = 0; row < _matrix.size1 (); ++row)
        {
          sum += _choice[row] * _matrix (row, column);
        }
        if (sum < -1 || sum > +1)
          return false;
      }
      return true;
    }

    bool choose_partition (size_t row = 0)
    {
      while (row < _choice.size ())
      {
        if (_choice[row])
        {
          _choice[row] = -1;
          if (choose_partition (row + 1))
            return true;
          _choice[row] = 1;
          if (choose_partition (row + 1))
            return true;
          return false;
        }
        else
          ++row;
      }
      return check_sum ();
    }

    bool choose_subset (size_t row = 0)
    {
      if (row < _choice.size ())
      {
        _choice[row] = 0;
        if (!choose_subset (row + 1))
          return false;
        _choice[row] = 1;
        if (!choose_subset (row + 1))
          return false;
        return true;
      }

      return choose_partition ();
    }

    inline bool check ()
    {
      return choose_subset ();
    }
  };

  /// Returns true, iff the given matrix is totally unimodular via ghouila-houri's criterion by enumeration of row subsets.

  bool ghouila_houri_is_totally_unimodular_enum_rows (const boost::numeric::ublas::matrix <int>& matrix)
  {
    ghouila_houri_enumerator enumerator (matrix);
    return enumerator.check ();
  }

  /// Returns true, iff the given matrix is totally unimodular via ghouila-houri's criterion by enumeration of column subsets.

  bool ghouila_houri_is_totally_unimodular_enum_columns (const boost::numeric::ublas::matrix <int>& matrix)
  {
    const boost::numeric::ublas::matrix <int> transposed = boost::numeric::ublas::trans (matrix);

    return ghouila_houri_is_totally_unimodular_enum_rows (transposed);
  }

  /// Returns true, iff the given matrix is totally unimodular via ghouila-houri's criterion.
  /// Calls either the row or column version, depending on the dimensions of the matrix.

  bool ghouila_houri_is_totally_unimodular (const boost::numeric::ublas::matrix <int>& matrix)
  {
    if (matrix.size1 () > matrix.size2 ())
    {
      return ghouila_houri_is_totally_unimodular_enum_columns (matrix);
    }
    else
    {
      return ghouila_houri_is_totally_unimodular_enum_rows (matrix);
    }
  }

}
