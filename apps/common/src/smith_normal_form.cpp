/*
 * smith_normal_form.cpp
 *
 *  Created on: Feb 3, 2011
 *      Author: xammy
 */

#include "smith_normal_form.hpp"

#include "matrix.hpp"
#include "matrix_permuted.hpp"

#include <assert.h>
#include <set>

#include <boost/numeric/ublas/io.hpp>

namespace unimod
{
  int gcd_impl(int a, int b, int& s, int& t)
  {
    assert(a >= 0 && b >= 0 && a >= b);

    if (b == 0)
    {
      s = 1;
      t = 0;
      return a;
    }

    int q = a / b;
    int r = a % b;
    int result = gcd_impl(b, r, t, s);
    t -= q * s;
    return result;
  }

  int gcd(int a, int b, int& s, int& t)
  {
    if (a >= 0 && b >= 0)
    {
      if (a >= b)
        return gcd_impl(a, b, s, t);
      else
        return gcd_impl(b, a, t, s);
    }
    else
    {
      bool a_neg = a < 0;
      bool b_neg = b < 0;
      a = a >= 0 ? a : -a;
      b = b >= 0 ? b : -b;
      int result;
      if (a >= b)
        result = gcd_impl(a, b, s, t);
      else
        result = gcd_impl(b, a, t, s);
      if (a_neg)
        s = -s;
      if (b_neg)
        t = -t;
      return result;
    }
  }

  int gcd(int a, int b)
  {
    int s, t;
    return gcd(a, b, s, t);
  }

  template <typename Matrix>
  bool find_smallest_nonzero_matrix_entry(Matrix matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond, size_t& row,
      size_t& column)
  {
    bool result = false;
    int current_value = 0;
    for (size_t r = row_first; r != row_beyond; ++r)
    {
      for (size_t c = column_first; c != column_beyond; ++c)
      {
        int value = matrix(r, c);
        if (value == 0)
          continue;

        value = value >= 0 ? value : -value;

        if (!result || value < current_value)
        {
          result = true;
          row = r;
          column = c;
          current_value = value;
        }
      }
    }
    return result;
  }

  template <typename Matrix>
  bool find_smallest_gcd_matrix_entry(Matrix matrix, size_t row_first, size_t row_beyond, size_t column_first, size_t column_beyond, int entry,
      size_t& row, size_t& column)
  {
    int best_gcd = entry;
    for (size_t r = row_first; r != row_beyond; ++r)
    {
      for (size_t c = column_first; c != column_beyond; ++c)
      {
        int value = matrix(r, c);
        if (value == 0)
          continue;

        int current_gcd = gcd(entry, value);

        if (current_gcd < best_gcd)
        {
          row = r;
          column = c;
          best_gcd = current_gcd;
        }
      }
    }
    return best_gcd < entry;
  }

  void smith_normal_form(const integer_matrix& input_matrix, std::vector <int>& diagonal)
  {
    integer_matrix matrix = input_matrix;
    matrix_permuted <integer_matrix> permuted_matrix(matrix);
    size_t handled = 0;
    size_t row = 0;
    size_t column = 0;
    diagonal.resize(matrix.size1() < matrix.size2() ? matrix.size1() : matrix.size2(), 0);

    while (find_smallest_nonzero_matrix_entry(permuted_matrix, handled, permuted_matrix.size1(), handled, permuted_matrix.size2(), row, column))
    {
      matrix_permute1(permuted_matrix, handled, row);
      matrix_permute2(permuted_matrix, handled, column);

      std::cout << "We now have smallest entry " << permuted_matrix(handled, handled) << " at " << handled << " x " << handled << std::endl;
      matrix_print(permuted_matrix);

      bool changed;
      do
      {
        changed = false;

        for (row = handled + 1; row < permuted_matrix.size1(); ++row)
        {

        }
      }
      while (changed);

      std::set <size_t> nonzero_rows;
      for (row = handled + 1; row < permuted_matrix.size1(); ++row)
      {
        if (permuted_matrix(row, handled))
          nonzero_rows.insert(row);
      }

      while (find_smallest_gcd_matrix_entry(permuted_matrix, handled + 1, permuted_matrix.size1(), handled + 1, permuted_matrix.size2(),
          permuted_matrix(handled, handled), row, column))
      {
        matrix_permute1(permuted_matrix, handled + 1, row);
        matrix_permute2(permuted_matrix, handled + 1, column);

        std::cout << "We can make the gcd smaller with " << permuted_matrix(handled + 1, handled + 1) << " at " << handled + 1 << " x " << handled
            + 1 << std::endl;
        matrix_print(permuted_matrix);
      }

      handled++;
    }

  }
}
