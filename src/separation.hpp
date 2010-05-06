/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef SEPARATION_HPP_
#define SEPARATION_HPP_

#include "../config.h"
#include <utility>
#include <cstddef>
#include <vector>
#include "nested_minor_sequence.hpp"
#include "matroid.hpp"
#include "total_unimodularity.hpp"

namespace tu {

  class separation
  {
  public:
    typedef std::pair <std::size_t, std::size_t> split_type;
    typedef std::pair <std::size_t, std::size_t> witness_type;

    separation ();
    separation (split_type split);
    separation (split_type split, witness_type witness1);
    separation (split_type split, witness_type witness1, witness_type witness2);
    separation (const separation& other);
    separation& operator= (const separation& other);

    ~separation ();

    separation transposed ();

    int upper_right_rank () const
    {
      return upper_right_rank_;
    }

    int lower_left_rank () const
    {
      return lower_left_rank_;
    }

    int rank () const
    {
      return lower_left_rank() + upper_right_rank();
    }

    const std::pair <size_t, size_t>& split () const
    {
      return split_;
    }

    bool is_valid () const
    {
      return split_.first > 0 || split_.second > 0;
    }

    const witness_type& witness (size_t id = 0) const
    {
      return witnesses_[id];
    }

    void set_special_swap (char type, size_t index)
    {
      special_swap_ = std::make_pair(type, index);
    }

    bool has_special_swap () const
    {
      return special_swap_.first != 0;
    }

    bool has_special_row_swap () const
    {
      return special_swap_.first == 'r';
    }

    size_t get_special_swap_index () const
    {
      return special_swap_.second;
    }

    template <typename MatroidType, typename MatrixType>
    void create_components (const MatroidType& matroid, const MatrixType& matrix, integer_matroid& upper_left_matroid,
        integer_matrix& upper_left_matrix, integer_matroid& lower_right_matroid, integer_matrix& lower_right_matrix) const
    {
      assert (is_valid());
      assert (matroid.size1() == matrix.size1());
      assert (matroid.size2() == matrix.size2());

      if (rank() == 0)
      {
        upper_left_matroid.resize(split_.first, split_.second);
        for (size_t row = 0; row < upper_left_matroid.size1(); ++row)
          upper_left_matroid.name1(row) = matroid.name1(row);
        for (size_t column = 0; column < upper_left_matroid.size2(); ++column)
          upper_left_matroid.name2(column) = matroid.name2(column);
        upper_left_matrix.resize(split_.first, split_.second, false);
        for (size_t row = 0; row < upper_left_matrix.size1(); ++row)
        {
          for (size_t column = 0; column < upper_left_matrix.size2(); ++column)
            upper_left_matrix(row, column) = matrix(row, column);
        }

        lower_right_matroid.resize(matroid.size1() - split_.first, matroid.size2() - split_.second);
        for (size_t row = 0; row < lower_right_matroid.size1(); ++row)
          lower_right_matroid.name1(row) = matroid.name1(split_.first + row);
        for (size_t column = 0; column < lower_right_matroid.size2(); ++column)
          lower_right_matroid.name2(column) = matroid.name2(split_.second + column);
        lower_right_matrix.resize(matrix.size1() - split_.first, matrix.size2() - split_.second, false);
        for (size_t row = 0; row < lower_right_matrix.size1(); ++row)
        {
          for (size_t column = 0; column < lower_right_matrix.size2(); ++column)
            lower_right_matrix(row, column) = matrix(split_.first + row, split_.second + column);
        }
      }
      else if (rank() == 1 && witness().first >= split_.first)
      {
        upper_left_matroid.resize(split_.first + 1, split_.second);
        for (size_t row = 0; row < split_.first; ++row)
          upper_left_matroid.name1(row) = matroid.name1(row);
        upper_left_matroid.name1(split_.first) = matroid.name1(witness().first);
        for (size_t column = 0; column < upper_left_matroid.size2(); ++column)
          upper_left_matroid.name2(column) = matroid.name2(column);

        upper_left_matrix.resize(split_.first + 1, split_.second, false);
        for (size_t column = 0; column < upper_left_matrix.size2(); ++column)
        {
          for (size_t row = 0; row < split_.first; ++row)
            upper_left_matrix(row, column) = matrix(row, column);
          upper_left_matrix(split_.first, column) = matrix(witness().first, column);
        }

        lower_right_matroid.resize(matroid.size1() - split_.first, matroid.size2() - split_.second + 1);
        for (size_t row = 0; row < lower_right_matroid.size1(); ++row)
          lower_right_matroid.name1(row) = matroid.name1(split_.first + row);
        for (size_t column = 1; column < lower_right_matroid.size2(); ++column)
          lower_right_matroid.name2(column) = matroid.name2(split_.second + column - 1);
        lower_right_matroid.name2(0) = matroid.name2(witness().second);

        lower_right_matrix.resize(matrix.size1() - split_.first, matrix.size2() - split_.second + 1, false);
        for (size_t row = 0; row < lower_right_matrix.size1(); ++row)
        {
          for (size_t column = 1; column < lower_right_matrix.size2(); ++column)
            lower_right_matrix(row, column) = matrix(split_.first + row, split_.second + column - 1);
          lower_right_matrix(row, 0) = matrix(split_.first + row, witness().second);
        }
      }
      else if (rank() == 1 && witness().second >= split_.second)
      {
        upper_left_matroid.resize(split_.first, split_.second + 1);
        for (size_t row = 0; row < split_.first; ++row)
          upper_left_matroid.name1(row) = matroid.name1(row);
        for (size_t column = 0; column < split_.second; ++column)
          upper_left_matroid.name2(column) = matroid.name2(column);
        upper_left_matroid.name2(split_.second) = matroid.name2(witness().second);

        upper_left_matrix.resize(split_.first, split_.second + 1, false);
        for (size_t row = 0; row < split_.first; ++row)
        {
          for (size_t column = 0; column < split_.second; ++column)
            upper_left_matrix(row, column) = matrix(row, column);
          upper_left_matrix(row, split_.second) = matrix(row, witness().second);
        }

        lower_right_matroid.resize(matroid.size1() - split_.first + 1, matroid.size2() - split_.second);
        for (size_t row = 1; row < lower_right_matroid.size1(); ++row)
          lower_right_matroid.name1(row) = matroid.name1(split_.first + row - 1);
        for (size_t column = 0; column < lower_right_matroid.size2(); ++column)
          lower_right_matroid.name2(column) = matroid.name2(split_.second + column);
        lower_right_matroid.name1(0) = matroid.name1(witness().first);

        lower_right_matrix.resize(matrix.size1() - split_.first + 1, matrix.size2() - split_.second, false);
        for (size_t column = 0; column < lower_right_matrix.size2(); ++column)
        {
          for (size_t row = 1; row < lower_right_matrix.size1(); ++row)
            lower_right_matrix(row, column) = matrix(split_.first + row - 1, split_.second + column);
          lower_right_matrix(0, column) = matrix(witness().first, split_.second + column);
        }
      }
      else
      {
        upper_left_matroid.resize(split_.first + 2, split_.second + 1);
        for (size_t row = 0; row < split_.first + 2; ++row)
          upper_left_matroid.name1(row) = matroid.name1(row);
        for (size_t column = 0; column < split_.second + 1; ++column)
          upper_left_matroid.name2(column) = matroid.name2(column);

        upper_left_matrix.resize(split_.first + 2, split_.second + 1, false);
        for (size_t row = 0; row < split_.first + 2; ++row)
        {
          for (size_t column = 0; column < split_.second + 1; ++column)
            upper_left_matrix(row, column) = matrix(row, column);
        }

        lower_right_matroid.resize(matroid.size1() - split_.first + 1, matroid.size2() - split_.second + 2);
        for (size_t row = 0; row < lower_right_matroid.size1(); ++row)
          lower_right_matroid.name1(row) = matroid.name1(split_.first + row - 1);
        for (size_t column = 0; column < lower_right_matroid.size2(); ++column)
          lower_right_matroid.name2(column) = matroid.name2(split_.second + column - 2);

        lower_right_matrix.resize(matrix.size1() - split_.first + 1, matrix.size2() - split_.second + 2, false);
        for (size_t column = 0; column < lower_right_matrix.size2(); ++column)
        {
          for (size_t row = 0; row < lower_right_matrix.size1(); ++row)
            lower_right_matrix(row, column) = matrix(split_.first + row - 1, split_.second + column - 2);
        }
      }
    }

  private:
    std::pair <size_t, size_t> split_;
    std::vector <std::pair <size_t, size_t> > witnesses_;
    int upper_right_rank_;
    int lower_left_rank_;
    std::pair <char, size_t> special_swap_;
  };

}

#endif /* SEPARATION_HPP_ */
