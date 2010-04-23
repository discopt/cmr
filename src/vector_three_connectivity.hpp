/*
 * vector_three_connectivity.hpp
 *
 *  Created on: Jan 15, 2010
 *      Author: xammy
 */

#ifndef VECTOR_THREE_CONNECTIVITY_HPP_
#define VECTOR_THREE_CONNECTIVITY_HPP_

#include "../config.h"
#include <vector>

namespace tu {

  template <typename MatrixType>
  class vector_three_connectivity
  {
  public:
    typedef MatrixType matrix_type;
    typedef std::pair <unsigned char, size_t> vector_data;

    static const unsigned char ZERO_VECTOR = 0;
    static const unsigned char UNIT_VECTOR = 1;
    static const unsigned char PARALLEL = 2;
    static const unsigned char OTHER = 3;

    vector_three_connectivity (const matrix_type& matrix, size_t dimension, size_t base) :
      matrix_ (matrix), data_ (matrix.size2 ())
    {
      reset (dimension, base);
    }

    vector_three_connectivity (const vector_three_connectivity <MatrixType>& other) :
      matrix_ (other.matrix_)
    {
      dimension_ = other.dimension_;
      base_ = other.base_;
      data_ = other.data_;
    }

    void reset (size_t dimension, size_t base)
    {
      dimension_ = dimension;
      base_ = base;

      for (size_t column = 0; column < base_; ++column)
      {
        data_[column] = std::make_pair (PARALLEL, column);
      }
      for (size_t column = base_; column < matrix_.size2 (); ++column)
      {
        size_t count = 0;
        size_t last_one = 0;
        for (size_t r = 0; r < dimension_; ++r)
        {
          if (matrix_ (r, column) != 0)
          {
            last_one = r;
            ++count;
          }
        }
        if (count == 0)
          data_[column] = std::make_pair (ZERO_VECTOR, 0);
        else if (count == 1)
          data_[column] = std::make_pair (UNIT_VECTOR, last_one);
        else
        {
          data_[column] = std::make_pair (OTHER, 0);
          for (size_t b = 0; b < base_; ++b)
          {
            if (are_equal (0, dimension_, column, b))
            {
              data_[column] = std::make_pair (PARALLEL, b);
            }
          }
        }
      }
    }

    ~vector_three_connectivity ()
    {

    }

    inline size_t base () const
    {
      return base_;
    }

    inline size_t dimension () const
    {
      return dimension_;
    }

    bool is_parallel (size_t index) const
    {
      return data_[index].first == PARALLEL;
    }

    bool is_zero (size_t index) const
    {
      return data_[index].first == ZERO_VECTOR;
    }

    bool is_unit (size_t index) const
    {
      return data_[index].first == UNIT_VECTOR;
    }

    bool is_other (size_t index) const
    {
      return data_[index].first == OTHER;
    }

    size_t get_referred (size_t index) const
    {
      return data_[index].second;
    }

    void swap_cross (size_t index1, size_t index2)
    {
      assert (index1 < dimension_);
      assert (index2 < dimension_);

      for (size_t column = base_; column < data_.size (); ++column)
      {
        if (is_unit (column))
        {
          if (get_referred (column) == index1)
            data_[column].second = index2;
          else if (get_referred (column) == index2)
            data_[column].second = index1;
        }
      }
    }

    void swap_vectors (size_t index1, size_t index2)
    {
      assert (index1 >= base_);
      assert (index2 >= base_);

      std::swap (data_[index1], data_[index2]);
    }

    void enlarge_base (size_t amount = 1)
    {
      while (amount--)
      {
        data_[base_] = std::make_pair (PARALLEL, base_);
        for (size_t column = base_ + 1; column < matrix_.size2 (); ++column)
        {
          if (are_equal (0, dimension_, column, base_))
          {
            data_[column] = std::make_pair (PARALLEL, base_);
          }
        }
        ++base_;
      }
    }

    void enlarge_dimension (size_t amount = 1)
    {
      //        std::cout << "Before:" << std::endl;
      //        for (size_t c = 0; c < 10; ++c)
      //        {
      //            std::cout << "Vector " << c << "(";
      //            for (size_t r = 0; r < dimension_; ++r)
      //            {
      //                std::cout << " " << matrix_ (r, c);
      //            }
      //            std::cout << " ) has data = " << data_[c] << std::endl;
      //        }

      for (size_t column = base_; column < matrix_.size2 (); ++column)
      {
        if (data_[column].first == OTHER)
          continue;

        size_t count = 0;
        size_t last_row = 0;
        for (size_t row = 0; row < amount; ++row)
        {
          if (matrix_ (dimension_ + row, column) != 0)
          {
            ++count;
            last_row = dimension_ + row;
          }
        }

        if (data_[column].first == ZERO_VECTOR && count == 1)
        {
          data_[column] = std::make_pair (UNIT_VECTOR, last_row);
        }
        else if ((data_[column].first == ZERO_VECTOR && count > 1) || (data_[column].first == UNIT_VECTOR && count > 0))
        {
          data_[column].first = OTHER;
          for (size_t b = 0; b < base_; ++b)
          {
            if (are_equal (0, dimension_ + amount, column, b))
            {
              data_[column] = std::make_pair (PARALLEL, b);
            }
          }
        }
        else if (data_[column].first == PARALLEL)
        {
          if (!are_equal (dimension_, dimension_ + amount, column, data_[column].second))
          {
            data_[column].first = OTHER;
          }
        }
      }

      dimension_ += amount;

      //        std::cout << "After:" << std::endl;
      //        for (size_t c = 0; c < 10; ++c)
      //        {
      //            std::cout << "Vector " << c << "(";
      //            for (size_t r = 0; r < dimension_; ++r)
      //            {
      //                std::cout << " " << matrix_ (r, c);
      //            }
      //            std::cout << " ) has data = " << data_[c] << std::endl;
      //        }
    }

    bool operator== (const vector_three_connectivity <MatrixType>& other)
    {
      if (other.base_ != base_)
        return false;
      if (other.dimension_ != dimension_)
        return false;
      for (size_t i = 0; i < data_.size (); ++i)
      {
        if (data_[i] != other.data_[i])
          return false;
      }
      return true;
    }

    std::ostream& print (std::ostream& stream)
    {
      stream << "base = " << base_ << ", dim = " << dimension_ << "\n";
      for (size_t i = 0; i < data_.size (); ++i)
      {
        stream << "data[" << i << "] = " << ((int) data_[i].first) << ", " << data_[i].second << std::endl;
      }
      return stream;
    }

  private:

    bool are_equal (size_t row_first, size_t row_beyond, size_t column1, size_t column2)
    {
      for (size_t row = row_first; row < row_beyond; ++row)
      {
        if (matrix_ (row, column1) != matrix_ (row, column2))
          return false;
      }
      return true;
    }

    const matrix_type& matrix_;
    size_t dimension_;
    size_t base_;
    std::vector <vector_data> data_;
  };

}

#endif /* VECTOR_THREE_CONNECTIVITY_HPP_ */
