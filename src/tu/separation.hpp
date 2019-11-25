#pragma once

#include <utility>
#include <cstddef>
#include <vector>
#include "nested_minor_sequence.hpp"
#include "matroid.hpp"
#include <tu/total_unimodularity.hpp>

namespace tu
{

  /**
   * Models different types of matroid separations.
   */

  class separation
  {
  public:
    typedef std::pair <std::size_t, std::size_t> split_type;
    typedef std::pair <std::size_t, std::size_t> witness_type;

    /**
     * Constructs a separation which is none.
     */

    separation();

    /**
     * Constructs a 1-separation.
     *
     * @param split Size of the first component, which must be at upper-left
     */

    separation(split_type split);

    /**
     * Constructs a 2-separation.
     *
     * @param split Size of the first component, which must be at upper-left
     * @param witness1 A witnessing one in the rank 1 part
     */

    separation(split_type split, witness_type witness1);

    /**
     * Constructs a 3-separation.
     *
     * @param split Size of the first component, which must be at upper-left
     * @param witness1 First witnessing one in the rank 2 part
     * @param witness2 Second witnessing one in the rank 2 part
     */

    separation(split_type split, witness_type witness1, witness_type witness2);

    /**
     * Copy constructor
     *
     * @param other Another separation
     */

    separation(const separation& other);

    /**
     * Assignment operator
     *
     * @param other Another separation
     * @return This separation
     */

    separation& operator=(const separation& other);

    /**
     * Destructor
     */

    ~separation();

    /**
     * @return A separation with the transposed details
     */

    separation transposed();

    /**
     * @return Rank in the upper right part
     */

    int upper_right_rank() const
    {
      return upper_right_rank_;
    }

    /**
     * @return Rank in the lower left part
     */

    int lower_left_rank() const
    {
      return lower_left_rank_;
    }

    /**
     * @return Total rank in lower left and upper right parts
     */

    int rank() const
    {
      return lower_left_rank() + upper_right_rank();
    }

    /**
     * @return The split
     */

    const std::pair <size_t, size_t>& split() const
    {
      return split_;
    }

    /**
     * @return true if and only if it is a real separation
     */

    bool is_valid() const
    {
      return split_.first > 0 || split_.second > 0;
    }

    /**
     * Returns a specific witness.
     *
     * @param index Index of the witness
     * @return The position of the witnessing one entry
     */

    const witness_type& witness(size_t index = 0) const
    {
      return witnesses_[index];
    }

    /**
     * Sets a special swap needed for 2-separations.
     *
     * @param type Type of the special swap. Can be 'r' and 'c'
     * @param index Row or column index of the swap
     */

    void set_special_swap(char type, size_t index)
    {
      special_swap_ = std::make_pair(type, index);
    }

    /**
     * @return true if and only if this separation has a special swap
     */

    bool has_special_swap() const
    {
      return special_swap_.first != 0;
    }

    /**
     * @return true if and only if this separation hsa a special row swap
     */

    bool has_special_row_swap() const
    {
      return special_swap_.first == 'r';
    }

    /**
     * @return The row / column index of a special swap
     */

    size_t get_special_swap_index() const
    {
      return special_swap_.second;
    }

    /**
     * Creates both components as matroids and representation matrices.
     *
     * @param matroid The original matroid
     * @param matrix Representation matrix of the given original matroid
     * @param upper_left_matroid New upper left component matroid
     * @param upper_left_matrix New upper left component representation matrix
     * @param lower_right_matroid New lower right component matroid
     * @param lower_right_matrix New lower right component representation matrix
     */

    template <typename MatroidType, typename MatrixType>
    void create_components(const MatroidType& matroid, const MatrixType& matrix, integer_matroid& upper_left_matroid,
        integer_matrix& upper_left_matrix, integer_matroid& lower_right_matroid, integer_matrix& lower_right_matrix) const
    {
      assert (is_valid());
      assert (matroid.size1() == matrix.size1());
      assert (matroid.size2() == matrix.size2());

      /// 1-separation
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
        /// 2-separation with lower left rank 1

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
        /// 2-separation with upper right rank 1

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
        /// 3-separation

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

} /* namespace tu */
