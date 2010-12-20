/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef COMPARATORS_HPP_
#define COMPARATORS_HPP_

namespace tu
{

  /**
   * This functor can be used to test a vector for being non-zero.
   * Repeated calls with parameter determine the result which
   * can be obtained by a call without a parameter.
   */

  struct is_non_zero
  {
    /**
     * Constructs the functor.
     */

    is_non_zero() :
      _valid(false)
    {

    }

    /**
     * Each call tests one entry of the vector.
     *
     * @param value A value to be tested
     */

    void operator()(int value)
    {
      if (value != 0)
        _valid = true;
    }

    /**
     * Resets the result to an empty vector.
     *
     * @return The result of the tests from the last call to this.
     */

    bool operator()()
    {
      bool result = _valid;
      _valid = false;
      return result;
    }

  private:
    bool _valid;
  };

  /**
   * This functor can be used to test a vector to contain only non-zeros.
   * Repeated calls with parameter determine the result which
   * can be obtained by a call without a parameter.
   */

  struct is_all_ones
  {
    /**
     * Constructs the functor.
     */

    is_all_ones() :
      _valid(true)
    {

    }

    /**
     * Each call tests one entry of the vector.
     *
     * @param value A value to be tested
     */

    void operator()(int value)
    {
      if (value == 0)
        _valid = false;
    }

    /**
     * Resets the result to an empty vector.
     *
     * @return The result of the tests from the last call to this.
     */

    bool operator()()
    {
      bool result = _valid;
      _valid = true;
      return result;
    }

  private:
    bool _valid;
  };

  /**
   * A wrapper for an element-comparator which can be used
   * to sort an indirection vector according to the data
   * it points to.
   */

  template <typename T, typename Less = std::less <T> >
  struct vector_less
  {
    /**
     * Constructs the comparator.
     *
     * @param data The vector
     * @param less An optional element-comparator
     */

    vector_less(const std::vector <T>& data, Less less = Less()) :
      _less(less), data_(data)
    {

    }

    /**
     * Compares two indices by comparing their elements.
     *
     * @param i First index
     * @param j Second index
     * @return true if and only if the element at the first index is smaller than the element at the second index
     */

    bool operator ()(size_t i, size_t j)
    {
      return _less(data_[i], data_[j]);
    }

  private:
    Less _less;
    const std::vector <T>& data_;
  };

}

#endif /* COMPARATORS_HPP_ */
