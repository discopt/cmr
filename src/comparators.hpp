/*
 * comparators.hpp
 *
 *  Created on: Jan 24, 2010
 *      Author: xammy
 */

#ifndef COMPARATORS_HPP_
#define COMPARATORS_HPP_

#include "../config.h"

namespace tu {

  struct is_non_zero
  {
    is_non_zero () :
      _valid (false)
    {

    }

    void operator() (int value)
    {
      if (value != 0)
        _valid = true;
    }

    bool operator() ()
    {
      bool result = _valid;
      _valid = false;
      return result;
    }

  private:
    bool _valid;
  };

  struct is_all_ones
  {
    is_all_ones () :
      _valid (true)
    {

    }

    void operator() (int value)
    {
      if (value == 0)
        _valid = false;
    }

    bool operator() ()
    {
      bool result = _valid;
      _valid = true;
      return result;
    }

  private:
    bool _valid;
  };

  template <typename T, typename Less = std::less <T> >
  struct vector_less
  {
    vector_less (const std::vector <T>& data, Less less = Less ()) :
      less_ (less), data_ (data)
    {

    }

    bool operator () (size_t i, size_t j)
    {
      return less_ (data_[i], data_[j]);
    }

  private:
    Less less_;
    const std::vector <T>& data_;
  };

}

#endif /* COMPARATORS_HPP_ */
