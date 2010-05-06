/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef MATROID_TRANSPOSED_HPP_
#define MATROID_TRANSPOSED_HPP_

#include "../config.h"
#include "matroid.hpp"

namespace tu {

  template <typename MatroidType>
  class matroid_transposed
  {
  public:
    typedef MatroidType matroid_type;
    typedef typename MatroidType::size_type size_type;
    typedef typename MatroidType::name_type name_type;
    typedef typename boost::mpl::if_ <boost::is_const <MatroidType>, typename MatroidType::const_reference_type, typename MatroidType::reference_type>::type
        reference_type;
    typedef const name_type& const_reference_type;
    typedef matroid_transposed <name_type> self_type;
    typedef permutation permutation_type;

    matroid_transposed (MatroidType& matroid) :
      matroid_(matroid)
    {

    }

    // Accessors
    inline size_type size1 () const
    {
      return matroid_.size2();
    }

    inline size_type size2 () const
    {
      return matroid_.size1();
    }

    // Element access
    inline reference_type name1 (size_type index)
    {
      return matroid_.name2(index);
    }

    inline const_reference_type name1 (size_type index) const
    {
      return matroid_.name2(index);
    }

    inline reference_type name2 (size_type index)
    {
      return matroid_.name1(index);
    }

    inline const_reference_type name2 (size_type index) const
    {
      return matroid_.name1(index);
    }

    inline matroid_type& data ()
    {
      return matroid_;
    }

  private:
    MatroidType& matroid_;

  };

  template <typename MatroidType>
  inline matroid_transposed <MatroidType> view_matroid_transposed (MatroidType& matroid)
  {
    return matroid_transposed <MatroidType> (matroid);
  }

  template <typename MatroidType>
  inline void matroid_permute1 (matroid_transposed <MatroidType>& matroid, size_t index1, size_t index2)
  {
    matroid_permute2(matroid.data(), index1, index2);
  }

  template <typename MatroidType>
  inline void matroid_permute2 (matroid_transposed <MatroidType>& matroid, size_t index1, size_t index2)
  {
    matroid_permute1(matroid.data(), index1, index2);
  }

  template <typename MatroidType>
  void matroid_binary_pivot (matroid_transposed <MatroidType>& matroid, size_t i, size_t j)
  {
    matroid_binary_pivot(matroid.data(), j, i);
  }

}

#endif /* MATROID_TRANSPOSED_HPP_ */
