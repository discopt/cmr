/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef MATROID_TRANSPOSED_HPP_
#define MATROID_TRANSPOSED_HPP_

#include "matroid.hpp"

namespace tu
{

  /**
   * A matroid proxy which exchanges the meaning of bases and cobases.
   */

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

    /**
     * Constructs the matroid proxy.
     *
     * @param matroid The original matroid
     */

    matroid_transposed(MatroidType& matroid) :
      _matroid(matroid)
    {

    }

    /**
     * @return Height, i.e. size of each base
     */

    inline size_type size1() const
    {
      return _matroid.size2();
    }

    /**
     * @return Width, i.e. size of each cobase
     */

    inline size_type size2() const
    {
      return _matroid.size1();
    }

    /**
     * @param index A row index
     * @return The corresponding matroid element
     */

    inline reference_type name1(size_type index)
    {
      return _matroid.name2(index);
    }

    /**
     * @param index A row index
     * @return The corresponding matroid element
     */

    inline const_reference_type name1(size_type index) const
    {
      return _matroid.name2(index);
    }

    /**
     * @param index A column  index
     * @return The corresponding matroid element
     */

    inline reference_type name2(size_type index)
    {
      return _matroid.name1(index);
    }

    /**
     * @param index A column  index
     * @return The corresponding matroid element
     */

    inline const_reference_type name2(size_type index) const
    {
      return _matroid.name1(index);
    }

    /**
     * @return Reference to the orginal matroid
     */

    inline matroid_type& data()
    {
      return _matroid;
    }

  private:
    MatroidType& _matroid;

  };

  /**
   * Creates a transpose proxy of a given matroid, i.e. the dual matroid.
   *
   * @param matroid The original matroid
   * @return A matroid_transposed proxy
   */

  template <typename MatroidType>
  inline matroid_transposed <MatroidType> make_transposed_matroid(MatroidType& matroid)
  {
    return matroid_transposed <MatroidType> (matroid);
  }

  /**
   * Free function to permute two rows.
   *
   * @param matroid The given matroid
   * @param index1 First row index
   * @param index2 Second row index
   */

  template <typename MatroidType>
  inline void matroid_permute1(matroid_transposed <MatroidType>& matroid, size_t index1, size_t index2)
  {
    matroid_permute2(matroid.data(), index1, index2);
  }

  /**
   * Free function to permute two columns.
   *
   * @param matroid The given matroid
   * @param index1 First column index
   * @param index2 Second column index
   */

  template <typename MatroidType>
  inline void matroid_permute2(matroid_transposed <MatroidType>& matroid, size_t index1, size_t index2)
  {
    matroid_permute1(matroid.data(), index1, index2);
  }

  /**
   * Free function to perform a pivot on a matroid.
   *
   * @param matroid The given matroid
   * @param i Row index
   * @param j Column index
   */

  template <typename MatroidType>
  void matroid_binary_pivot(matroid_transposed <MatroidType>& matroid, size_t i, size_t j)
  {
    matroid_binary_pivot(matroid.data(), j, i);
  }

}

#endif /* MATROID_TRANSPOSED_HPP_ */
