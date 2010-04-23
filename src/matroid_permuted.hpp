/*
 * matroid_permuted.hpp
 *
 *  Created on: Dec 20, 2009
 *      Author: xammy
 */

#ifndef MATROID_PERMUTED_HPP_
#define MATROID_PERMUTED_HPP_

#include "../config.h"
#include "matroid.hpp"

namespace tu {

  template <typename MatroidType>
  class matroid_permuted
  {
  public:
    typedef MatroidType matroid_type;
    typedef typename MatroidType::size_type size_type;
    typedef typename MatroidType::name_type name_type;
    typedef typename boost::mpl::if_ <boost::is_const <MatroidType>, typename MatroidType::const_reference_type, typename MatroidType::reference_type>::type
        reference_type;
    typedef const name_type& const_reference_type;
    typedef matroid_permuted <name_type> self_type;
    typedef permutation permutation_type;

    matroid_permuted (MatroidType& data) :
      _data (data), _perm1 (data.size1 ()), _perm2 (data.size2 ())
    {

    }

    // Accessors
    inline size_type size1 () const
    {
      return _perm1.size ();
    }

    inline size_type size2 () const
    {
      return _perm2.size ();
    }

    inline const permutation_type& perm1 () const
    {
      return _perm1;
    }

    inline permutation_type& perm1 ()
    {
      return _perm1;
    }

    inline const permutation_type& perm2 () const
    {
      return _perm2;
    }

    inline permutation_type& perm2 ()
    {
      return _perm2;
    }

    // Element access
    inline reference_type name1 (size_type index)
    {
      return _data.name1 (_perm1 (index));
    }

    inline const_reference_type name1 (size_type index) const
    {
      return _data.name1 (_perm1 (index));
    }

    inline reference_type name2 (size_type index)
    {
      return _data.name2 (_perm2 (index));
    }

    inline const_reference_type name2 (size_type index) const
    {
      return _data.name2 (_perm2 (index));
    }

    inline matroid_type& data ()
    {
      return _data;
    }

  private:
    MatroidType& _data;
    permutation_type _perm1;
    permutation_type _perm2;

  };

  template <typename MatroidType>
  inline void matroid_permute1 (matroid_permuted <MatroidType>& matroid, size_t index1, size_t index2)
  {
    matroid.perm1 ().swap (index1, index2);
  }

  template <typename MatroidType>
  inline void matroid_permute2 (matroid_permuted <MatroidType>& matroid, size_t index1, size_t index2)
  {
    matroid.perm2 ().swap (index1, index2);
  }

  template <typename MatroidType>
  void matroid_binary_pivot (matroid_permuted <MatroidType>& matroid, size_t i, size_t j)
  {
    matroid_binary_pivot (matroid.data (), matroid.perm1 () (i), matroid.perm2 () (j));
  }

//template<typename MatroidType, typename MatrixType, typename ElementLess>
//inline void matroid_reorder_rows (matroid_permuted<MatroidType>& permuted_matroid,
//        matrix_permuted<MatrixType>& permuted_matrix, size_t row_first, size_t row_beyond, size_t column_first,
//        size_t column_beyond, ElementLess element_less)
//{
//    matrix_reorder_rows (permuted_matrix, row_first, row_beyond, column_first, column_beyond, element_less);
//    permuted_matroid.perm1 () = permuted_matrix.perm1 ();
//}
//
//template<typename MatroidType, typename MatrixType, typename ElementLess>
//inline void matroid_reorder_columns (matroid_permuted<MatroidType>& permuted_matroid,
//        matrix_permuted<MatrixType>& permuted_matrix, size_t row_first, size_t row_beyond, size_t column_first,
//        size_t column_beyond, ElementLess element_less)
//{
//    matrix_reorder_columns (permuted_matrix, row_first, row_beyond, column_first, column_beyond, element_less);
//    permuted_matroid.perm2 () = permuted_matrix.perm2 ();
//}

}

#endif /* MATROID_PERMUTED_HPP_ */
