/*
 * matrix_permuted.hpp
 *
 *  Created on: Dec 20, 2009
 *      Author: xammy
 */

#ifndef MATRIX_PERMUTED_HPP_
#define MATRIX_PERMUTED_HPP_

#include "../config.h"
#include "permutations.hpp"

#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/detail/iterator.hpp>

namespace tu {

  template <typename M>
  class matrix_permuted: public boost::numeric::ublas::matrix_expression <matrix_permuted <M> >
  {
  public:
    typedef matrix_permuted <M> self_type;
    typedef M matrix_type;
    typedef typename M::size_type size_type;
    typedef typename M::difference_type difference_type;
    typedef typename M::value_type value_type;
    typedef typename M::const_reference const_reference;
    typedef typename boost::mpl::if_ <boost::is_const <M>, typename M::const_reference, typename M::reference>::type reference;
    typedef typename boost::mpl::if_ <boost::is_const <M>, typename M::const_closure_type, typename M::closure_type>::type matrix_closure_type;
    typedef const self_type const_closure_type;
    typedef self_type closure_type;
    typedef permutation permutation_type;
    typedef typename M::orientation_category orientation_category;
    typedef typename M::storage_category storage_category;

  private:
    matrix_type& _data;
    permutation_type _perm1;
    permutation_type _perm2;

  public:

    matrix_permuted (M& matrix) :
      _data (matrix), _perm1 (matrix.size1 ()), _perm2 (matrix.size2 ())
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
    inline const_reference operator () (size_type i, size_type j) const
    {
      return _data (_perm1 (i), _perm2 (j));
    }

    inline reference operator () (size_type i, size_type j)
    {
      return _data (_perm1 (i), _perm2 (j));
    }

    inline matrix_type& data ()
    {
      return _data;
    }

    inline const matrix_type& data () const
    {
      return _data;
    }

    typedef boost::numeric::ublas::detail::indexed_iterator1 <self_type, typename matrix_type::iterator1::iterator_category> iterator1;
    typedef boost::numeric::ublas::detail::indexed_iterator2 <self_type, typename matrix_type::iterator2::iterator_category> iterator2;
    typedef boost::numeric::ublas::detail::indexed_const_iterator1 <self_type, typename matrix_type::const_iterator1::iterator_category>
        const_iterator1;
    typedef boost::numeric::ublas::detail::indexed_const_iterator2 <self_type, typename matrix_type::const_iterator2::iterator_category>
        const_iterator2;
  };

  template <typename MatrixType>
  inline void matrix_set_value (matrix_permuted <MatrixType>& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {
    matrix (row, column) = value;
  }

  template <typename MatrixType>
  inline void matrix_set_value (matrix_permuted <const MatrixType>& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {

  }

  template <typename MatrixType>
  inline void matrix_permute1 (matrix_permuted <MatrixType>& matrix, size_t index1, size_t index2)
  {
    matrix.perm1 ().swap (index1, index2);
  }

  template <typename MatrixType>
  inline void matrix_permute2 (matrix_permuted <MatrixType>& matrix, size_t index1, size_t index2)
  {
    matrix.perm2 ().swap (index1, index2);
  }

  template <typename MatrixType>
  inline void matrix_binary_pivot (matrix_permuted <MatrixType>& matrix, size_t i, size_t j)
  {
    matrix_binary_pivot (matrix.data (), matrix.perm1 () (i), matrix.perm2 () (j));
  }

}

#endif /* MATRIX_PERMUTED_HPP_ */
