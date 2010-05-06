/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef MATRIX_TRANSPOSED_HPP_
#define MATRIX_TRANSPOSED_HPP_

#include "../config.h"
#include <boost/numeric/ublas/matrix_expression.hpp>

namespace tu {

  namespace detail {
    template <typename OrientationCategory>
    struct transpose_orientation
    {
      typedef boost::numeric::ublas::unknown_orientation_tag orientation_category;
    };

    template <>
    struct transpose_orientation <boost::numeric::ublas::row_major_tag>
    {
      typedef boost::numeric::ublas::column_major_tag orientation_category;
    };

    template <>
    struct transpose_orientation <boost::numeric::ublas::column_major_tag>
    {
      typedef boost::numeric::ublas::row_major_tag orientation_category;
    };
  }

  template <typename M>
  class matrix_transposed: public boost::numeric::ublas::matrix_expression <matrix_transposed <M> >
  {
  public:
    typedef matrix_transposed <M> self_type;
    typedef M matrix_type;
    typedef typename M::size_type size_type;
    typedef typename M::difference_type difference_type;
    typedef typename M::value_type value_type;
    typedef typename M::const_reference const_reference;
    typedef typename boost::mpl::if_ <boost::is_const <M>, typename M::const_reference, typename M::reference>::type reference;
    typedef typename boost::mpl::if_ <boost::is_const <M>, typename M::const_closure_type, typename M::closure_type>::type matrix_closure_type;
    typedef const self_type const_closure_type;
    typedef self_type closure_type;
    typedef typename detail::transpose_orientation <typename M::orientation_category>::orientation_category orientation_category;
    typedef typename M::storage_category storage_category;

  private:
    matrix_type& _data;

  public:

    matrix_transposed (M& matrix) :
      _data(matrix)
    {

    }

    // Accessors
    inline size_type size1 () const
    {
      return _data.size2();
    }

    inline size_type size2 () const
    {
      return _data.size1();
    }

    inline M& data ()
    {
      return _data;
    }

    // Element access
    inline const_reference operator () (size_type i, size_type j) const
    {
      return _data(j, i);
    }

    inline reference operator () (size_type i, size_type j)
    {
      return _data(j, i);
    }

    typedef boost::numeric::ublas::detail::indexed_iterator1 <self_type, typename matrix_type::iterator2::iterator_category> iterator1;
    typedef boost::numeric::ublas::detail::indexed_iterator2 <self_type, typename matrix_type::iterator1::iterator_category> iterator2;
    typedef boost::numeric::ublas::detail::indexed_const_iterator1 <self_type, typename matrix_type::const_iterator2::iterator_category>
        const_iterator1;
    typedef boost::numeric::ublas::detail::indexed_const_iterator2 <self_type, typename matrix_type::const_iterator1::iterator_category>
        const_iterator2;
  };

  template <typename MatrixType>
  inline matrix_transposed <MatrixType> view_matrix_transposed (MatrixType& matrix)
  {
    return matrix_transposed <MatrixType> (matrix);
  }

  template <typename MatrixType>
  inline void matrix_set_value (matrix_transposed <MatrixType>& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {
    matrix(row, column) = value;
  }

  template <typename MatrixType>
  inline void matrix_set_value (matrix_transposed <const MatrixType>& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {

  }

  template <typename MatrixType>
  inline void matrix_permute1 (matrix_transposed <MatrixType>& matrix, size_t index1, size_t index2)
  {
    matrix_permute2(matrix.data(), index1, index2);
  }

  template <typename MatrixType>
  inline void matrix_permute2 (matrix_transposed <MatrixType>& matrix, size_t index1, size_t index2)
  {
    matrix_permute1(matrix.data(), index1, index2);
  }

  template <typename MatrixType>
  inline void matrix_binary_pivot (matrix_transposed <MatrixType>& matrix, size_t i, size_t j)
  {
    matrix_binary_pivot(matrix.data(), j, i);
  }

}

#endif /* MATRIX_TRANSPOSED_HPP_ */
