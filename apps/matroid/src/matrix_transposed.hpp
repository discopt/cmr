/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef MATRIX_TRANSPOSED_HPP_
#define MATRIX_TRANSPOSED_HPP_

#include <boost/numeric/ublas/matrix_expression.hpp>

namespace unimod
{

  template <typename MatrixType>
  void matrix_binary_pivot(MatrixType& matrix, size_t i, size_t j);

  namespace detail
  {

    /// Helper struct to manage orientation tags

    template <typename OrientationCategory>
    struct transpose_orientation
    {
      typedef boost::numeric::ublas::unknown_orientation_tag orientation_category;
    };

    /// Helper struct to manage orientation tags

    template <>
    struct transpose_orientation <boost::numeric::ublas::row_major_tag>
    {
      typedef boost::numeric::ublas::column_major_tag orientation_category;
    };

    /// Helper struct to manage orientation tags

    template <>
    struct transpose_orientation <boost::numeric::ublas::column_major_tag>
    {
      typedef boost::numeric::ublas::row_major_tag orientation_category;
    };
  }

  /**
   * A matrix proxy which exchanges the meaning of rows and columns.
   */

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

    /**
     * Constructs the matrix proxy.
     *
     * @param matrix The original matrix
     */

    matrix_transposed(matrix_type& matrix) :
      _data(matrix)
    {

    }

    /**
     * @return Height of the matrix
     */

    inline size_type size1() const
    {
      return _data.size2();
    }

    /**
     * @return Width of the matrix
     */

    inline size_type size2() const
    {
      return _data.size1();
    }

    /**
     * @return Reference to the original matrix
     */

    inline matrix_type& data()
    {
      return _data;
    }

    /**
     * Read-only access operator
     *
     * @param i Row index
     * @param j Column index
     * @return original(column, row)
     */

    inline const_reference operator ()(size_type i, size_type j) const
    {
      return _data(j, i);
    }

    /**
     * Access operator
     *
     * @param i Row index
     * @param j Column index
     * @return original(column, row)
     */

    inline reference operator ()(size_type i, size_type j)
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

  /**
   * Creates a transpose proxy of a given matrix.
   *
   * @param matrix The original matrix
   * @return A matrix_transposed proxy
   */

  template <typename MatrixType>
  inline matrix_transposed <MatrixType> make_transposed_matrix(MatrixType& matrix)
  {
    return matrix_transposed <MatrixType> (matrix);
  }

  /**
   * Free function to set a matrix value of a permuted matrix.
   *
   * @param matrix The transposed matrix
   * @param row Row index
   * @param column Column index
   * @param value New value
   */

  template <typename MatrixType>
  inline void matrix_set_value(matrix_transposed <MatrixType>& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {
    matrix(row, column) = value;
  }

  /**
   * Dummy function for the version with writable orignal matrix.
   */

  template <typename MatrixType>
  inline void matrix_set_value(matrix_transposed <const MatrixType>& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {

  }

  /**
   * Free function to permute two rows of a transposed matrix.
   *
   * @param matrix The transposed matrix
   * @param index1 First index
   * @param index2 Second index
   */

  template <typename MatrixType>
  inline void matrix_permute1(matrix_transposed <MatrixType>& matrix, size_t index1, size_t index2)
  {
    matrix_permute2(matrix.data(), index1, index2);
  }

  /**
   * Free function to permute two columns of a transposed matrix.
   *
   * @param matrix The transposed matrix
   * @param index1 First index
   * @param index2 Second index
   */

  template <typename MatrixType>
  inline void matrix_permute2(matrix_transposed <MatrixType>& matrix, size_t index1, size_t index2)
  {
    matrix_permute1(matrix.data(), index1, index2);
  }

  /**
   * Free function to perform a pivot on a transposed matrix.
   *
   * @param matrix The transposed matrix
   * @param i Row index
   * @param j Column index
   */

  template <typename MatrixType>
  inline void matrix_binary_pivot(matrix_transposed <MatrixType>& matrix, size_t i, size_t j)
  {
    matrix_binary_pivot(matrix.data(), j, i);
  }

}

#endif /* MATRIX_TRANSPOSED_HPP_ */
