#pragma once

#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/detail/iterator.hpp>

#include "permutations.hpp"

namespace tu
{

  /**
   * A matrix proxy with permuted rows and columns.
   */

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

    /**
     * Constructs the matrix proxy.
     *
     * @param matrix The original matrix
     */

    matrix_permuted(matrix_type& matrix) :
      _data(matrix), _perm1(matrix.size1()), _perm2(matrix.size2())
    {

    }

    /**
     * @return Height of the matrix
     */

    inline size_type size1() const
    {
      return _perm1.size();
    }

    /**
     * @return Width of the matrix
     */

    inline size_type size2() const
    {
      return _perm2.size();
    }

    /**
     * @return Reference to the row permutation
     */

    inline const permutation_type& perm1() const
    {
      return _perm1;
    }

    /**
     * @return Reference to the row permutation
     */

    inline permutation_type& perm1()
    {
      return _perm1;
    }

    /**
     * @return Reference to the column permutation
     */

    inline const permutation_type& perm2() const
    {
      return _perm2;
    }

    /**
     * @return Reference to the column permutation
     */

    inline permutation_type& perm2()
    {
      return _perm2;
    }

    /**
     * Read-only access operator
     *
     * @param i Row index
     * @param j Column index
     * @return original(row-permutation(row), column-permutation(column))
     */

    inline const_reference operator ()(size_type i, size_type j) const
    {
      return _data(_perm1(i), _perm2(j));
    }

    /**
     * Access operator
     *
     * @param i Row index
     * @param j Column index
     * @return original(row-permutation(row), column-permutation(column))
     */

    inline reference operator ()(size_type i, size_type j)
    {
      return _data(_perm1(i), _perm2(j));
    }

    /**
     * @return Reference to original matrix
     */

    inline matrix_type& data()
    {
      return _data;
    }

    /**
     * @return Read-only reference to original matrix
     */

    inline const matrix_type& data() const
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

  /**
   * Free function to set a matrix value of a permuted matrix.
   *
   * @param matrix The permuted matrix
   * @param row Row index
   * @param column Column index
   * @param value New value
   */

  template <typename MatrixType>
  inline void matrix_set_value(matrix_permuted <MatrixType>& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {
    matrix(row, column) = value;
  }

  /**
   * Dummy function for the version with writable orignal matrix.
   */

  template <typename MatrixType>
  inline void matrix_set_value(matrix_permuted <const MatrixType>& matrix, size_t row, size_t column, typename MatrixType::value_type value)
  {
    CMR_UNUSED(matrix);
    CMR_UNUSED(row);
    CMR_UNUSED(column);
    CMR_UNUSED(value);
    assert (false);
  }

  /**
   * Free function to permute two rows of a permuted matrix.
   *
   * @param matrix The permuted matrix
   * @param index1 First index
   * @param index2 Second index
   */

  template <typename MatrixType>
  inline void matrix_permute1(matrix_permuted <MatrixType>& matrix, size_t index1, size_t index2)
  {
    matrix.perm1().swap(index1, index2);
  }

  /**
   * Free function to permute two columns of a permuted matrix.
   *
   * @param matrix The permuted matrix
   * @param index1 First index
   * @param index2 Second index
   */

  template <typename MatrixType>
  inline void matrix_permute2(matrix_permuted <MatrixType>& matrix, size_t index1, size_t index2)
  {
    matrix.perm2().swap(index1, index2);
  }

  /**
   * Free function to perform a pivot on a permuted matrix.
   *
   * @param matrix The permuted matrix
   * @param i Row index
   * @param j Column index
   */

  template <typename MatrixType>
  inline void matrix_binary_pivot(matrix_permuted <MatrixType>& matrix, size_t i, size_t j)
  {
    matrix_binary_pivot(matrix.data(), matrix.perm1()(i), matrix.perm2()(j));
  }

  template <typename MatrixType>
  inline permutation matrix_get_perm1(matrix_permuted <MatrixType>& matrix)
  {
    return matrix.perm1();
  }

  template <typename MatrixType>
  inline permutation matrix_get_perm2(matrix_permuted <MatrixType>& matrix)
  {
    return matrix.perm2();
  }

  template <typename MatrixType>
  inline permutation matrix_get_perm1(matrix_transposed <MatrixType>& matrix)
  {
    return matrix_get_perm2(matrix.data());
  }

  template <typename MatrixType>
  inline permutation matrix_get_perm2(matrix_transposed <MatrixType>& matrix)
  {
    return matrix_get_perm1(matrix.data());
  }

  template <typename MatrixType>
  inline void matrix_set_perm1(matrix_permuted <MatrixType>& matrix, const permutation& permutation)
  {
    matrix.perm1() = permutation;
  }

  template <typename MatrixType>
  inline void matrix_set_perm2(matrix_permuted <MatrixType>& matrix, const permutation& permutation)
  {
    matrix.perm2() = permutation;
  }

  template <typename MatrixType>
  inline void matrix_set_perm1(matrix_transposed <MatrixType>& matrix, const permutation& permutation)
  {
    matrix_set_perm2(matrix.data(), permutation);
  }

  template <typename MatrixType>
  inline void matrix_set_perm2(matrix_transposed <MatrixType>& matrix, const permutation& permutation)
  {
    matrix_set_perm1(matrix.data(), permutation);
  }
} /* namespace tu */
