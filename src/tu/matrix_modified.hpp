#pragma once

#include <boost/numeric/ublas/matrix_expression.hpp>

namespace tu
{

  /**
   * A matrix proxy which can modify arbitrary entries by attaching a separator functor to it.
   */

  template <typename MatrixType, typename ModifierType>
  class matrix_modified: public boost::numeric::ublas::matrix_expression <matrix_modified <MatrixType, ModifierType> >
  {
  public:
    typedef matrix_modified <MatrixType, ModifierType> self_type;
    typedef MatrixType matrix_type;
    typedef typename MatrixType::size_type size_type;
    typedef typename MatrixType::difference_type difference_type;
    typedef typename ModifierType::value_type value_type;
    typedef typename MatrixType::const_reference const_reference;
    typedef typename MatrixType::const_reference reference;
    typedef typename MatrixType::const_closure_type matrix_closure_type;
    typedef const self_type const_closure_type;
    typedef self_type closure_type;

    typedef ModifierType modifier_type;

  private:
    const matrix_type& _data;
    modifier_type& _modifier;

  public:

    /**
     * Constructs the modified matrix.
     *
     * @param matrix The orginal matrix
     * @param modifier Modifier to calculate each entry
     */

    matrix_modified(const MatrixType& matrix, modifier_type& modifier) :
      _data(matrix), _modifier(modifier)
    {

    }

    /**
     * @return Height of the matrix
     */

    inline size_type size1() const
    {
      return _data.size1();
    }

    /**
     * @return Width of the matrix
     */

    inline size_type size2() const
    {
      return _data.size2();
    }

    /**
     * @return Reference to the proxied matrix
     */

    inline MatrixType& data()
    {
      return _data;
    }

    /**
     * @return The modifier instance
     */

    inline modifier_type& modifier()
    {
      return _modifier;
    }

    /**
     * Access operator which calls the modifier to calculate the new value.
     *
     * @param i Row index
     * @param j Column index
     * @return Modified entry
     */

    inline value_type operator ()(size_type i, size_type j) const
    {
      return _modifier(i, j, _data(i, j));
    }

  };

} /* namespace tu */
