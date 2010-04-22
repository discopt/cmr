/*
 * matrix_modified.hpp
 *
 *  Created on: Dec 20, 2009
 *      Author: xammy
 */

#ifndef MATRIX_MODIFIED_HPP_
#define MATRIX_MODIFIED_HPP_

#include "../config.h"

#include <boost/numeric/ublas/matrix_expression.hpp>

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

  matrix_modified (const MatrixType& matrix, modifier_type& modifier) :
    _data (matrix), _modifier (modifier)
  {

  }

  // Accessors
  inline size_type size1 () const
  {
    return _data.size1 ();
  }

  inline size_type size2 () const
  {
    return _data.size2 ();
  }

  inline MatrixType& data ()
  {
    return _data;
  }

  inline modifier_type& modifier ()
  {
    return _modifier;
  }

  // Element access
  inline value_type operator () (size_type i, size_type j) const
  {
    return _modifier (i, j, _data (i, j));
  }

};

#endif /* MATRIX_MODIFIED_HPP_ */
