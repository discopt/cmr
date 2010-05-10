/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef MATROID_HPP_
#define MATROID_HPP_

#include <iomanip>
#include <iostream>
#include <set>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "permutations.hpp"

namespace tu {

  /**
   * A matroid class which maintains the index-to-element relation.
   * It should mostly be associated with a representation matrix.
   */

  template <typename NameType>
  class matroid
  {
  public:
    typedef NameType name_type;
    typedef size_t size_type;
    typedef std::vector <name_type> name_vector_type;
    typedef name_type& reference_type;
    typedef const name_type& const_reference_type;
    typedef matroid <name_type> self_type;

    /**
     * Constructs a matroid by copying the matroid elements from the row/column vectors.
     *
     * @param names1 Row names
     * @param names2 Column names
     */

    matroid (const name_vector_type& names1, const name_vector_type& names2) :
      _names1(names1), _names2(names2)
    {

    }

    /**
     * Constructs a matroid with given height and width. The names are set in default manner.
     *
     * @param size1 Size of a base
     * @param size2 Size of a cobase
     */

    matroid (size_t size1 = 0, size_t size2 = 0)
    {
      resize(size1, size2);
    }

    /**
     * Resizes a matroid.
     *
     * @param size1 New height, i.e. size of each base
     * @param size2 New width, i.e. size of each cobase
     */

    void resize (size_t size1, size_t size2)
    {
      _names1.resize(size1);
      for (size_t i = 0; i < size1; ++i)
        _names1[i] = -(i + 1);
      _names2.resize(size2);
      for (size_t i = 0; i < size2; ++i)
        _names2[i] = (i + 1);
    }

    /**
     * @return Height, i.e. size of each base
     */

    inline size_t size1 () const
    {
      return _names1.size();
    }

    /**
     * @return Width, i.e. size of each cobase
     */

    inline size_t size2 () const
    {
      return _names2.size();
    }

    /**
     * @param index A row index
     * @return The corresponding matroid element
     */

    inline name_type& name1 (size_t index)
    {
      return _names1[index];
    }

    /**
     * @param index A row index
     * @return The corresponding matroid element
     */

    inline const name_type& name1 (size_t index) const
    {
      return _names1[index];
    }

    /**
     * @param index A column  index
     * @return The corresponding matroid element
     */

    inline name_type& name2 (size_t index)
    {
      return _names2[index];
    }

    /**
     * @param index A column index
     * @return The corresponding matroid element
     */

    inline const name_type& name2 (size_t index) const
    {
      return _names2[index];
    }

    /**
     * @return A set of all matroid elements
     */

    inline std::set <NameType> get_elements () const
    {
      std::set <NameType> result;

      std::copy(_names1.begin(), _names1.end(), std::inserter(result, result.end()));
      std::copy(_names2.begin(), _names2.end(), std::inserter(result, result.end()));

      return result;
    }

  private:
    name_vector_type _names1;
    name_vector_type _names2;
  };

  /// A matroid with integer names
  typedef matroid <int> integer_matroid;

  /**
   * Free function to permute two rows.
   *
   * @param matroid The given matroid
   * @param index1 First row index
   * @param index2 Second row index
   */

  template <typename NameType>
  inline void matroid_permute1 (matroid <NameType>& matroid, size_t index1, size_t index2)
  {
    std::swap(matroid.name1(index1), matroid.name1(index2));
  }

  /**
   * Free function to permute two rows of a matroid and its representation matrix.
   *
   * @param matroid The given matroid
   * @param matrix Representation matrix of the given matroid
   * @param index1 First row index
   * @param index2 Second row index
   */

  template <typename MatroidType, typename MatrixType>
  inline void matroid_permute1 (MatroidType& matroid, MatrixType& matrix, size_t index1, size_t index2)
  {
    matroid_permute1(matroid, index1, index2);
    matrix_permute1(matrix, index1, index2);
  }

  /**
   * Free function to permute two columns.
   *
   * @param matroid The given matroid
   * @param index1 First column index
   * @param index2 Second column index
   */

  template <typename NameType>
  inline void matroid_permute2 (matroid <NameType>& matroid, size_t index1, size_t index2)
  {
    std::swap(matroid.name2(index1), matroid.name2(index2));
  }

  /**
   * Free function to permute two columns of a matroid and its representation matrix.
   *
   * @param matroid The given matroid
   * @param matrix Representation matrix of the given matroid
   * @param index1 First column index
   * @param index2 Second column index
   */

  template <typename MatroidType, typename MatrixType>
  inline void matroid_permute2 (MatroidType& matroid, MatrixType& matrix, size_t index1, size_t index2)
  {
    matroid_permute2(matroid, index1, index2);
    matrix_permute2(matrix, index1, index2);
  }

  /**
   * Free function to perform a pivot on a matroid.
   *
   * @param matroid The given matroid
   * @param i Row index
   * @param j Column index
   */

  template <typename NameType>
  void matroid_binary_pivot (matroid <NameType>& matroid, size_t i, size_t j)
  {
    std::swap(matroid.name1(i), matroid.name2(j));
  }

  /**
   * Free function to perform a pivot on a matroid and its representation matrix.
   *
   * @param matroid The given matroid
   * @param matrix Representation matrix of the given matroid
   * @param i Row index
   * @param j Column index
   */

  template <typename MatroidType, typename MatrixType>
  void matroid_binary_pivot (MatroidType& matroid, MatrixType& matrix, size_t i, size_t j)
  {
    matroid_binary_pivot(matroid, i, j);
    matrix_binary_pivot(matrix, i, j);
  }

  /**
   * Prints a matroid with its representation matrix.
   *
   * @param matroid The given matrix
   * @param matrix Representation matrix of the given matroid
   */

  template <typename MatroidType, typename MatrixType>
  void matroid_print (const MatroidType& matroid, const MatrixType& matrix)
  {
    assert (matroid.size1() == matrix.size1());
    assert (matroid.size2() == matrix.size2());

    std::cout << "     ";
    for (size_t column = 0; column < matroid.size2(); ++column)
    {
      std::cout << std::setw(3) << matroid.name2(column);
    }
    std::cout << "\n    ";
    for (size_t column = 0; column < matroid.size2(); ++column)
    {
      std::cout << "--";
    }
    for (size_t row = 0; row < matroid.size1(); ++row)
    {
      std::cout << "\n" << std::setw(3) << matroid.name1(row) << " |";
      for (size_t column = 0; column < matroid.size2(); ++column)
      {
        std::cout << std::setw(3) << matrix(row, column);
      }
    }
    std::cout << std::endl;
  }

  /**
   * Prints a matroid minor with its representation matrix.
   *
   * @param matroid The given matrix
   * @param matrix Representation matrix of the given matroid
   * @param height Number of rows of the minor
   * @param width Number of columns of the minor
   */

  template <typename MatroidType, typename MatrixType>
  void matroid_print_minor (const MatroidType& matroid, const MatrixType& matrix, size_t height, size_t width)
  {
    assert (matroid.size1() == matrix.size1());
    assert (matroid.size2() == matrix.size2());

    std::cout << "     ";
    for (size_t column = 0; column < width; ++column)
    {
      std::cout << std::setw(3) << matroid.name2(column);
    }
    std::cout << "\n    ";
    for (size_t column = 0; column < width; ++column)
    {
      std::cout << "--";
    }
    for (size_t row = 0; row < height; ++row)
    {
      std::cout << "\n" << std::setw(3) << matroid.name1(row) << " |";
      for (size_t column = 0; column < width; ++column)
      {
        std::cout << std::setw(3) << matrix(row, column);
      }
    }
    std::cout << std::endl;
  }

  /**
   * Returns all matroid elements.
   *
   * @param matroid The given matroid
   * @return A set of matroid elements
   */

  template <typename MatroidType>
  std::set <int> matroid_elements (const MatroidType& matroid)
  {
    std::set <int> elements;
    for (size_t row = 0; row < matroid.size1(); ++row)
      elements.insert(matroid.name1(row));
    for (size_t column = 0; column < matroid.size2(); ++column)
      elements.insert(matroid.name2(column));

    return elements;
  }

}

#endif /* MATROID_HPP_ */
