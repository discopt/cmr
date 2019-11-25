/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#include <tu/total_unimodularity.hpp>

#include "algorithm.hpp"
#include "matroid.hpp"
#include "violator_search.hpp"
#include "signing.hpp"
#include "logger.hpp"

namespace tu
{

  /**
   * Returns a decomposition of a given binary matroid into a k-sum-decomposition (k=1,2,3)
   * in graphic, cographic, R10 and maybe irregular components.
   *
   * @param matrix Representation matrix of a binary matroid
   * @param level Log level
   * @return Root of decomposition tree
   */

  decomposed_matroid* decompose_binary_matroid(const integer_matrix& matrix, log_level level)
  {
    logger log(level);

    if (!is_zero_one_matrix(matrix))
      return NULL;

    integer_matrix worker_matrix(matrix);
    integer_matroid worker_matroid(worker_matrix.size1(), worker_matrix.size2());

    return decompose_binary_matroid(worker_matroid, worker_matrix, matroid_element_set(), true, log).second;
  }

  /**
   * Tests for total unimodularity without certificates.
   * A matrix is totally unimodular if and only if every square submatrix
   * has a determinant of -1, 0 or +1.
   *
   * @param matrix The matrix to be tested
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  bool is_totally_unimodular(const integer_matrix& matrix, log_level level)
  {
    logger log(level);

    if (log.is_progressive())
    {
      log.line() << "(" << matrix.size1() << " x " << matrix.size2() << ")";
      std::cout << log;
    }

    /// Test for being a -1,0,+1 matrix
    if (!is_zero_plus_minus_one_matrix(matrix))
    {
      if (log.is_progressive())
      {
        log.line() << " NOT -1/0/+1";
        std::cout << log << std::endl;
      }
      else if (log.is_verbose())
      {
        std::cout << "Given " << matrix.size1() << " x " << matrix.size2() << " matrix does not contain only -1,0 and +1 entries." << std::endl;
      }

      return false;
    }

    if (log.is_progressive())
    {
      log.line() << " -1/0/+1 OK";
      std::cout << log;
    }
    else if (log.is_verbose())
    {
      std::cout << "Given " << matrix.size1() << " x " << matrix.size2() << " matrix contains only -1,0 and +1 entries." << std::endl;
    }

    /// Signing test
    if (!is_signed_matrix(matrix))
    {
      if (log.is_progressive())
      {
        log.line() << ", SIGNING FAILED\n";
        std::cout << log << std::endl;
      }
      else if (log.is_verbose())
      {
        std::cout << "The matrix is not its signed version." << std::endl;
      }

      return false;
    }

    if (log.is_progressive())
    {
      log.clear();
      std::cout << ", SIGNING OK" << std::endl;
    }
    else if (log.is_verbose())
    {
      std::cout << "The matrix is its signed version.\n" << std::endl;
    }

    /// Decomposition of matroid represented by support matrix
    integer_matrix worker_matrix;
    if (matrix.size1() < matrix.size2())
      worker_matrix = matrix;
    else
    {
      if (log.is_verbose())
      {
        std::cout << "Working on transposed matrix. Graphs and cographs are interchanged!\n" << std::endl;
      }
      else if (log.is_progressive())
      {
         log.line() << "(" << matrix.size1() << " x " << matrix.size2() << ") TRANSPOSING (graphs <-> cographs)" << std::endl;
         std::cout << log;
      }
      worker_matrix = make_transposed_matrix(matrix);
    }

    integer_matroid worker_matroid(worker_matrix.size1(), worker_matrix.size2());
    support_matrix(worker_matrix);

    std::pair <bool, decomposed_matroid*> result = decompose_binary_matroid(worker_matroid, worker_matrix, matroid_element_set(), false, log);

    assert (result.second == NULL);

    return result.first;
  }

  /**
   * Tests for total unimodularity with a positive certificate.
   * A matrix is totally unimodular if and only if every square submatrix
   * has a determinant of -1, 0 or +1.
   * If the matrix has this property, the routine returns a decomposition of
   * the underlying binary matroid into a k-sum-decomposition (k=1,2,3)
   * in graphic, cographic, R10 and maybe irregular components.
   *
   * @param matrix The matrix to be tested
   * @param decomposition Returns the root of the decomposition tree
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  bool is_totally_unimodular(const integer_matrix& matrix, decomposed_matroid*& decomposition, log_level level)
  {
    logger log(level);

    if (log.is_progressive())
    {
      log.line() << "(" << matrix.size1() << " x " << matrix.size2() << ")";
      std::cout << log;
    }

    /// Test for being a -1,0,+1 matrix
    if (!is_zero_plus_minus_one_matrix(matrix))
    {
      if (log.is_progressive())
      {
        log.line() << " NOT -1/0/+1";
        std::cout << log << std::endl;
      }
      else if (log.is_verbose())
      {
        std::cout << "Given " << matrix.size1() << " x " << matrix.size2() << " matrix does not contain only -1,0 and +1 entries." << std::endl;
      }

      decomposition = NULL;
      return false;
    }

    if (log.is_progressive())
    {
      log.line() << " -1/0/+1 OK";
      std::cout << log;
    }
    else if (log.is_verbose())
    {
      std::cout << "Given " << matrix.size1() << " x " << matrix.size2() << " matrix contains only -1,0 and +1 entries." << std::endl;
    }

    /// Signing test
    if (!is_signed_matrix(matrix))
    {
      if (log.is_progressive())
      {
        log.line() << ", SIGNING FAILED\n";
        std::cout << log << std::endl;
      }
      else if (log.is_verbose())
      {
        std::cout << "The matrix is not its signed version." << std::endl;
      }

      decomposition = NULL;
      return false;
    }

    if (log.is_progressive())
    {
      log.clear();
      std::cout << ", SIGNING OK" << std::endl;
    }
    else if (log.is_verbose())
    {
      std::cout << "The matrix is its signed version.\n" << std::endl;
    }

    /// Decomposition of matroid represented by support matrix
    integer_matrix worker_matrix;
    if (matrix.size1() < matrix.size2())
      worker_matrix = matrix;
    else
    {
      if (log.is_verbose())
      {
        std::cout << "Working on transposed matrix. Graphs and cographs are interchanged!\n" << std::endl;
      }
      worker_matrix = make_transposed_matrix(matrix);
    }

    integer_matroid worker_matroid(worker_matrix.size1(), worker_matrix.size2());
    support_matrix(worker_matrix);

    std::pair <bool, decomposed_matroid*> result = decompose_binary_matroid(worker_matroid, worker_matrix, matroid_element_set(), true, log);
    decomposition = result.second;

    return result.first;
  }

  /**
   * Tests for total unimodularity with certificates.
   * A matrix is totally unimodular if and only if every square submatrix
   * has a determinant of -1, 0 or +1.
   * If the matrix has this property, the routine returns a decomposition of
   * the underlying binary matroid into a k-sum-decomposition (k=1,2,3)
   * in graphic, cographic, R10 and maybe irregular components.
   * If the matrix does not have this property,the routine returns the indices
   * of a violating submatrix.
   *
   * @param matrix The matrix to be tested
   * @param decomposition Returns the root of the decomposition tree.
   * @param violator Returns violator indices
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  bool is_totally_unimodular(const integer_matrix& matrix, decomposed_matroid*& decomposition, submatrix_indices& violator, log_level level)
  {
    logger log(level);

    if (log.is_progressive())
    {
      log.line() << "(" << matrix.size1() << " x " << matrix.size2() << ")";
      std::cout << log;
    }

    /// Test for being a -1,0,+1 matrix
    std::pair <integer_matrix::size_type, integer_matrix::size_type> entry;
    if (!is_zero_plus_minus_one_matrix(matrix, entry))
    {
      if (log.is_progressive())
      {
        log.line() << " NOT -1/0/+1";
        std::cout << log << std::endl;
      }
      else if (log.is_verbose())
      {
        std::cout << "Given " << matrix.size1() << " x " << matrix.size2() << " matrix does not contain only -1,0 and +1 entries." << std::endl;
      }

      violator.rows = submatrix_indices::indirect_array_type(1);
      violator.rows[0] = entry.first;
      violator.columns = submatrix_indices::indirect_array_type(1);
      violator.columns[0] = entry.second;
      decomposition = NULL;
      return false;
    }

    if (log.is_progressive())
    {
      log.line() << " -1/0/+1 OK";
      std::cout << log;
    }
    else if (log.is_verbose())
    {
      std::cout << "Given " << matrix.size1() << " x " << matrix.size2() << " matrix contains only -1,0 and +1 entries." << std::endl;
    }

    /// Signing test
    if (!is_signed_matrix(matrix, violator))
    {
      if (log.is_progressive())
      {
        log.line() << ", SIGNING FAILED\n";
        std::cout << log << std::endl;
      }
      else if (log.is_verbose())
      {
        std::cout << "The matrix is not its signed version." << std::endl;
      }

      decomposition = NULL;
      assert (violator.rows.size() == violator.columns.size());
      return false;
    }

    if (log.is_progressive())
    {
      log.clear();
      std::cout << ", SIGNING OK" << std::endl;
    }
    else if (log.is_verbose())
    {
      std::cout << "The matrix is its signed version.\n" << std::endl;
    }

    /// Decomposition of matroid represented by support matrix
    integer_matrix worker_matrix;
    if (matrix.size1() < matrix.size2())
      worker_matrix = matrix;
    else
    {
      if (log.is_verbose())
      {
        std::cout << "Working on transposed matrix. Graphs and cographs are interchanged!\n" << std::endl;
      }
      worker_matrix = make_transposed_matrix(matrix);
    }

    integer_matroid worker_matroid(worker_matrix.size1(), worker_matrix.size2());
    support_matrix(worker_matrix);
    bool is_tu;
    boost::tie(is_tu, decomposition) = decompose_binary_matroid(worker_matroid, worker_matrix, matroid_element_set(), true, log);
    if (is_tu)
      return true;

    matroid_element_set rows, columns;
    for (std::size_t r = 0; r < matrix.size1(); ++r)
      rows.insert(-1 - r);
    for (std::size_t c = 0; c < matrix.size2(); ++c)
      columns.insert(1 + c);

    detail::violator_strategy* strategy = new detail::greedy_violator_strategy(matrix, rows, columns, log);

    strategy->search();
    strategy->create_matrix(violator);

    delete strategy;

    assert (violator.rows.size() == violator.columns.size());

    return false;
  }

  /**
   * Tests for total unimodularity with a negative certificate.
   * A matrix is totally unimodular if and only if every square submatrix
   * has a determinant of -1, 0 or +1.
   * If the matrix does not have this property,the routine returns the indices
   * of a violating submatrix.
   *
   * @param matrix The matrix to be tested
   * @param violator Returns violator indices
   * @param level Log level
   * @return true if and only if the matrix is totally unimodular
   */

  bool is_totally_unimodular(const integer_matrix& matrix, submatrix_indices& violator, log_level level)
  {
    decomposed_matroid* decomposition;

    bool is_tu = is_totally_unimodular(matrix, decomposition, violator, level);

    delete decomposition;

    return is_tu;
  }

  /**
   * Tests if a given matrix is a signed version of its support matrix already.
   * Running time: O(height * width * min(height, width))
   *
   * @param matrix The matrix to be tested
   * @return true if and only if the support matrix can be signed to the orignal
   */

  bool is_signed_matrix(const integer_matrix& matrix)
  {
    if (matrix.size2() > matrix.size1())
    {
      const matrix_transposed <const integer_matrix> transposed(matrix);
      return sign_matrix(transposed, NULL);
    }
    else
    {
      return sign_matrix(matrix, NULL);
    }
  }

  /**
   * Tests if a given matrix is a signed version of its support matrix already,
   * returning the indices of a violating submatrix if this is not the case.
   * Running time: O(height * width * min(height, width))
   *
   * @param matrix The matrix to be tested
   * @param violator Returns violator indices
   * @return true if and only if the support matrix can be signed to the original
   */

  bool is_signed_matrix(const integer_matrix& matrix, submatrix_indices& violator)
  {
    if (matrix.size2() > matrix.size1())
    {
      const matrix_transposed <const integer_matrix> transposed(matrix);
      bool result = sign_matrix(transposed, &violator);
      std::swap(violator.rows, violator.columns);
      return result;
    }
    else
    {
      return sign_matrix(matrix, &violator);
    }
  }

  /**
   * Signes a given matrix to be a signed version of its support matrix.
   * Running time: O(height * width * min(height, width))
   *
   * @param matrix The matrix to be signed
   * @return true if and only if any change was necessary
   */

  bool sign_matrix(integer_matrix& matrix)
  {
    if (matrix.size2() > matrix.size1())
    {
      matrix_transposed <integer_matrix> transposed(matrix);
      return sign_matrix(transposed, NULL);
    }
    else
    {
      return sign_matrix(matrix, NULL);
    }
  }

  /**
   * Makes the matrix its own support matrix.
   *
   * @param matrix The given matrix
   */

  void support_matrix(integer_matrix& matrix)
  {
    for (size_t i = 0; i < matrix.size1(); ++i)
    {
      for (size_t j = 0; j < matrix.size2(); ++j)
      {
        const int value = matrix(i, j);
        if (value != 0)
        {
          matrix(i, j) = 1;
        }
      }
    }
  }

}
