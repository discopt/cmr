/*
 * unimodularity.cpp
 *
 *  Created on: Feb 3, 2011
 *      Author: xammy
 */

#include <tu/unimodularity.hpp>

#include <tu/common.hpp>
#include <tu/linear_algebra.hpp>
#include <tu/smith_normal_form.hpp>
#include <tu/total_unimodularity.hpp>

namespace tu
{

  template <typename Matrix>
  bool test_k_modularity(const Matrix& matrix, size_t& rank, unsigned int* pk, bool enforce_unimodularity, log_level level)
  {
    typedef boost::numeric::ublas::matrix_indirect <const integer_matrix, submatrix_indices::indirect_array_type> indirect_matrix_t;

    integer_matrix transformed;
    std::vector <size_t> basis;

    rank = matrix_find_column_basis_and_transform_integral(matrix, transformed, basis);
    assert(rank == transformed.size1());
    assert(rank == basis.size());

    bool found_identity_matrix = true;
    for (size_t i = 0; i < basis.size(); ++i)
    {
      size_t c = basis[i];
      for (size_t r = 0; r < transformed.size1(); ++r)
      {
        if ((r == i && transformed(r, c) != 1) || ((r != i) && transformed(r, c) != 0))
        {
          found_identity_matrix = false;
          break;
        }
      }
      if (!found_identity_matrix)
        break;
    }

    if (!found_identity_matrix)
    {
      if (pk)
        *pk = 0;
      return false;
    }

    if (pk)
    {
      /// In case we need k, let's compute Smith Normal Form of submatrix under column basis.
      submatrix_indices basis_matrix_indices =
      { submatrix_indices::indirect_array_type(matrix.size1()), submatrix_indices::indirect_array_type(rank) };

      for (size_t i = 0; i < rank; ++i)
        basis_matrix_indices.columns[i] = basis[i];
      for (size_t r = 0; r < matrix.size1(); ++r)
        basis_matrix_indices.rows[r] = r;

      const indirect_matrix_t basis_matrix(matrix, basis_matrix_indices.rows, basis_matrix_indices.columns);

      /// Assert that Basis * Transformed = Original
      assert(equals(boost::numeric::ublas::prod(basis_matrix, transformed), matrix));

      /// Compute Smith Normal Form and check that the gcd of all full-rank submatrices is 1.
      std::vector <int> smith_diagonal;
      smith_normal_form_diagonal(basis_matrix, smith_diagonal);
      *pk = 1;
      for (size_t i = 0; i < smith_diagonal.size(); ++i)
      {
        *pk *= smith_diagonal[i];
      }

      if (enforce_unimodularity && (*pk != 1))
        return false;
    }

    /// Create submatrix of non-basis.
    submatrix_indices nonbasis_matrix_indices =
    { submatrix_indices::indirect_array_type(rank), submatrix_indices::indirect_array_type(transformed.size2() - rank) };

    size_t c = 0;
    for (size_t j = 0; j < transformed.size2(); ++j)
    {
      if (std::find(basis.begin(), basis.end(), j) == basis.end())
      {
        nonbasis_matrix_indices.columns[c++] = j;
      }
    }
    for (size_t r = 0; r < rank; ++r)
      nonbasis_matrix_indices.rows[r] = r;

    const indirect_matrix_t nonbasis_transformed(transformed, nonbasis_matrix_indices.rows, nonbasis_matrix_indices.columns);

    return is_totally_unimodular(nonbasis_transformed, LOG_QUIET);
  }

  /**
   * Tests for unimodularity without certificates.
   * A matrix of rank r is unimodular if and only if for every
   * column basis B the gcd (greatest common divisor) of the
   * determinants of regular r x r submatrices of B is 1.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank k
   * @param level Log level
   * @return true if and only if the matrix is unimodular
   */

  bool is_unimodular(const integer_matrix& matrix, size_t& rank, log_level level)
  {
    unsigned int k;
    bool result = test_k_modularity(matrix, rank, &k, true, level);
    return result && k == 1;
  }

  /**
   * Tests for strong unimodularity without certificates.
   * A matrix is strongly unimodular if and only if
   * it and its transpose are unimodular.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank r of the matrix
   * @param level Log level
   * @return true if and only if the matrix is strongly unimodular
   */

  bool is_strongly_unimodular(const integer_matrix& matrix, size_t& rank, log_level level)
  {
    size_t rank2;
    unsigned int k1, k2;
    bool result;

    result = test_k_modularity(matrix, rank, &k1, true, level);
    if (!result || k1 != 1)
      return false;

    const matrix_transposed <const integer_matrix> transposed(matrix);
    result = test_k_modularity(transposed, rank2, &k2, true, level);
    assert(k1 == k2);
    assert(rank == rank2);
    return result;
  }

  /**
   * Tests for k-modularity without certificates.
   * A matrix of rank r is k-modular if and only if for every
   * column basis B the gcd (greatest common divisor) of the
   * determinants of regular r x r submatrices of B has the same
   * absolute value.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank r
   * @param level Log level
   * @return true if and only if the matrix is k-modular
   */

  bool is_k_modular(const integer_matrix& matrix, size_t& rank, log_level level)
  {
    try
    {
      bool result = test_k_modularity(matrix, rank, NULL, false, level);
      return result;
    }
    catch (std::runtime_error &e)
    {
      std::cerr << "[Error: " << e.what() << " - Check for k-modularity failed.]\n" << std::flush;
      return false;
    }
  }

  /**
   * Tests for k-modularity without certificates.
   * A matrix of rank r is k-modular if and only if for every
   * column basis B the gcd (greatest common divisor) of the
   * determinants of regular r x r submatrices of B has the same
   * absolute value k.
   * It also computes k.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank r
   * @param k The common absolute value
   * @param level Log level
   * @return true if and only if the matrix is k-modular
   */

  bool is_k_modular(const integer_matrix& matrix, size_t& rank, unsigned int& k, log_level level)
  {
    bool result;
    try
    {
      result = test_k_modularity(matrix, rank, &k, false, level);
      return result;
    }
    catch (std::runtime_error &e)
    {
      std::cerr << "[Error: " << e.what() << " - Check for k-modularity failed.]\n" << std::flush;
      result = false;
      k = 0;
    }
    if (!result)
      k = 0;
    return result;
  }

  /**
   * Tests for strong k-modularity without certificates.
   * A matrix of rank r is k-modular if and only if
   * it and its transpose are k-modular.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank r
   * @param level Log level
   * @return true if and only if the matrix is k-modular
   */

  bool is_strongly_k_modular(const integer_matrix& matrix, size_t& rank, log_level level)
  {
    size_t rank2;
    bool result;

    if (!test_k_modularity(matrix, rank, NULL, false, level))
      return false;

    const matrix_transposed <const integer_matrix> transposed(matrix);
    result = test_k_modularity(transposed, rank2, NULL, false, level);
    assert(rank == rank2);
    return result;
  }

  /**
   * Tests for strong k-modularity without certificates.
   * A matrix of rank r is k-modular if and only if
   * it and its transpose are k-modular.
   * It also computes k.
   *
   * @param matrix The matrix to be tested
   * @param rank Returns the rank r
   * @param k The common absolute value
   * @param level Log level
   * @return true if and only if the matrix is k-modular
   */

  bool is_strongly_k_modular(const integer_matrix& matrix, size_t& rank, unsigned int& k, log_level level)
  {
    size_t rank2;
    unsigned int k1, k2;
    bool result;

    result = test_k_modularity(matrix, rank, &k1, false, level);
    if (!result)
    {
      k = 0;
      return false;
    }

    const matrix_transposed <const integer_matrix> transposed(matrix);
    result = test_k_modularity(transposed, rank2, &k2, false, level);
    assert(k1 == k2);
    assert(rank == rank2);
    k = result ? k1 : 0;
    return result;
  }

  /**
   * In case a matrix A is k-modular, it may lead to q-integrality of the
   * polyhedron A*x = b, x >= 0 if B*x = b*q for some basis B of A. This method
   * finds a minimal integer q, assuming that A is k-modular.
   * If q = 1, the polyhedron is integral, if q = 2 it is half-integral, etc.
   *
   * @param matrix The matrix A, being k-modular.
   * @param rhs The rhs, given as a m x 1 matrix.
   * @return Minimal q as defined above
   */

  unsigned int get_k_modular_integrality(const integer_matrix& matrix, const integer_matrix& rhs)
  {
    // TODO:
    return 0;
  }

  /**
   * In case a matrix A is k-modular, it may lead to integrality of the
   * polyhedron A*x = b, x >= 0 if B*x = b for some basis B of A.
   * This method tests of that property.
   *
   * @param matrix The matrix A, having the Dantzig property.
   * @param rhs The rhs, given as a m x 1 matrix.
   * @return true if and only if the polyhedron is integral
   */

  bool is_k_modular_integral(const integer_matrix& matrix, const integer_matrix& rhs)
  {
    return get_k_modular_integrality(matrix, rhs) == 1;
  }

  bool is_complement_total_unimodular(const integer_matrix& matrix, std::size_t& complementedRow, std::size_t& complementedColumn, log_level level)
  {
    integer_matrix row_complemented(matrix.size1(), matrix.size2());
    integer_matrix column_complemented(matrix.size1(), matrix.size2());
    std::pair <integer_matrix::size_type, integer_matrix::size_type> position;
    if (!is_zero_one_matrix(matrix, position))
    {
      complementedRow = matrix.size1();
      complementedColumn = matrix.size2();
      return false;
    }

    for (std::size_t crow = 0; crow <= matrix.size1(); ++crow)
    {
      if (crow != matrix.size1())
      {
        // Complement crow.

        for (std::size_t r = 0; r < matrix.size1(); ++r)
        {
          for (std::size_t c = 0; c < matrix.size2(); ++c)
          {
            if (r != crow && matrix(crow,c) != 0)
              row_complemented(r,c) = 1 - matrix(r,c);
            else
              row_complemented(r,c) = matrix(r,c);
          }
        }
      }

      for (std::size_t ccolumn = 0; ccolumn <= matrix.size2(); ++ccolumn)
      {
        if (ccolumn != matrix.size2())
        {
          // Complement ccolumn.

          for (std::size_t r = 0; r < matrix.size1(); ++r)
          {
            for (std::size_t c = 0; c < matrix.size2(); ++c)
            {
              if (c != ccolumn && row_complemented(r,ccolumn) != 0)
                column_complemented(r,c) = 1 - row_complemented(r,c);
              else
                column_complemented(r,c) = row_complemented(r,c);
            }
          }
        }

        if (crow == matrix.size1() && ccolumn == matrix.size2())
          continue;

        if (!is_totally_unimodular(column_complemented))
        {
          complementedRow = crow;
          complementedColumn = ccolumn;
          return false;
        }
      }
    }

    complementedRow = matrix.size1();
    complementedColumn = matrix.size2();
    return true;
  }
}
