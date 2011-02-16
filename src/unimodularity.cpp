/*
 * unimodularity.cpp
 *
 *  Created on: Feb 3, 2011
 *      Author: xammy
 */

#include "unimodularity.hpp"

#include "common.hpp"
#include "linear_algebra.hpp"
#include "smith_normal_form.hpp"
#include "total_unimodularity.hpp"

namespace unimod
{

  bool is_unimodular(const integer_matrix& matrix, log_level level)
  {
    size_t rank;
    return is_unimodular(matrix, rank, level);
  }

  bool is_unimodular(const integer_matrix& matrix, size_t& rank, log_level level)
  {
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
      return false;

    /// Create submatrix for column basis.
    submatrix_indices basis_matrix_indices =
    { submatrix_indices::indirect_array_type(matrix.size1()), submatrix_indices::indirect_array_type(rank) };

    for (size_t i = 0; i < rank; ++i)
      basis_matrix_indices.columns[i] = basis[i];
    for (size_t r = 0; r < matrix.size1(); ++r)
      basis_matrix_indices.rows[r] = r;

    typedef boost::numeric::ublas::matrix_indirect <const integer_matrix, submatrix_indices::indirect_array_type> indirect_matrix_t;
    const indirect_matrix_t basis_matrix(matrix, basis_matrix_indices.rows, basis_matrix_indices.columns);

    /// Assert that Basis * Transformed = Original
    assert(equals(boost::numeric::ublas::prod(basis_matrix, transformed), matrix));

    /// Compute Smith Normal Form and check that the gcd of all full-rank submatrices is 1.
    std::vector <int> smith_diagonal;
    smith_normal_form_diagonal(basis_matrix, smith_diagonal);
    for (size_t i = 0; i < smith_diagonal.size(); ++i)
    {
      if (smith_diagonal[i] != 1)
        return false;
    }

    /// Create submatrix of non-basis.
    submatrix_indices nonbasis_matrix_indices =
    { submatrix_indices::indirect_array_type(rank), submatrix_indices::indirect_array_type(transformed.size2() - rank) };

    size_t i = 0;
    size_t c = 0;
    for (size_t j = 0; j < transformed.size2(); ++j)
    {
      if (basis[i] == j)
        i++;
      else
        nonbasis_matrix_indices.columns[c++] = basis[i];
    }
    for (size_t r = 0; r < rank; ++r)
      nonbasis_matrix_indices.rows[r] = r;

    const indirect_matrix_t nonbasis_transformed(transformed, nonbasis_matrix_indices.rows, nonbasis_matrix_indices.columns);

    return is_totally_unimodular(nonbasis_transformed, level);
  }

  bool has_dantzig_property(const integer_matrix& matrix, log_level level)
  {
    size_t rank;
    return has_dantzig_property(matrix, rank, level);
  }

  bool has_dantzig_property(const integer_matrix& matrix, size_t& rank, log_level level)
  {
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
      return false;

    /// Create submatrix for column basis.
    submatrix_indices basis_matrix_indices =
    { submatrix_indices::indirect_array_type(matrix.size1()), submatrix_indices::indirect_array_type(rank) };

    for (size_t i = 0; i < rank; ++i)
      basis_matrix_indices.columns[i] = basis[i];
    for (size_t r = 0; r < matrix.size1(); ++r)
      basis_matrix_indices.rows[r] = r;

    typedef boost::numeric::ublas::matrix_indirect <const integer_matrix, submatrix_indices::indirect_array_type> indirect_matrix_t;
    const indirect_matrix_t basis_matrix(matrix, basis_matrix_indices.rows, basis_matrix_indices.columns);

    /// Assert that Basis * Transformed = Original
    assert(equals(boost::numeric::ublas::prod(basis_matrix, transformed), matrix));

    /// Create submatrix of non-basis.
    submatrix_indices nonbasis_matrix_indices =
    { submatrix_indices::indirect_array_type(rank), submatrix_indices::indirect_array_type(transformed.size2() - rank) };

    size_t i = 0;
    size_t c = 0;
    for (size_t j = 0; j < transformed.size2(); ++j)
    {
      if (basis[i] == j)
        i++;
      else
        nonbasis_matrix_indices.columns[c++] = basis[i];
    }
    for (size_t r = 0; r < rank; ++r)
      nonbasis_matrix_indices.rows[r] = r;

    const indirect_matrix_t nonbasis_transformed(transformed, nonbasis_matrix_indices.rows, nonbasis_matrix_indices.columns);

    return is_totally_unimodular(nonbasis_transformed, level);
  }

  unsigned int get_dantzig_integrality(const integer_matrix& matrix, const integer_matrix& rhs)
  {
    return 0;
  }
}
