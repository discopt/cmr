#include <iostream>

#include <polymake/Integer.h>
#include <polymake/Matrix.h>
#include <polymake/client.h>

#include "total_unimodularity.hpp"
#include "unimodularity.hpp"

namespace polymake
{
  namespace common
  {

    bool is_totally_unimodular_plain(const Matrix <Integer>& matrix)
    {
      unimod::integer_matrix input_matrix(matrix.rows(), matrix.cols());
      for (size_t r = 0; r < input_matrix.size1(); ++r)
      {
        for (size_t c = 0; c < input_matrix.size2(); ++c)
        {
          input_matrix(r, c) = matrix(r, c);
        }
      }

      return unimod::is_totally_unimodular(input_matrix);
    }

    bool is_totally_unimodular_violator(const Matrix <Integer>& matrix, Vector <Integer>& rows, Vector <Integer>& columns)
    {
      unimod::integer_matrix input_matrix(matrix.rows(), matrix.cols());
      for (size_t r = 0; r < input_matrix.size1(); ++r)
      {
        for (size_t c = 0; c < input_matrix.size2(); ++c)
          input_matrix(r, c) = matrix(r, c);
      }

      unimod::submatrix_indices indices;

      bool result = unimod::is_totally_unimodular(input_matrix, indices);
      if (!result)
      {
        rows.resize(indices.rows.size());
        for (size_t i = 0; i < indices.rows.size(); ++i)
          rows[i] = (int) indices.rows[i];
        columns.resize(indices.columns.size());
        for (size_t i = 0; i < indices.columns.size(); ++i)
          columns[i] = (int) indices.columns[i];
      }
      return result;
    }

    bool is_unimodular_plain(const Matrix <Integer>& matrix)
    {
      size_t rank;
      unimod::integer_matrix input_matrix(matrix.rows(), matrix.cols());
      for (size_t r = 0; r < input_matrix.size1(); ++r)
      {
        for (size_t c = 0; c < input_matrix.size2(); ++c)
        {
          input_matrix(r, c) = matrix(r, c);
        }
      }

      return unimod::is_unimodular(input_matrix, rank);
    }

    bool is_strongly_unimodular_plain(const Matrix <Integer>& matrix)
    {
      size_t rank;
      unimod::integer_matrix input_matrix(matrix.rows(), matrix.cols());
      for (size_t r = 0; r < input_matrix.size1(); ++r)
      {
        for (size_t c = 0; c < input_matrix.size2(); ++c)
        {
          input_matrix(r, c) = matrix(r, c);
        }
      }

      return unimod::is_strongly_unimodular(input_matrix, rank);
    }

    perl::ListReturn is_k_modular_compute(const Matrix <Integer>& matrix)
    {
      size_t rank;
      unsigned int k;
      unimod::integer_matrix input_matrix(matrix.rows(), matrix.cols());
      for (size_t r = 0; r < input_matrix.size1(); ++r)
      {
        for (size_t c = 0; c < input_matrix.size2(); ++c)
        {
          input_matrix(r, c) = matrix(r, c);
        }
      }

      perl::ListReturn result;
      if (unimod::is_k_modular(input_matrix, rank, k))
        result << (int) k;
      return result;
    }

    perl::ListReturn is_strongly_k_modular_compute(const Matrix <Integer>& matrix)
    {
      size_t rank;
      unsigned int k;
      unimod::integer_matrix input_matrix(matrix.rows(), matrix.cols());
      for (size_t r = 0; r < input_matrix.size1(); ++r)
      {
        for (size_t c = 0; c < input_matrix.size2(); ++c)
        {
          input_matrix(r, c) = matrix(r, c);
        }
      }

      if (unimod::is_strongly_k_modular(input_matrix, rank, k))
      {

      }
      perl::ListReturn result;
      if (unimod::is_strongly_k_modular(input_matrix, rank, k))
        result << (int) k;
      return result;
    }

  /// total unimodularity

  UserFunction4perl  ("# Tests a given //matrix// for total unimodularity, without certificates."
      "# @param Matrix<Integer> matrix"
      "# @return Bool",
      &is_totally_unimodular_plain,"is_totally_unimodular(Matrix<Integer>)");

  /// total unimodularity with violator

  UserFunction4perl("# Tests a given //matrix// for total unimodularity, setting row/column indices to a violating submatrix if the result is negative."
      "# @param Matrix<Integer> matrix"
      "# @param Vector<Integer> row_indices of violating submatrix."
      "# @param Vector<Integer> column_indices of violating submatrix."
      "# @return Bool",
      &is_totally_unimodular_violator,"is_totally_unimodular(Matrix<Integer>,Vector<Integer>,Vector<Integer>)");

  /// unimodularity

  UserFunction4perl("# Tests a given //matrix// for unimodularity, without certificates."
      "# @param Matrix<Integer> matrix"
      "# @return Bool",
      &is_unimodular_plain,"is_unimodular(Matrix<Integer>)");

  /// strong unimodularity

  UserFunction4perl("# Tests a given //matrix// for strong unimodularity, without certificates."
      "# @param Matrix<Integer> matrix"
      "# @return Bool",
      &is_strongly_unimodular_plain,"is_strongly_unimodular(Matrix<Integer>)");

  /// k-modularity

  UserFunction4perl("# Tests a given //matrix// for k-modularity, without certificates. It also computes k."
      "# @param Matrix<Integer> matrix"
      "# @return ListReturn k or undefined",
      &is_k_modular_compute,"is_k_modular(Matrix<Integer>)");

  /// strong k-modularity

  UserFunction4perl("# Tests a given //matrix// for strong k-modularity, without certificates."
      "# @param Matrix<Integer> matrix"
      "# @return ListReturn k or undefined",
      &is_strongly_k_modular_compute,"is_strongly_k_modular(Matrix<Integer>)");
}}

