/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef ALGORITHM_HPP_
#define ALGORITHM_HPP_

#include "matroid_decomposition.hpp"
#include "find_wheel_minor.hpp"
#include "matroid_separate.hpp"
#include "nested_minor_sequence.hpp"
#include "find_minor_sequence.hpp"
#include "graphicness.hpp"
#include "r10.hpp"
#include "enumeration.hpp"
#include "logger.hpp"

namespace unimod
{

  /// Forward declaration

  template <typename MatroidType, typename MatrixType>
  std::pair<bool, decomposed_matroid*> decompose_binary_matroid(MatroidType& matroid, MatrixType& matrix,
      matroid_element_set extra_elements, bool construct_decomposition, logger& log);

  template <typename MatroidType, typename MatrixType>
  std::pair<bool, decomposed_matroid*> decompose_with_minor_sequence(matroid_permuted<MatroidType>& permuted_matroid,
      matrix_permuted<MatrixType>& permuted_matrix, nested_minor_sequence& nested_minors,
      matroid_element_set extra_elements, bool construct_decomposition, logger& log)
  {
    separation sep = find_minor_sequence(permuted_matroid, permuted_matrix, nested_minors, extra_elements, log);
    if (sep.is_valid())
    {
      if (log.is_progressive())
      {
        log.line() << " --> " << (sep.rank() + 1) << "-SEP";
        std::cout << log << std::endl;
        log.clear();
        log.indent();
      }
      else if (log.is_verbose())
      {
        std::cout << "Found a " << (sep.rank() + 1) << "-separation instead." << std::endl;
      }

      integer_matroid upper_left_matroid;
      integer_matrix upper_left_matrix;
      integer_matroid lower_right_matroid;
      integer_matrix lower_right_matrix;

      sep.create_components(permuted_matroid, permuted_matrix, upper_left_matroid, upper_left_matrix,
          lower_right_matroid, lower_right_matrix);

      if (log.is_verbose())
      {
        std::cout << "Summands are " << upper_left_matrix.size1() << " x " << upper_left_matrix.size2() << " and "
            << lower_right_matrix.size1() << " x " << lower_right_matrix.size2() << "." << std::endl;
      }

      /// Filter or copy extra_elements depending on type of separation
      matroid_element_set upper_left_extra_elements, lower_right_extra_elements;
      if (sep.rank() == 0)
      {
        matroid_element_set upper_left_elements = upper_left_matroid.get_elements();
        matroid_element_set lower_right_elements = lower_right_matroid.get_elements();

        std::set_difference(extra_elements.begin(), extra_elements.end(), lower_right_elements.begin(),
            lower_right_elements.end(), std::inserter(upper_left_extra_elements, upper_left_extra_elements.end()));
        std::set_difference(extra_elements.begin(), extra_elements.end(), upper_left_elements.begin(),
            upper_left_elements.end(), std::inserter(lower_right_extra_elements, lower_right_extra_elements.end()));
      }
      else
      {
        std::copy(extra_elements.begin(), extra_elements.end(), std::inserter(upper_left_extra_elements,
            upper_left_extra_elements.end()));
        std::copy(extra_elements.begin(), extra_elements.end(), std::inserter(lower_right_extra_elements,
            lower_right_extra_elements.end()));
      }

      matroid_permuted<integer_matroid> permuted_upper_left_matroid(upper_left_matroid);
      matrix_permuted<integer_matrix> permuted_upper_left_matrix(upper_left_matrix);

      if (sep.has_special_swap())
      {
        if (sep.has_special_row_swap())
          matroid_permute1(permuted_upper_left_matroid, permuted_upper_left_matrix, upper_left_matroid.size1() - 1,
              sep.get_special_swap_index());
        else
          matroid_permute2(permuted_upper_left_matroid, permuted_upper_left_matrix, upper_left_matroid.size2() - 1,
              sep.get_special_swap_index());
      }

      if (log.is_progressive())
      {
        log.line() << "(" << upper_left_matrix.size1() << " x " << upper_left_matrix.size2() << ") W3";
        std::cout << log;
      }
      else if (log.is_verbose())
      {
        std::cout << "Proceeding search for a sequence of nested minors in binary "
            << permuted_upper_left_matroid.size1() << " x " << permuted_upper_left_matrix.size2() << " matroid."
            << std::endl;
      }

      std::pair<bool, decomposed_matroid*> upper_left_result = decompose_with_minor_sequence(
          permuted_upper_left_matroid, permuted_upper_left_matrix, nested_minors, upper_left_extra_elements,
          construct_decomposition, log);

      if (!construct_decomposition && !upper_left_result.first)
      {
        log.unindent();
        return std::pair<bool, decomposed_matroid*>(false, NULL);
      }

      std::pair<bool, decomposed_matroid*> lower_right_result = decompose_binary_matroid(lower_right_matroid,
          lower_right_matrix, lower_right_extra_elements, construct_decomposition, log);

      if (log.is_progressive())
      {
        log.unindent();
      }

      if (construct_decomposition)
      {
        int type = sep.rank() == 0 ? decomposed_matroid_separator::ONE_SEPARATION
            : decomposed_matroid_separator::TWO_SEPARATION;

        return std::pair<bool, decomposed_matroid*>(lower_right_result.first && upper_left_result.first,
            new decomposed_matroid_separator(upper_left_result.second, lower_right_result.second, type,
                matroid_elements(permuted_matroid), extra_elements));
      }
      else
        return std::pair<bool, decomposed_matroid*>(lower_right_result.first, NULL);
    }

    if (log.is_verbose())
    {
      std::cout << "Constructed sequence of " << (nested_minors.size() + 1)
          << " 3-connected nested minors starting with W3." << std::endl;
    }

    size_t largest_graphic_minor = 0;
    matroid_graph* graph = construct_matroid_graph(permuted_matroid, permuted_matrix, nested_minors,
        largest_graphic_minor);

    if (log.is_progressive())
    {
      log.line() << (graph == NULL ? ", NON-GRAPHIC" : ", GRAPHIC");
      std::cout << log;
    }
    else if (log.is_verbose())
      std::cout << "Matroid is " << (graph == NULL ? "not " : "") << "graphic." << std::endl;

    if (!construct_decomposition && graph)
    {
      if (log.is_progressive())
      {
        log.clear();
        std::cout << " --> REGULAR" << std::endl;
      }
      else if (log.is_verbose())
      {
        std::cout << "Graphic matroids are always regular." << std::endl;
      }

      delete graph;
      return std::make_pair(true, (decomposed_matroid*) NULL);
    }

    size_t largest_cographic_minor = 0;
    matroid_graph* cograph = construct_matroid_graph(make_transposed_matroid(permuted_matroid), make_transposed_matrix(
        permuted_matrix), make_transposed_nested_minor_sequence(nested_minors), largest_cographic_minor);

    if (log.is_progressive())
    {
      log.line() << (cograph == NULL ? ", NON-COGRAPHIC" : ", COGRAPHIC");
      std::cout << log;
    }
    else if (log.is_verbose())
    {
      std::cout << "Matroid is " << (cograph == NULL ? "not " : "") << "cographic." << std::endl;
    }

    if (!construct_decomposition && cograph)
    {
      if (log.is_progressive())
      {
        log.clear();
        std::cout << " --> REGULAR" << std::endl;
      }
      else if (log.is_verbose())
      {
        std::cout << "Cographic matroids are always regular." << std::endl;
      }

      delete cograph;
      return std::make_pair(true, (decomposed_matroid*) NULL);
    }

    if (construct_decomposition && (graph || cograph))
    {
      if (log.is_progressive())
      {
        log.clear();
        std::cout << std::endl;
      }
      else if (log.is_verbose())
      {
        if (graph && cograph)
          std::cout << "Planar matroids are always regular." << std::endl;
        else if (graph && cograph == NULL)
          std::cout << "Planar matroids are always regular." << std::endl;
        else if (graph == NULL && cograph)
          std::cout << "Planar matroids are always regular." << std::endl;

      }

      return std::make_pair(true, new decomposed_matroid_leaf(graph, cograph, false,
          matroid_elements(permuted_matroid), extra_elements));
    }

    if (is_r10(permuted_matrix))
    {
      if (log.is_progressive())
      {
        log.line() << ", R10 --> REGULAR";
        std::cout << log << std::endl;
        log.clear();
      }
      else if (log.is_verbose())
      {
        std::cout << "Matroid is isomorphic to R10 and thus regular." << std::endl;
      }

      if (construct_decomposition)
      {
        return std::make_pair(true, new decomposed_matroid_leaf(NULL, NULL, true, matroid_elements(permuted_matroid),
            extra_elements));
      }
      else
        return std::make_pair(true, (decomposed_matroid*) NULL);
    }
    else
    {
      if (log.is_progressive())
      {
        log.line() << ", NOT R10";
        std::cout << log;
      }
      else if (log.is_verbose())
        std::cout << "Matroid is not isomorphic to R10." << std::endl;
    }

    size_t new_size = largest_graphic_minor > largest_cographic_minor ? largest_graphic_minor : largest_cographic_minor;
    if (log.is_progressive())
    {
      log.line() << ", (CO)GRAPHIC LEN: " << new_size;
      std::cout << log;
    }
    else if (log.is_verbose())
      std::cout << "Sequence is (co)graphic until N_" << new_size << "." << std::endl;

    nested_minors.resize(new_size);

    sep = enumerate_separations(permuted_matroid, permuted_matrix, nested_minors, extra_elements, log);
    if (sep.is_valid())
    {
      if (log.is_progressive())
      {
        log.line() << " --> " << (sep.rank() + 1) << "-SEP";
        std::cout << log << std::endl;
        log.clear();
        log.indent();
      }
      else if (log.is_verbose())
      {
        std::cout << "Found a (3|4)-separation." << std::endl;
      }

      integer_matroid upper_left_matroid;
      integer_matrix upper_left_matrix;
      integer_matroid lower_right_matroid;
      integer_matrix lower_right_matrix;

      sep.create_components(permuted_matroid, permuted_matrix, upper_left_matroid, upper_left_matrix,
          lower_right_matroid, lower_right_matrix);

      if (log.is_verbose())
      {
        std::cout << "Summands are " << upper_left_matrix.size1() << " x " << upper_left_matrix.size2() << " and "
            << lower_right_matrix.size1() << " x " << lower_right_matrix.size2() << "." << std::endl;
      }

      matroid_permuted<integer_matroid> permuted_upper_left_matroid(upper_left_matroid);
      matrix_permuted<integer_matrix> permuted_upper_left_matrix(upper_left_matrix);

      std::pair<bool, decomposed_matroid*> upper_left_result = decompose_binary_matroid(permuted_upper_left_matroid,
          permuted_upper_left_matrix, extra_elements, construct_decomposition, log);

      if (!construct_decomposition && !upper_left_result.first)
      {
        if (log.is_progressive())
        {
          log.unindent();
        }
        return std::pair<bool, decomposed_matroid*>(false, NULL);
      }

      std::pair<bool, decomposed_matroid*> lower_right_result = decompose_binary_matroid(lower_right_matroid,
          lower_right_matrix, extra_elements, construct_decomposition, log);

      if (log.is_progressive())
      {
        log.unindent();
      }

      if (construct_decomposition)
        return std::make_pair(lower_right_result.first && upper_left_result.first,
            (decomposed_matroid *) (new decomposed_matroid_separator(upper_left_result.second,
                lower_right_result.second, decomposed_matroid_separator::THREE_SEPARATION, matroid_elements(
                    permuted_matroid), extra_elements)));
      else
        return std::pair<bool, decomposed_matroid*>(lower_right_result.first, NULL);
    }

    if (construct_decomposition)
      return std::make_pair(false, new decomposed_matroid_leaf(NULL, NULL, false, matroid_elements(permuted_matroid),
          extra_elements));
    else
      return std::make_pair(false, (decomposed_matroid*) NULL);
  }

  /**
   * Create the graph for 2 x w or h x 2 matrices.
   *
   * @param matroid
   * @param matrix
   * @return Created graph
   */

  template <typename MatroidType, typename MatrixType>
  matroid_graph* construct_small_matroid_graph(MatroidType& matroid, MatrixType& matrix)
  {
    assert(matroid.size1() <= 2 || matroid.size2() <= 2);

    matroid_graph* graph = new matroid_graph(matroid.size1() + 1);

    if (matroid.size1() == 0)
    {
      for (size_t column = 0; column < matroid.size2(); ++column)
      {
        boost::add_edge(boost::vertex(0, *graph), boost::vertex(0, *graph), matroid.name2(0), *graph);
      }
    }
    else if (matroid.size2() == 0)
    {
      for (size_t row = 0; row < matroid.size1(); ++row)
      {
        boost::add_edge(boost::vertex(row, *graph), boost::vertex(row + 1, *graph), matroid.name1(row), *graph);
      }
    }
    else if (matroid.size1() >= 3 && matroid.size2() == 1)
    {
      size_t current_edge_vertex = 0;
      size_t current_free_vertex = matroid.size1();

      for (size_t row = 0; row < matroid.size1(); ++row)
      {
        if (matrix(row, 0) == 0)
        {
          boost::add_edge(boost::vertex(current_free_vertex - 1, *graph), boost::vertex(current_free_vertex, *graph),
              matroid.name1(row), *graph);
          --current_free_vertex;
        }
        else
        {
          boost::add_edge(boost::vertex(current_edge_vertex, *graph), boost::vertex(current_edge_vertex + 1, *graph),
              matroid.name1(row), *graph);
          ++current_edge_vertex;
        }
      }

      assert(current_edge_vertex == current_free_vertex);

      boost::add_edge(boost::vertex(0, *graph), boost::vertex(current_edge_vertex, *graph), matroid.name2(0), *graph);
    }
    else if (matroid.size1() >= 3 && matroid.size2() == 2)
    {
      size_t count[4] =
      { 0, 0, 0, 0 };
      for (size_t row = 0; row < matroid.size1(); ++row)
      {
        count[((matrix(row, 0) != 0) ? 2 : 0) + ((matrix(row, 1) != 0) ? 1 : 0)]++;
      }

      size_t vertex[4];
      vertex[0] = 0;
      vertex[1] = count[1];
      vertex[2] = vertex[1] + count[3];
      vertex[3] = vertex[2] + count[2];

      boost::add_edge(boost::vertex(vertex[0], *graph), boost::vertex(vertex[2], *graph), matroid.name2(1), *graph);
      boost::add_edge(boost::vertex(vertex[1], *graph), boost::vertex(vertex[3], *graph), matroid.name2(0), *graph);

      for (size_t row = 0; row < matroid.size1(); ++row)
      {
        int element = matroid.name1(row);
        switch (((matrix(row, 0) != 0) ? 2 : 0) + ((matrix(row, 1) != 0) ? 1 : 0))
        {
        case 1: /// 0 1
          boost::add_edge(boost::vertex(vertex[0], *graph), boost::vertex(vertex[0] + 1, *graph), element, *graph);
          ++vertex[0];
        break;
        case 2: /// 1 0
          boost::add_edge(boost::vertex(vertex[2], *graph), boost::vertex(vertex[2] + 1, *graph), element, *graph);
          ++vertex[2];
        break;
        case 3: /// 1 1
          boost::add_edge(boost::vertex(vertex[1], *graph), boost::vertex(vertex[1] + 1, *graph), element, *graph);
          ++vertex[1];
        break;
        default:
          assert(matrix(row, 0) == 0 && matrix(row, 1) == 0);
          boost::add_edge(boost::vertex(vertex[3], *graph), boost::vertex(vertex[3] + 1, *graph), element, *graph);
          ++vertex[3];
          break;
        }
      }
    }
    else
    {
      assert(matroid.size1() <= 2);

      for (size_t row = 0; row < matroid.size1(); ++row)
      {
        boost::add_edge(boost::vertex(row, *graph), boost::vertex(row + 1, *graph), matroid.name1(row), *graph);
      }

      for (size_t column = 0; column < matroid.size2(); ++column)
      {
        size_t end_vertex, start_vertex = 0;
        while (start_vertex < matroid.size1() && matrix(start_vertex, column) == 0)
        {
          ++start_vertex;
        }
        end_vertex = start_vertex;
        while (end_vertex < matroid.size1() && matrix(end_vertex, column) == 1)
        {
          ++end_vertex;
        }

        boost::add_edge(boost::vertex(start_vertex, *graph), boost::vertex(end_vertex, *graph), matroid.name2(column),
            *graph);
      }
    }

    return graph;
  }

  /**
   * Decomposes a given binary matroid.
   *
   * @param matroid
   * @param matrix
   * @param construct_decomposition
   * @return
   */

  template <typename MatroidType, typename MatrixType>
  std::pair<bool, decomposed_matroid*> decompose_binary_matroid(MatroidType& matroid, MatrixType& matrix,
      matroid_element_set extra_elements, bool construct_decomposition, logger& log)
  {
    assert(is_zero_one_matrix(matrix));

    if (log.is_progressive())
    {
      log.clear();
      log.line() << "(" << matrix.size1() << " x " << matrix.size2() << ")";
      std::cout << log;
    }
    else if (log.is_verbose())
    {
      std::cout << "Decomposing binary " << matrix.size1() << " x " << matrix.size2() << " matroid." << std::endl;
    }

    if (matroid.size1() <= 2 || matroid.size2() <= 2)
    {
      if (log.is_progressive())
      {
        log.line() << " TRIVIAL --> REGULAR";
        std::cout << log << std::endl;
        log.clear();
      }
      else if (log.is_verbose())
      {
        std::cout << "The matroid is trivial and thus regular." << std::endl;
      }

      if (construct_decomposition)
      {
        matroid_transposed<MatroidType> transposed_matroid(matroid);
        matrix_transposed<MatrixType> transposed_matrix(matrix);
        matroid_graph* g = construct_small_matroid_graph(matroid, matrix);
        matroid_graph* c = construct_small_matroid_graph(transposed_matroid, transposed_matrix);

        return std::make_pair(true, new decomposed_matroid_leaf(g, c, false, matroid_elements(matroid), extra_elements));
      }
      else
      {
        return std::make_pair(true, (decomposed_matroid*) NULL);
      }
    }

    typedef matroid_permuted<MatroidType> permuted_matroid_type;
    typedef matrix_permuted<MatrixType> permuted_marix_type;

    permuted_matroid_type permuted_matroid(matroid);
    permuted_marix_type permuted_matrix(matrix);

    if (log.is_progressive())
    {
      log.line() << " W3";
      std::cout << log;
    }
    else if (log.is_verbose())
    {
      std::cout << "Searching for a W3 minor." << std::endl;
    }

    /// Identifies a W_3 minor in the upper left corner or finds a separation

    separation sep = find_wheel_minor(permuted_matroid, permuted_matrix, extra_elements);

    if (sep.is_valid())
    {
      /// Separation case
      if (log.is_progressive())
      {
        log.line() << " --> " << (sep.rank() + 1) << "-separation";
        std::cout << log << std::endl;
        log.clear();
        log.indent();
      }
      else if (log.is_verbose())
      {
        std::cout << "Found a " << (sep.rank() + 1) << "-separation instead." << std::endl;
      }

      integer_matroid upper_left_matroid;
      integer_matrix upper_left_matrix;
      integer_matroid lower_right_matroid;
      integer_matrix lower_right_matrix;

      sep.create_components(permuted_matroid, permuted_matrix, upper_left_matroid, upper_left_matrix,
          lower_right_matroid, lower_right_matrix);

      /// Filter or copy extra_elements depending on type of separation
      matroid_element_set upper_left_extra_elements, lower_right_extra_elements;
      if (sep.rank() == 0 && false)
      {
        matroid_element_set upper_left_elements = upper_left_matroid.get_elements();
        matroid_element_set lower_right_elements = lower_right_matroid.get_elements();

        std::set_difference(extra_elements.begin(), extra_elements.end(), lower_right_elements.begin(),
            lower_right_elements.end(), std::inserter(upper_left_extra_elements, upper_left_extra_elements.end()));
        std::set_difference(extra_elements.begin(), extra_elements.end(), upper_left_elements.begin(),
            upper_left_elements.end(), std::inserter(lower_right_extra_elements, lower_right_extra_elements.end()));
      }
      else
      {
        std::copy(extra_elements.begin(), extra_elements.end(), std::inserter(upper_left_extra_elements,
            upper_left_extra_elements.end()));
        std::copy(extra_elements.begin(), extra_elements.end(), std::inserter(lower_right_extra_elements,
            lower_right_extra_elements.end()));
      }

      if (log.is_verbose())
      {
        std::cout << "Summands are " << upper_left_matrix.size1() << " x " << upper_left_matrix.size2() << " and "
            << lower_right_matrix.size1() << " x " << lower_right_matrix.size2() << "." << std::endl;
      }

      std::pair<bool, decomposed_matroid*> upper_left_result = decompose_binary_matroid(upper_left_matroid,
          upper_left_matrix, upper_left_extra_elements, construct_decomposition, log);

      if (!construct_decomposition && !upper_left_result.first)
      {
        if (log.is_progressive())
        {
          log.unindent();
        }
        return std::pair<bool, decomposed_matroid*>(false, NULL);
      }

      std::pair<bool, decomposed_matroid*> lower_right_result = decompose_binary_matroid(lower_right_matroid,
          lower_right_matrix, lower_right_extra_elements, construct_decomposition, log);

      if (log.is_progressive())
      {
        log.unindent();
      }

      if (construct_decomposition)
      {
        int type = sep.rank() == 0 ? decomposed_matroid_separator::ONE_SEPARATION
            : decomposed_matroid_separator::TWO_SEPARATION;

        return std::pair<bool, decomposed_matroid*>(lower_right_result.first && upper_left_result.first,
            new decomposed_matroid_separator(upper_left_result.second, lower_right_result.second, type,
                matroid_elements(permuted_matroid), extra_elements));
      }
      else
        return std::pair<bool, decomposed_matroid*>(lower_right_result.first && upper_left_result.first, NULL);
    }

    if (log.is_verbose())
    {
      std::cout << "Searching for a sequence of nested minors in binary " << permuted_matrix.size1() << " x "
          << permuted_matrix.size2() << " matroid." << std::endl;
    }

    nested_minor_sequence nested_minors;

    return decompose_with_minor_sequence(permuted_matroid, permuted_matrix, nested_minors, extra_elements,
        construct_decomposition, log);
  }

}

#endif /* ALGORITHM_HPP_ */
