/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef FIND_MINOR_SEQUENCE_HPP_
#define FIND_MINOR_SEQUENCE_HPP_

#include "../config.h"
#include "nested_minor_sequence.hpp"
#include "vector_three_connectivity.hpp"
#include "separation.hpp"
#include "matroid_transposed.hpp"
#include "matrix_transposed.hpp"
#include "matrix_modified.hpp"
#include "bipartite_graph_bfs.hpp"
#include "comparators.hpp"
#include "logger.hpp"

namespace tu {

  template <typename MatroidType, typename MatrixType, typename NestedMinorSequenceType, typename RowThreeConnectivity,
      typename ColumnThreeConnectivity>
  bool find_simple_row_extension (MatroidType& matroid, MatrixType& matrix, NestedMinorSequenceType& nested_minors,
      RowThreeConnectivity& row_three_connectivity, ColumnThreeConnectivity& column_three_connectivity)
  {
    for (size_t row = nested_minors.height(); row < matroid.size1(); ++row)
    {
      if (row_three_connectivity.is_other(row)) // Neither parallel, nor zero- / unit-vector
      {
        matroid_permute1(matroid, matrix, nested_minors.height(), row);
        row_three_connectivity.swap_vectors(nested_minors.height(), row);
        nested_minors.push(nested_minor_sequence::ONE_ROW);
        row_three_connectivity.enlarge_base();
        column_three_connectivity.enlarge_dimension();
        return true;
      }
    }

    return false;
  }

  template <typename MatroidType, typename MatrixType, typename NestedMinorSequenceType, typename RowThreeConnectivity,
      typename ColumnThreeConnectivity>
  char find_parallel_or_unit_vector (MatroidType& matroid, MatrixType& matrix, NestedMinorSequenceType& nested_minors,
      RowThreeConnectivity& row_three_connectivity, ColumnThreeConnectivity& column_three_connectivity, size_t& index)
  {
    bool found_column = false;
    for (size_t row = nested_minors.height(); row < matroid.size1(); ++row)
    {
      //      std::cout << "row = " << row << std::endl;
      if (row_three_connectivity.is_parallel(row))
      {
        //        std::cout << "parallel" << std::endl;
        index = row_three_connectivity.get_referred(row);
        return 'r';
      }
      else if (row_three_connectivity.is_unit(row))
      {
        //        std::cout << "unit" << std::endl;
        index = row_three_connectivity.get_referred(row);
        found_column = true;
      }
    }

    for (size_t column = nested_minors.width(); column < matroid.size2(); ++column)
    {
      //      std::cout << "column = " << column << std::endl;
      if (column_three_connectivity.is_unit(column))
      {
        //        std::cout << "unit" << std::endl;
        index = column_three_connectivity.get_referred(column);
        return 'r';
      }
      else if (column_three_connectivity.is_parallel(column))
      {
        //        std::cout << "parallel" << std::endl;
        index = column_three_connectivity.get_referred(column);
        found_column = true;
      }
    }
    return found_column ? 'c' : 0;
  }

  struct elaborate_extension_matrix_modifier
  {
    typedef int value_type;
    typedef unsigned char index_type;

    static const index_type BLOCK = 0;
    static const index_type ZERO = 1;
    static const index_type START = 2;
    static const index_type END0 = 3;
    static const index_type END1 = 4;

    elaborate_extension_matrix_modifier (const std::vector <index_type>& row_types, const std::vector <index_type>& column_types) :
      row_types_(row_types), column_types_(column_types)
    {

    }

    int operator () (size_t i, size_t j, int value)
    {
      switch (5 * row_types_[i] + column_types_[j])
      {
      case 5 * ZERO + ZERO:
      case 5 * ZERO + START:
      case 5 * ZERO + END0:
      case 5 * ZERO + END1:
      case 5 * START + ZERO:
      case 5 * END0 + ZERO:
      case 5 * END0 + START:
      case 5 * END0 + END0:
      case 5 * END0 + END1:
      case 5 * END1 + ZERO:
      case 5 * END1 + START:
      case 5 * END1 + END0:
      case 5 * END1 + END1:
        return value;
      case 5 * START + END0:
        return value;
      case 5 * START + END1:
        return 1 - value;
      default:
        return 0;
      }
    }

  private:
    const std::vector <index_type>& row_types_;
    const std::vector <index_type>& column_types_;
  };

  template <typename MatroidType, typename MatrixType, typename NestedMinorSequenceType, typename RowThreeConnectivity,
      typename ColumnThreeConnectivity>
  separation find_elaborate_extension (MatroidType& matroid, MatrixType& matrix, NestedMinorSequenceType& nested_minors,
      RowThreeConnectivity& row_three_connectivity, ColumnThreeConnectivity& column_three_connectivity, size_t index,
      matroid_element_set& extra_elements)
  {
    //    std::cout << "find_elaborate_extension on matroid with already found sequence of size " << nested_minors.height () << ","
    //        << nested_minors.width () << " on index " << index << "\n" << std::endl;
    //    matroid_print (matroid, matrix);

    // initialize bfs
    bipartite_graph_dimensions dim(matroid.size1(), matroid.size2());
    std::vector <size_t> start_nodes, end_nodes;
    start_nodes.reserve(matroid.size1());
    end_nodes.reserve(matroid.size1());

    // bfs-matrix 
    std::vector <elaborate_extension_matrix_modifier::index_type> row_types(matroid.size1());
    for (size_t row = 0; row < matroid.size1(); ++row)
    {
      if (row < nested_minors.height())
        row_types[row] = elaborate_extension_matrix_modifier::BLOCK;
      else if (row_three_connectivity.is_parallel(row) && row_three_connectivity.get_referred(row) == index)
      {
        row_types[row] = elaborate_extension_matrix_modifier::START;
        start_nodes.push_back(dim.row_to_index(row));
      }
      else if (row_three_connectivity.is_zero(row))
        row_types[row] = elaborate_extension_matrix_modifier::ZERO;
      else
      {
        row_types[row] = elaborate_extension_matrix_modifier::END0;
        end_nodes.push_back(dim.row_to_index(row));
      }
      //        std::cout << "row " << row << " is of type " << int(row_types[row]) << std::endl;
    }

    std::vector <elaborate_extension_matrix_modifier::index_type> column_types(matroid.size2());
    for (size_t column = 0; column < matroid.size2(); ++column)
    {
      if (column < nested_minors.width())
        column_types[column] = elaborate_extension_matrix_modifier::BLOCK;
      else if (column_three_connectivity.is_unit(column) && column_three_connectivity.get_referred(column) == index)
      {
        column_types[column] = elaborate_extension_matrix_modifier::START;
        start_nodes.push_back(dim.column_to_index(column));
      }
      else if (column_three_connectivity.is_zero(column))
        column_types[column] = elaborate_extension_matrix_modifier::ZERO;
      else
      {
        column_types[column] = matrix(index, column) == 1 ? elaborate_extension_matrix_modifier::END1 : elaborate_extension_matrix_modifier::END0;
        end_nodes.push_back(dim.column_to_index(column));
      }

      //        std::cout << "column " << column << " is of type " << int(column_types[column]) << std::endl;
    }

    elaborate_extension_matrix_modifier modifier(row_types, column_types);
    matrix_modified <MatrixType, elaborate_extension_matrix_modifier> modified_matrix(matrix, modifier);

    //    matroid_print (matroid, modified_matrix);

    std::vector <bipartite_graph_bfs_node> bfs_result;
    bool found_path = bipartite_graph_bfs(modified_matrix, dim, start_nodes, end_nodes, false, bfs_result);

    //    for (size_t i = 0; i < bfs_result.size (); ++i)
    //    {
    //      std::cout << "node " << i << " is ";
    //      if (dim.is_row (i))
    //        std::cout << "row " << dim.index_to_row (i);
    //      else
    //        std::cout << "column " << dim.index_to_column (i);
    //      std::cout << " with distance " << bfs_result[i].distance << " and pred " << bfs_result[i].predecessor << std::endl;
    //    }
    //    std::cout << "found_path = " << int(found_path) << std::endl;

    if (found_path)
    {
      size_t nearest_end = 0;
      for (std::vector <bipartite_graph_dimensions::index_type>::const_iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter)
      {
        if (bfs_result[*iter].is_reachable())
          nearest_end = *iter;
      }
      assert (bfs_result[nearest_end].is_reachable());

      size_t nearest_distance = bfs_result[nearest_end].distance;

      //        std::cout << "Before pivots:\n";
      //        matroid_print (matroid, matrix);
      //        std::cout << "nearest end = " << nearest_end << " with a path of length " << nearest_distance << std::endl;

      size_t connecting_start = 0;
      size_t first_connected = 0;
      size_t last_index = nearest_end;
      size_t current_index = bfs_result[nearest_end].predecessor;

      size_t count = 0;
      while (last_index != current_index)
      {
        std::pair <size_t, size_t> coords = dim.indexes_to_coordinates(current_index, last_index);

        //            std::cout << "I am at " << coords.first << "," << coords.second << std::endl;
        if ((row_types[coords.first] == elaborate_extension_matrix_modifier::ZERO) && (column_types[coords.second]
            == elaborate_extension_matrix_modifier::ZERO))
        {
          //                std::cout << "SWAP ALLOWED" << std::endl;
          //                std::cout << "row type is " << int(row_types[coords.first]) << ", column type is "
          //                        << int(column_types[coords.second]) << std::endl;

          //            if ((count >= 1) && (bfs_result[current_index].distance >= 2))
          //            {
          //                std::cout << "Would pivot here" << std::endl;
          if (count % 2 == 0)
          {
            matroid_binary_pivot(matroid, matrix, coords.first, coords.second);
            extra_elements.insert(matroid.name1(coords.first));
            extra_elements.insert(matroid.name2(coords.second));
            nearest_distance -= 2;
          }
          count++;
        }
        //            ++count;

        if (bfs_result[current_index].distance == 0)
        {
          first_connected = last_index;
          connecting_start = current_index;
        }

        last_index = current_index;
        current_index = bfs_result[current_index].predecessor;
      }

      //        std::cout << "connecting start = " << connecting_start << std::endl;
      //        std::cout << "first connected = " << first_connected << std::endl;
      //        std::cout << "After pivots:\n";
      //        matroid_print (matroid, matrix);

      if (nearest_distance == 1)
      {
        if (dim.is_column(nearest_end))
          std::swap(nearest_end, connecting_start);

        assert (dim.is_row(nearest_end));
        assert (dim.is_column(connecting_start));

        matroid_permute1(matroid, matrix, dim.index_to_row(nearest_end), nested_minors.height());
        matroid_permute2(matroid, matrix, dim.index_to_column(connecting_start), nested_minors.width());
        nested_minors.push(nested_minor_sequence::ONE_ROW_ONE_COLUMN);
        return separation();
      }
      else if (nearest_distance == 2)
      {
        if (connecting_start > nearest_end)
          std::swap(connecting_start, nearest_end);

        if (dim.is_row(nearest_end))
        {
          assert (dim.is_row(connecting_start));

          matroid_permute1(matroid, matrix, dim.index_to_row(connecting_start), nested_minors.height());
          matroid_permute1(matroid, matrix, dim.index_to_row(nearest_end), nested_minors.height() + 1);
          matroid_permute2(matroid, matrix, nested_minors.width(), dim.index_to_column(first_connected));
          nested_minors.push(nested_minor_sequence::TWO_ROWS_ONE_COLUMN);
          return separation();
        }
        else
        {
          assert (dim.is_column(connecting_start));
          assert (dim.is_column(nearest_end));

          matroid_permute2(matroid, matrix, dim.index_to_column(connecting_start), nested_minors.width());
          matroid_permute2(matroid, matrix, dim.index_to_column(nearest_end), nested_minors.width() + 1);
          matroid_permute1(matroid, matrix, nested_minors.height(), dim.index_to_row(first_connected));
          nested_minors.push(nested_minor_sequence::ONE_ROW_TWO_COLUMNS);
          return separation();
        }
      }

      return separation();
    }
    else
    {
      // this implies a 2-separation

      //        std::cout << "NO PATH FOUND!!!" << std::endl;

      std::pair <size_t, size_t> split(0, 0);

      std::vector <int> reachable(matrix.size1());
      for (size_t row = 0; row < matrix.size1(); ++row)
      {
        const bipartite_graph_bfs_node& node = bfs_result[dim.row_to_index(row)];
        int value = row < nested_minors.height() ? row - nested_minors.height() : (node.is_reachable() ? 2 : 1);
        reachable[row] = value;
        if (value < 2)
          split.first++;
      }
      vector_less <int, std::less <int> > less(reachable, std::less <int>());

      permutation p(matrix.size1());
      sort(p, less);
      matroid_apply_row_permutation(matroid, matrix, p);

      reachable.resize(matrix.size2());
      for (size_t column = 0; column < matrix.size2(); ++column)
      {
        const bipartite_graph_bfs_node& node = bfs_result[dim.column_to_index(column)];
        int value = column < nested_minors.width() ? column - nested_minors.width() : (node.is_reachable() ? 2 : 1);
        reachable[column] = value;
        if (value < 2)
          split.second++;
      }

      p.reset(matrix.size2());
      sort(p, less);
      matroid_apply_column_permutation(matroid, matrix, p);

      split.first--;

      //        std::cout << "split = " << split.first << ", " << split.second << std::endl;

      matroid_permute1(matroid, matrix, split.first, index);

      //        matroid_print (matroid, matrix);

      std::pair <size_t, size_t> witness(split.first, 0);
      while (matrix(witness.first, witness.second) == 0)
      {
        assert (witness.second < nested_minors.width());
        witness.second++;
      }

      separation result(split, witness);
      result.set_special_swap('r', index);
      return result;
    }
  }

  template <typename MatroidType, typename MatrixType>
  separation find_minor_sequence (MatroidType& matroid, MatrixType& matrix, nested_minor_sequence& nested_minors,
      matroid_element_set& extra_elements, logger& log)
  {
    size_t cut = log.size();

    //    std::cout << "Constructing a sequence of nested minors..." << std::flush;

    vector_three_connectivity <MatrixType> column_three_connectivity(matrix, nested_minors.height(), nested_minors.width());
    vector_three_connectivity <matrix_transposed <MatrixType> > row_three_connectivity(view_matrix_transposed(matrix), nested_minors.width(),
        nested_minors.height());

    matroid_transposed <MatroidType> transposed_matroid(matroid);
    matrix_transposed <MatrixType> transposed_matrix(matrix);
    nested_minor_sequence_transposed transposed_nested_minors(nested_minors);

    while (nested_minors.height() < matroid.size1() || nested_minors.width() != matroid.size2())
    {
      if (log.is_updating())
      {
        log.erase(cut);
        log.line() << " + " << nested_minors.size() << " EXT";
        std::cout << log;
      }

      //      std::cout << "\nNew ext ension round... " << nested_minors.height () << " x " << nested_minors.width () << std::endl;

      assert (nested_minors.height() == row_three_connectivity.base());
      assert (nested_minors.height() == column_three_connectivity.dimension());
      assert (nested_minors.width() == row_three_connectivity.dimension());
      assert (nested_minors.width() == column_three_connectivity.base());

      // Simple row extension
      if (find_simple_row_extension(matroid, matrix, nested_minors, row_three_connectivity, column_three_connectivity))
      {
        //        std::cout << "Found a non-parallel non-unit/zero-vector row" << std::endl;
        continue;
      }

      // Simple row extension
      if (find_simple_row_extension(transposed_matroid, transposed_matrix, transposed_nested_minors, column_three_connectivity,
          row_three_connectivity))
      {
        //        std::cout << "Found a non-parallel non-unit/zero-vector column" << std::endl;
        continue;
      }

      //      std::cout << "There are only parallel / unit and zero vectors beyond " << nested_minors.height () << " x " << nested_minors.width ()
      //          << std::endl;
      //      matroid_print (matroid, matrix);

      size_t the_index = 0;
      char type = find_parallel_or_unit_vector(matroid, matrix, nested_minors, row_three_connectivity, column_three_connectivity, the_index);
      if (type == 0)
      {
        //        std::cout << "<< found a 1-separation instead >>" << std::endl;
        return separation(std::make_pair(nested_minors.height(), nested_minors.width()));
      }
      else if (type == 'r')
      {
        //            std::cout << "Working on real matroid " << std::endl;

        separation sep = find_elaborate_extension(matroid, matrix, nested_minors, row_three_connectivity, column_three_connectivity, the_index,
            extra_elements);
        if (sep.is_valid())
        {
          return sep;
        }
      }
      else
      {
        assert (type == 'c');

        //            std::cout << "Working on transposed matroid " << std::endl;

        separation sep = find_elaborate_extension(transposed_matroid, transposed_matrix, transposed_nested_minors, column_three_connectivity,
            row_three_connectivity, the_index, extra_elements);
        if (sep.is_valid())
          return sep.transposed();
      }

      row_three_connectivity.reset(nested_minors.width(), nested_minors.height());
      column_three_connectivity.reset(nested_minors.height(), nested_minors.width());
    }

    if (log.is_updating())
    {
      log.erase(cut);
      log.line() << " + " << nested_minors.size() << " EXT";
      std::cout << log;
    }

    return separation();
  }

}

#endif /* FIND_MINOR_SEQUENCE_HPP_ */
