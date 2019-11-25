/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef SIGNING_HPP_
#define SIGNING_HPP_

#include <tu/total_unimodularity.hpp>

#include <set>

#include <boost/type_traits/is_const.hpp>

#include <tu/matrix_transposed.hpp>
#include <tu/matrix_permuted.hpp>
#include "matrix_reorder.hpp"
#include <tu/matrix.hpp>
#include "bipartite_graph_bfs.hpp"

namespace unimod
{

  /**
   * Generic function to find a non-zero column and swap it to a given position.
   *
   * @param matrix The given matrix
   * @param column_first First index of a column range
   * @param column_beyond Beyond index of a column range
   * @param row_first First index of a row range
   * @param row_beyond Beyond index of a row range
   * @param target_column Column to be swapped to
   * @return Whether a non-zero column was found.
   */

  template <typename MatrixType>
  bool find_nonzero_column(MatrixType& matrix, size_t column_first, size_t column_beyond, size_t row_first, size_t row_beyond, size_t target_column)
  {
    for (size_t column = column_first; column < column_beyond; ++column)
    {
      for (size_t row = row_first; row < row_beyond; ++row)
      {
        if (matrix(row, column) != 0)
        {
          matrix_permute2(matrix, target_column, column);
          return true;
        }

      }
    }
    return false;
  }

  /**
   * Takes a spanning tree in the bipartite graph of a matrix and a set of nodes to be signed.
   *
   * @param matrix The given matrix
   * @param spanning_tree A spanning tree
   * @param dim Index-mapping for bipartite graph
   * @param nodes Set of nodes
   * @param current_index Current index in the spanning tree
   * @param column The column to be signed
   * @param changes Necessary changes to the entries in the column
   */

  template <typename MatrixType>
  void check_sign(const MatrixType& matrix, const std::vector <bipartite_graph_bfs_node>& spanning_tree, const bipartite_graph_dimensions& dim,
      const std::set <size_t>& nodes, size_t current_index, size_t column, std::map <size_t, bool>& changes)
  {
    /// Root does never change.
    if (spanning_tree[current_index].predecessor == current_index)
    {
      changes[dim.index_to_row(current_index)] = false;
      return;
    }

    /// Search for ancestors until reaching one of the given nodes.
    int value = matrix(dim.index_to_row(current_index), column);
    size_t last, index = current_index;
    do
    {
      last = index;
      index = spanning_tree[index].predecessor;
      std::pair <size_t, size_t> coords = dim.indexes_to_coordinates(index, last);
      value += matrix(coords.first, coords.second);
    }
    while (nodes.find(index) == nodes.end());

    /// If the ancestor is not yet processed, we recurse.
    if (changes.find(dim.index_to_row(index)) == changes.end())
    {
      check_sign(matrix, spanning_tree, dim, nodes, index, column, changes);
    }

    value += matrix(dim.index_to_row(index), column);
    if (changes[dim.index_to_row(index)])
    {
      value += 2;
    }

    value = (value >= 0 ? value : -value) % 4;
    /// If sum (modulo 4) is not 0, we'd like to change the current one
    changes[dim.index_to_row(current_index)] = (value == 2);

    if (value != 0 && value != 2)
    {
      throw std::logic_error("Signing procedure: modulo-sum of cycle was neither 0, nor 2!");
    }
  }

  /**
   * A functor which compares the absolute values.
   */

  template <typename T>
  struct abs_greater
  {
    /**
     * Comparison function
     *
     * @param first First value
     * @param second Second value
     * @return true if and only if the |first| > |second|
     */

    bool operator()(const T& first, const T& second)
    {
      T abs_first = first >= 0 ? first : -first;
      T abs_second = second >= 0 ? second : -second;
      return abs_first > abs_second;
    }
  };

  /**
   * Generic signing function of a matrix, which might also be const.
   * Running time: O (height * width^2)
   *
   * @param matrix The given matrix
   * @param violator Pointer to violator indices to be filled.
   * @return true if and only if the matrix is signed already.
   */

  template <typename M>
  bool sign_matrix(M& matrix, submatrix_indices* violator)
  {
    bool result = true;
    matrix_permuted <M> permuted(matrix);
    size_t handled_rows = 0;

    /// Go trough column by column.
    for (size_t handled_columns = 0; handled_columns < permuted.size2(); ++handled_columns)
    {
      if (find_nonzero_column(permuted, handled_columns, permuted.size2(), 0, handled_rows, handled_columns))
      {
        /// There is a non-zero column right of the already-handled submatrix.

        std::set <size_t> start_nodes;
        std::set <size_t> end_nodes;
        std::set <size_t> all_nodes;

        bipartite_graph_dimensions dim(handled_rows, handled_columns);
        for (size_t row = 0; row < handled_rows; ++row)
        {
          if (permuted(row, handled_columns) != 0)
          {
            size_t index = dim.row_to_index(row);
            if (start_nodes.empty())
              start_nodes.insert(index);
            else
              end_nodes.insert(index);
            all_nodes.insert(index);
          }
        }

        /// Start a BFS on bipartite graph of the submatrix and look for shortest paths from first 1 to all others

        std::vector <bipartite_graph_bfs_node> bfs_result;
        if (!bipartite_graph_bfs(permuted, dim, start_nodes, end_nodes, true, bfs_result))
          throw std::logic_error("Signing procedure: Did not reach all nodes via bfs!");

        /// Evaluate matrix-entries on the shortest paths
        std::map <size_t, bool> changes;
        for (typename std::set <size_t>::const_iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter)
        {
          check_sign(permuted, bfs_result, dim, all_nodes, *iter, handled_columns, changes);
        }

        /// Checking changes
        for (std::map <size_t, bool>::iterator iter = changes.begin(); iter != changes.end(); ++iter)
        {
          if (!iter->second)
            continue;

          if (boost::is_const <M>::value)
          {
            if (violator)
            {
              /// Find the violator, going along the path
              std::set <size_t> violator_rows, violator_columns;

              size_t index = iter->first;
              do
              {
                if (dim.is_row(index))
                  violator_rows.insert(permuted.perm1()(dim.index_to_row(index)));
                else
                  violator_columns.insert(permuted.perm2()(dim.index_to_column(index)));

                index = bfs_result[index].predecessor;
              }
              while (all_nodes.find(index) == all_nodes.end());
              violator_rows.insert(permuted.perm1()(dim.index_to_row(index)));
              violator_columns.insert(permuted.perm2()(handled_columns));

              /// Fill violator data
              violator->rows = submatrix_indices::indirect_array_type(violator_rows.size());
              violator->columns = submatrix_indices::indirect_array_type(violator_columns.size());
              size_t i = 0;
              for (std::set <size_t>::const_iterator iter = violator_rows.begin(); iter != violator_rows.end(); ++iter)
                violator->rows[i++] = *iter;
              i = 0;
              for (std::set <size_t>::const_iterator iter = violator_columns.begin(); iter != violator_columns.end(); ++iter)
                violator->columns[i++] = *iter;
            }
            return false;
          }
          else
          {
            /// We are not just testing, so swap the sign on a one.
            size_t real_row = permuted.perm1()(dim.index_to_row(iter->first));
            size_t real_column = permuted.perm2()(handled_columns);
            matrix_set_value(matrix, real_row, real_column, -matrix(real_row, real_column));

            result = false;
          }
        }

        matrix_reorder_rows(permuted, handled_rows, permuted.size1(), handled_columns, permuted.size2(), abs_greater <int> ());

        /// Augment submatrix by rows with 1 in the new column.
        while (handled_rows < permuted.size1())
        {
          if (permuted(handled_rows, handled_columns) == 0)
            break;
          else
            ++handled_rows;
        }
      }
      else
      {
        /// Handled upper-left submatrix and lower-right submatrix are disconnected
        for (size_t column = handled_columns; column < permuted.size2(); ++column)
        {
          size_t count = 0;
          for (size_t row = handled_rows; row < permuted.size1(); ++row)
          {
            if (permuted(row, column) != 0)
              ++count;
          }

          /// A zero column can be skipped, as it is handled by definition.
          if (count > 0)
          {
            /// Found a nonzero column and swap ones to the top.
            matrix_reorder_rows(permuted, handled_rows, permuted.size1(), handled_columns, permuted.size2(), abs_greater <int> ());
            while (handled_rows < permuted.size1())
            {
              if (permuted(handled_rows, handled_columns) == 0)
                break;
              else
                ++handled_rows;
            }

            break;
          }
        }
      }
    }

    return result;
  }

}

#endif /* SIGNING_HPP_ */
