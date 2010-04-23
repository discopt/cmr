/*
 * signing.cpp
 *
 *  Created on: Dec 9, 2009
 *      Author: xammy
 */

#include "total_unimodularity.hpp"

#include "../config.h"
#include "matrix_transposed.hpp"
#include "matrix_permuted.hpp"
#include "matrix_reorder.hpp"
#include "matrix.hpp"
#include "bipartite_graph_bfs.hpp"

#include <set>
#include <boost/type_traits/is_const.hpp>

namespace tu {

  /**
   * Generic function to find a non-zero column and swap it to a given position.
   * 
   * @param matrix Given matrix
   * @param column_first begin index of a column range.
   * @param column_beyond beyond index of a column range
   * @param row_first begin index of a row range
   * @param row_beyond beyond index of a row range
   * @param target_column column to be swapped to
   * @return Whether a non-zero column was found.
   */

  template <typename MatrixType>
  bool find_nonzero_column (MatrixType& matrix, size_t column_first, size_t column_beyond, size_t row_first, size_t row_beyond, size_t target_column)
  {
    for (size_t column = column_first; column < column_beyond; ++column)
    {
      for (size_t row = row_first; row < row_beyond; ++row)
      {
        if (matrix (row, column) != 0)
        {
          matrix_permute2 (matrix, target_column, column);
          return true;
        }

      }
    }
    return false;
  }

  /**
   * 
   * 
   * @param matrix
   * @param spanning_tree
   * @param dim
   * @param nodes
   * @param current_index
   * @param column
   * @param changes
   */

  template <typename MatrixType>
  void check_sign (const MatrixType& matrix, const std::vector <bipartite_graph_bfs_node>& spanning_tree, const bipartite_graph_dimensions& dim,
      const std::set <size_t>& nodes, size_t current_index, size_t column, std::map <size_t, bool>& changes)
  {
    // Root does not change.
    if (spanning_tree[current_index].predecessor == current_index)
    {
      //        std::cout << "root at " << dim.index_to_row (current_index) << " does not change" << std::endl;
      changes[dim.index_to_row (current_index)] = false;
      return;
    }

    //    std::cout << "current_index = " << current_index << std::endl;

    // Search for ancestors until reaching one of the given nodes.
    int value = matrix (dim.index_to_row (current_index), column);
    //    std::cout << "value at current_index = " << value << std::endl;
    size_t last, index = current_index;
    do
    {
      last = index;
      index = spanning_tree[index].predecessor;
      std::pair <size_t, size_t> coords = dim.indexes_to_coordinates (index, last);
      value += matrix (coords.first, coords.second);

      //        std::cout << "adding value from " << coords.first << ", " << coords.second << std::endl;
    }
    while (nodes.find (index) == nodes.end ());

    // Ancestor not yet processed
    if (changes.find (dim.index_to_row (index)) == changes.end ())
    {
      check_sign (matrix, spanning_tree, dim, nodes, index, column, changes);
    }

    //    std::cout << "adding value at " << dim.index_to_row (index) << std::endl;

    value += matrix (dim.index_to_row (index), column);
    if (changes[dim.index_to_row (index)])
    {
      //        std::cout << "adding 2" << std::endl;
      value += 2;
    }

    //    std::cout << "value for row " << dim.index_to_row (current_index) << " is " << value << std::endl;

    value = (value >= 0 ? value : -value) % 4;
    // If sum (modulo 4) is not 0, we'd like to change the current one
    changes[dim.index_to_row (current_index)] = (value == 2);

    if (value != 0 && value != 2)
    {
      throw std::logic_error ("Signing procedure: modulo-sum of cycle was neither 0, nor 2!");
    }
  }

  template <typename T>
  struct abs_greater
  {
    bool operator() (const T& first, const T& second)
    {
      T abs_first = first >= 0 ? first : -first;
      T abs_second = second >= 0 ? second : -second;
      return abs_first > abs_second;
    }
  };

  /// Generic signing function. M might be const.
  /// Running time: O(height * width^2)

  template <typename M>
  bool sign_matrix (M& matrix, submatrix_indices* violator)
  {
    matrix_permuted <M> permuted (matrix);
    size_t handled_rows = 0;
    for (size_t handled_columns = 0; handled_columns < permuted.size2 (); ++handled_columns)
    {
      //        std::cout << "handled upper-left = " << handled_rows << " x " << handled_columns << std::endl;
      //        matrix_print (permuted);
      //
      //        std::cout << "Original is:\n";
      //        matrix_print (permuted.data ());

      if (find_nonzero_column (permuted, handled_columns, permuted.size2 (), 0, handled_rows, handled_columns))
      {
        // There is a non-zero column right of the already-handled submatrix. 

        std::set <size_t> start_nodes;
        std::set <size_t> end_nodes;
        std::set <size_t> all_nodes;

        bipartite_graph_dimensions dim (handled_rows, handled_columns);
        for (size_t row = 0; row < handled_rows; ++row)
        {
          if (permuted (row, handled_columns) != 0)
          {
            size_t index = dim.row_to_index (row);
            if (start_nodes.empty ())
              start_nodes.insert (index);
            else
              end_nodes.insert (index);
            all_nodes.insert (index);
          }
        }

        // Start a BFS on BFS(submatrix), look for shortest paths from first one to all others

        std::vector <bipartite_graph_bfs_node> bfs_result;
        if (!bipartite_graph_bfs (permuted, dim, start_nodes, end_nodes, true, bfs_result))
          throw std::logic_error ("Signing procedure: Did not reach all nodes via bfs!");

        // Evaluate matrix-entries on the shortest paths
        std::map <size_t, bool> changes;
        for (typename std::set <size_t>::const_iterator iter = end_nodes.begin (); iter != end_nodes.end (); ++iter)
        {
          check_sign (permuted, bfs_result, dim, all_nodes, *iter, handled_columns, changes);
        }

        //            std::cout << "Real columns are";
        //            for (size_t i = 0; i < handled_columns; i++)
        //            {
        //                std::cout << " " << permuted.perm2 () (i);
        //            }
        //            std::cout << std::endl;

        // Checking changes
        for (std::map <size_t, bool>::iterator iter = changes.begin (); iter != changes.end (); ++iter)
        {
          if (!iter->second)
            continue;

          if (boost::is_const <M>::value)
          {
            //                    std::cout << "Found a violation at (permuted) current column " << handled_columns << std::endl;

            if (violator)
            {
              // Find the violator, going along the path
              std::set <size_t> violator_rows, violator_columns;

              size_t index = iter->first;
              do
              {
                if (dim.is_row (index))
                  violator_rows.insert (permuted.perm1 () (dim.index_to_row (index)));
                else
                  violator_columns.insert (permuted.perm2 () (dim.index_to_column (index)));

                index = bfs_result[index].predecessor;
              }
              while (all_nodes.find (index) == all_nodes.end ());
              violator_rows.insert (permuted.perm1 () (dim.index_to_row (index)));
              violator_columns.insert (permuted.perm2 () (handled_columns));

              // Fill violator data
              violator->rows = submatrix_indices::indirect_array_type (violator_rows.size ());
              violator->columns = submatrix_indices::indirect_array_type (violator_columns.size ());
              size_t i = 0;
              for (std::set <size_t>::const_iterator iter = violator_rows.begin (); iter != violator_rows.end (); ++iter)
                violator->rows[i++] = *iter;
              i = 0;
              for (std::set <size_t>::const_iterator iter = violator_columns.begin (); iter != violator_columns.end (); ++iter)
                violator->columns[i++] = *iter;
            }
            return false;
          }
          else
          {
            //                    std::cout << "Really changing at " << (dim.index_to_row (iter->first)) << "," << handled_columns
            //                            << std::endl;
            //                    matrix_print (permuted);

            // We are not just testing, so swap the sign on a one.

            size_t real_row = permuted.perm1 () (dim.index_to_row (iter->first));
            size_t real_column = permuted.perm2 () (handled_columns);
            matrix_set_value (matrix, real_row, real_column, -matrix (real_row, real_column));

            //                    std::cout << "Original after change is\n";
            //                    matrix_print (permuted.data ());
          }
        }

        matrix_reorder_rows (permuted, handled_rows, permuted.size1 (), handled_columns, permuted.size2 (), abs_greater <int> ());

        //            std::cout << "handled rows is increased from " << handled_rows;

        while (handled_rows < permuted.size1 ())
        {
          if (permuted (handled_rows, handled_columns) == 0)
            break;
          else
            ++handled_rows;
        }

        //            std::cout << " to " << handled_rows << std::endl;
      }
      else
      {
        //            std::cout << "Disconnected: increasing handled_rows from " << handled_rows;

        // Handled upper-left submatrix and lower-right submatrix are disconnected 

        for (size_t column = handled_columns; column < permuted.size2 (); ++column)
        {
          size_t count = 0;
          for (size_t row = handled_rows; row < permuted.size1 (); ++row)
          {
            if (permuted (row, column) != 0)
              ++count;
          }

          // A zero column can be skipped, as it is handled by definition.

          if (count > 0)
          {
            // Found a nonzero column and swap ones to the top.

            //                    std::cout << "reordering for disconnected case starting at " << handled_rows << ","
            //                            << handled_columns << std::endl;

            matrix_reorder_rows (permuted, handled_rows, permuted.size1 (), handled_columns, permuted.size2 (), abs_greater <int> ());

            while (handled_rows < permuted.size1 ())
            {
              if (permuted (handled_rows, handled_columns) == 0)
                break;
              else
                ++handled_rows;
            }

            break;
          }
        }

        //            std::cout << " to " << handled_rows << std::endl;
      }
    }

    return true;
  }

  /// Returns true, iff the given matrix is already a signed version of its support matrix.
  /// Running time: O(height * width * min (height, width) )

  bool is_signed_matrix (const boost::numeric::ublas::matrix <int>& matrix)
  {
    if (matrix.size2 () > matrix.size1 ())
    {
      matrix_transposed <const boost::numeric::ublas::matrix <int> > transposed (matrix);
      return sign_matrix (transposed, NULL);
    }
    else
    {
      return sign_matrix (matrix, NULL);
    }
  }

  /// Returns true, iff the given matrix is already a signed version of its support matrix.
  /// If not, violator describes a violating submatrix.
  /// Running time: O(height * width * min (height, width) )

  bool is_signed_matrix (const boost::numeric::ublas::matrix <int>& matrix, submatrix_indices& violator)
  {
    if (matrix.size2 () > matrix.size1 ())
    {
      matrix_transposed <const boost::numeric::ublas::matrix <int> > transposed (matrix);
      return sign_matrix (transposed, &violator);
    }
    else
    {
      return sign_matrix (matrix, &violator);
    }
  }

  /// Returns true, iff the given matrix is already a signed version of its support matrix.
  /// If not, the given matrix is changed to be such a signed version.
  /// Running time: O(height * width * min (height, width) )

  bool sign_matrix (boost::numeric::ublas::matrix <int>& matrix)
  {
    if (matrix.size2 () > matrix.size1 ())
    {
      matrix_transposed <boost::numeric::ublas::matrix <int> > transposed (matrix);
      return sign_matrix (transposed, NULL);
    }
    else
    {
      return sign_matrix (matrix, NULL);
    }
  }

  /// Drops all signs in the given matrix.

  void support_matrix (boost::numeric::ublas::matrix <int>& matrix)
  {
    for (size_t i = 0; i < matrix.size1 (); ++i)
    {
      for (size_t j = 0; j < matrix.size2 (); ++j)
      {
        const int value = matrix (i, j);
        if (value != 0)
        {
          matrix (i, j) = 1;
        }
      }
    }
  }

}
