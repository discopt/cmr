/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef FIND_WHEEL_MINOR_HPP_
#define FIND_WHEEL_MINOR_HPP_

#include "matrix.hpp"
#include "matrix_permuted.hpp"
#include "matrix_reorder.hpp"
#include "matroid.hpp"
#include "matroid_permuted.hpp"
#include "matroid_reorder.hpp"
#include "separation.hpp"
#include "matrix_modified.hpp"
#include "bipartite_graph_bfs.hpp"
#include "comparators.hpp"

namespace tu {

  /**
   * A matrix proxy which makes an upper-left part zero.
   */

  struct zero_block_matrix_modifier
  {
    typedef int value_type;

    /**
     * Constructs the modifier
     *
     * @param height Number of rows in the zero-part
     * @param width Number of columns in the zero-part
     */

    zero_block_matrix_modifier (size_t height, size_t width) :
      height_(height), width_(width)
    {

    }

    /**
     * Applies the modification
     *
     * @param i row
     * @param j column
     * @param value Original value
     * @return Modified value
     */

    int operator () (size_t i, size_t j, int value)
    {
      if (i < height_ && j < width_)
        return 0;
      else
        return value;
    }

  private:
    size_t height_;
    size_t width_;
  };

  /**
   * Searches for a W3 minor of a given matroid. It might find a 1- or 2-separation instead.
   * The matroid bust be at least 3 x 3.
   *
   * @param permuted_matroid Given matroid
   * @param permuted_matrix Representation matrix for the given matroid
   * @param extra_elements Set of matroid-elements to be filled with pivot elements
   * @return Separation, if one was found
   */

  template <typename MatroidType, typename MatrixType>
  separation find_wheel_minor (matroid_permuted <MatroidType>& permuted_matroid, matrix_permuted <MatrixType>& permuted_matrix,
      matroid_element_set& extra_elements)
  {
    assert (permuted_matrix.size1() >= 3 && permuted_matroid.size2() >= 3);

    matroid_reorder_columns(permuted_matroid, permuted_matrix, 0, 1, 0, permuted_matroid.size2(), std::greater <int>());
    size_t count_first_row_ones = matrix_count_property_column_series(permuted_matrix, 0, 1, 0, permuted_matrix.size2(), is_non_zero());

    /// Trivial 1-separation
    if (count_first_row_ones == 0)
    {
      return separation(std::make_pair(1, 0));
    }

    matroid_reorder_rows(permuted_matroid, permuted_matrix, 1, permuted_matroid.size1(), 0, 1, std::greater <int>());
    size_t count_first_column_ones = matrix_count_property_row_series(permuted_matrix, 0, permuted_matroid.size1(), 0, 1, is_non_zero());

    /// 1- or 2-separation
    if (count_first_row_ones == 1)
    {
      if (count_first_column_ones == 0)
      {
        return separation(std::make_pair(1, 1));
      }
      else
      {
        return separation(std::make_pair(1, 1), std::make_pair(1, 0));
      }
    }
    else if (count_first_column_ones == 1)
    {
      return separation(std::make_pair(1, 1), std::make_pair(0, 1));
    }

    assert ((permuted_matrix(0,0) == 1) && (permuted_matrix(1,0) == 1) && (permuted_matrix(0,1) == 1));

    /// Ensure we have a 2x2 block of ones
    if (permuted_matrix(1, 1) != 1)
    {
      matroid_binary_pivot(permuted_matroid, permuted_matrix, 0, 0);
      extra_elements.insert(permuted_matroid.name1(0));
      extra_elements.insert(permuted_matroid.name2(0));
    }

    assert ((permuted_matrix(0,0) == 1) && (permuted_matrix(1,0) == 1) && (permuted_matrix(0,1) == 1) && (permuted_matrix(1,1) == 1));

    /// Grow the block to a set-maximal one.
    matroid_reorder_columns(permuted_matroid, permuted_matrix, 0, 2, 2, permuted_matroid.size2(), std::greater <int>());
    size_t block_width = 2 + matrix_count_property_column_series(permuted_matrix, 0, 2, 2, permuted_matrix.size2(), is_all_ones());

    matroid_reorder_rows(permuted_matroid, permuted_matrix, 2, permuted_matroid.size1(), 0, block_width, std::greater <int>());

    size_t block_height = 2 + matrix_count_property_row_series(permuted_matrix, 2, permuted_matrix.size1(), 0, block_width, is_all_ones());

    /// Search for a path in BFS graph
    zero_block_matrix_modifier modifier(block_height, block_width);
    matrix_modified <matrix_permuted <MatrixType> , zero_block_matrix_modifier> modified_matrix(permuted_matrix, modifier);
    bipartite_graph_dimensions dim(permuted_matrix.size1(), permuted_matrix.size2());
    std::vector <bipartite_graph_dimensions::index_type> start_nodes(block_height);
    for (size_t i = 0; i < block_height; i++)
      start_nodes[i] = dim.row_to_index(i);
    std::vector <bipartite_graph_dimensions::index_type> end_nodes(block_width);
    for (size_t i = 0; i < block_width; i++)
      end_nodes[i] = dim.column_to_index(i);

    std::vector <bipartite_graph_bfs_node> bfs_result;
    bool found_path = bipartite_graph_bfs(modified_matrix, dim, start_nodes, end_nodes, false, bfs_result);

    /// There must be a 2-separation.
    if (!found_path)
    {
      std::pair <size_t, size_t> split(0, 0);

      /// Swap unreachable rows to top
      std::vector <int> reachable(permuted_matrix.size1());
      for (size_t i = 0; i < permuted_matrix.size1(); ++i)
      {
        const bipartite_graph_bfs_node& node = bfs_result[dim.row_to_index(i)];
        int value = (node.is_reachable() ? (node.distance > 0 ? 2 : 1) : 0);
        reachable[permuted_matrix.perm1()(i)] = value;
        if (value == 0)
          split.first++;
      }
      vector_less <int, std::less <int> > less(reachable, std::less <int>());

      sort(permuted_matrix.perm1(), less);
      permuted_matroid.perm1() = permuted_matrix.perm1();

      /// Swap unreachable columns to left
      reachable.resize(permuted_matrix.size2());
      for (size_t i = 0; i < permuted_matrix.size2(); ++i)
      {
        const bipartite_graph_bfs_node& node = bfs_result[dim.column_to_index(i)];
        int value = (node.is_reachable() ? 2 : (node.distance == -2 ? 1 : 0));
        reachable[permuted_matrix.perm2()(i)] = value;
        if (value < 2)
          split.second++;
      }

      sort(permuted_matrix.perm2(), less);
      permuted_matroid.perm2() = permuted_matrix.perm2();

      return separation(split, std::make_pair(split.first, split.second - 1));
    }

    /// Go back along the path to find the nearest endpoint.
    bipartite_graph_dimensions::index_type nearest_end = 0;
    for (std::vector <bipartite_graph_dimensions::index_type>::const_iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter)
    {
      if (bfs_result[*iter].is_reachable())
        nearest_end = *iter;
    }

    assert (bfs_result[nearest_end].is_reachable());

    size_t w3_one_column = dim.index_to_column(nearest_end);
    size_t nearest_distance = bfs_result[nearest_end].distance + 1;

    assert (nearest_distance % 2 == 0);

    size_t last_index = nearest_end;
    size_t current_index = bfs_result[last_index].predecessor;

    /// Find indices for basic W3-parts.
    size_t w3_one_row = 0;
    size_t w3_path_column = 0;
    size_t w3_path_row = dim.index_to_row(current_index);
    size_t w3_zero_column = 0;
    while (permuted_matrix(w3_path_row, w3_zero_column) != 0)
    {
      w3_zero_column++;
      assert (w3_zero_column < block_width);
    }

    /// Apply path shortening technique.
    while (last_index != current_index)
    {
      std::pair <size_t, size_t> coords = dim.indexes_to_coordinates(current_index, last_index);

      if ((bfs_result[current_index].distance % 2 == 0) && (bfs_result[current_index].distance >= 2) && (bfs_result[current_index].distance + 2
          < (int) nearest_distance))
      {
        matroid_binary_pivot(permuted_matroid, permuted_matrix, coords.first, coords.second);
        extra_elements.insert(permuted_matroid.name1(coords.first));
        extra_elements.insert(permuted_matroid.name2(coords.second));
      }

      if (bfs_result[current_index].distance == 1)
      {
        assert (dim.is_column(current_index));

        w3_path_column = dim.index_to_column(current_index);
      }
      else if (bfs_result[current_index].distance == 0)
      {
        assert (dim.is_row(current_index));

        w3_one_row = dim.index_to_row(current_index);
      }

      last_index = current_index;
      current_index = bfs_result[current_index].predecessor;
    }

    /// Collect more W3-indices.
    size_t w3_zero_row = 0;
    while (permuted_matrix(w3_zero_row, w3_path_column) != 0)
    {
      w3_zero_row++;
      assert (w3_zero_row< block_height);
    }

    assert (permuted_matrix (w3_one_row, w3_one_column) == 1);
    assert (permuted_matrix (w3_one_row, w3_zero_column) == 1);
    assert (permuted_matrix (w3_one_row, w3_path_column) == 1);
    assert (permuted_matrix (w3_zero_row, w3_one_column) == 1);
    assert (permuted_matrix (w3_zero_row, w3_zero_column) == 1);
    assert (permuted_matrix (w3_zero_row, w3_path_column) == 0);
    assert (permuted_matrix (w3_path_row, w3_one_column) == 1);
    assert (permuted_matrix (w3_path_row, w3_zero_column) == 0);
    assert (permuted_matrix (w3_path_row, w3_path_column) == 1);

    /// Some normalization
    if (w3_zero_row > w3_one_row)
    {
      matroid_permute1(permuted_matroid, permuted_matrix, w3_one_row, w3_zero_row);
      std::swap(w3_one_row, w3_zero_row);
    }
    if (w3_one_row > w3_path_row)
    {
      matroid_permute1(permuted_matroid, permuted_matrix, w3_path_row, w3_one_row);
      std::swap(w3_path_row, w3_one_row);
    }

    if (w3_zero_column > w3_one_column)
    {
      matroid_permute2(permuted_matroid, permuted_matrix, w3_one_column, w3_zero_column);
      std::swap(w3_one_column, w3_zero_column);
    }
    if (w3_one_column > w3_path_column)
    {
      matroid_permute2(permuted_matroid, permuted_matrix, w3_path_column, w3_one_column);
      std::swap(w3_path_column, w3_one_column);
    }

    matroid_permute1(permuted_matroid, permuted_matrix, 0, w3_zero_row);
    matroid_permute1(permuted_matroid, permuted_matrix, 1, w3_one_row);
    matroid_permute1(permuted_matroid, permuted_matrix, 2, w3_path_row);

    matroid_permute2(permuted_matroid, permuted_matrix, 0, w3_zero_column);
    matroid_permute2(permuted_matroid, permuted_matrix, 1, w3_one_column);
    matroid_permute2(permuted_matroid, permuted_matrix, 2, w3_path_column);

    return separation();
  }

}

#endif /* FIND_WHEEL_MINOR_HPP_ */
