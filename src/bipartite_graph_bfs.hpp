/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef BIPARTITE_GRAPH_BFS_HPP_
#define BIPARTITE_GRAPH_BFS_HPP_

#include <vector>
#include <list>

namespace unimod
{

  /**
   * A node for a breadth first search on
   * bipartite graphs.
   */

  struct bipartite_graph_bfs_node
  {
    /// Combinatorial distance to the search root or -1 if not reachable
    int distance;

    /// Index of the predecessor in the search tree
    size_t predecessor;

    /**
     * @return true if and only if this node is reachable
     */

    inline bool is_reachable() const
    {
      return distance >= 0;
    }
  };

  /**
   * A helper struct for breadth first search on
   * bipartite graphs to map row- and column-indices
   * into on index-range.
   */

  struct bipartite_graph_dimensions
  {
    typedef size_t index_type;

    /**
     * Constructs the helper object with the given
     * number of rows and columns.
     *
     * @param height Number of rows
     * @param width Number of columns
     */

    bipartite_graph_dimensions(size_t height, size_t width) :
      _height(height), _width(width)
    {

    }

    /**
     * @return Number of rows
     */

    inline size_t height() const
    {
      return _height;
    }

    /**
     * @return Number of columns
     */

    inline size_t width() const
    {
      return _width;
    }

    /**
     * @return Number of all node-indices
     */

    inline size_t size() const
    {
      return _width + _height;
    }

    /**
     * @param index A given node-index
     * @return true if it indexes a row
     */

    inline bool is_row(index_type index) const
    {
      return index < _height;
    }

    /**
     * @param index A given node-index
     * @return true if it indexes a column
     */

    inline bool is_column(index_type index) const
    {
      return index >= _height;
    }

    /**
     * @param index A given node-index which indexes a row
     * @return Corresponding row index
     */

    inline size_t index_to_row(index_type index) const
    {
      assert(is_row(index));
      return index;
    }

    /**
     * @param index A given node-index which indexes a column
     * @return Corresponding column index
     */

    inline size_t index_to_column(index_type index) const
    {
      assert(is_column(index));
      return index - _height;
    }

    /**
     * @param row A given row-index
     * @return Correspoding node-index
     */

    inline index_type row_to_index(size_t row) const
    {
      return row;
    }

    /**
     * @param column A given column-index
     * @return Corresponding node-index
     */

    inline index_type column_to_index(size_t column) const
    {
      return column + _height;
    }

    /**
     * Maps two node-indices at once and puts the resulting
     * row and column and indices in this fixed ording into
     * the resulting pair.
     *
     * @param first A first node-index
     * @param second A second node-index
     * @return A pair consisting of the corresponding column-index and row-index.
     */

    inline std::pair <size_t, size_t> indexes_to_coordinates(index_type first, index_type second) const
    {
      if (first < _height)
        return std::make_pair(first, second - _height);
      else
        return std::make_pair(second, first - _height);
    }

  private:
    size_t _height;
    size_t _width;
  };

  /**
   * Does a breadth-first-search on the bipartite graph defined by a submatrix of a given matrix.
   * Running time: O(height * width)
   *
   * @param matrix Given matrix
   * @param dimensions Object for row/column-to-index calculations
   * @param start_nodes Set of starting nodes for search
   * @param end_nodes Set of nodes to be found
   * @param reach_all Whether we'd like to find all end nodes, or just one
   * @param result Whether all resp. one end node was found.
   * @return true iff enough of the target nodes were found.
   */

  template <typename MatrixType, typename StartNodesContainerType, typename EndNodesContainerType>
  inline bool bipartite_graph_bfs(const MatrixType& matrix, const bipartite_graph_dimensions& dimensions, const StartNodesContainerType& start_nodes,
      const EndNodesContainerType& end_nodes, bool reach_all, std::vector <bipartite_graph_bfs_node>& result)
  {
    typedef std::vector <bipartite_graph_bfs_node> node_vector_type;
    typedef std::list <size_t> node_list_type;

    const size_t size = dimensions.size();
    const size_t height = dimensions.height();

    assert(height <= matrix.size1());
    assert(dimensions.width() <= matrix.size2());

    /// Initialize datastructure
    int needed_end_nodes = (reach_all ? end_nodes.size() : 1);
    node_list_type unprocessed;
    result.resize(size);
    for (typename node_vector_type::iterator iter = result.begin(); iter != result.end(); ++iter)
    {
      iter->distance = -1;
      iter->predecessor = 0;
    }
    for (typename StartNodesContainerType::const_iterator iter = start_nodes.begin(); iter != start_nodes.end(); ++iter)
    {
      result[*iter].distance = 0;
      result[*iter].predecessor = *iter;
      unprocessed.push_back(*iter);
    }
    for (typename EndNodesContainerType::const_iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter)
    {
      if (result[*iter].distance == 0)
        --needed_end_nodes;
      result[*iter].distance = -2;
    }

    if (needed_end_nodes <= 0)
      return true;

    /// Start bfs
    while (!unprocessed.empty())
    {
      size_t current_index = unprocessed.front();
      int current_distance = result[current_index].distance;
      unprocessed.pop_front();

      if (dimensions.is_row(current_index))
      {
        /// row -> columns
        const size_t row = current_index;
        for (size_t neighbor_index = height; neighbor_index < size; ++neighbor_index)
        {
          const size_t column = neighbor_index - height;
          if (!result[neighbor_index].is_reachable() && matrix(row, column))
          {
            if (result[neighbor_index].distance == -2)
              --needed_end_nodes;

            result[neighbor_index].distance = current_distance + 1;
            result[neighbor_index].predecessor = current_index;

            if (needed_end_nodes <= 0)
              return true;

            unprocessed.push_back(neighbor_index);
          }
        }
      }
      else
      {
        assert(dimensions.is_column(current_index));

        /// column -> rows
        const size_t column = current_index - height;
        for (size_t neighbor_index = 0; neighbor_index < height; ++neighbor_index)
        {
          const size_t row = neighbor_index;
          if (!result[neighbor_index].is_reachable() && matrix(row, column))
          {
            if (result[neighbor_index].distance == -2)
              --needed_end_nodes;

            result[neighbor_index].distance = current_distance + 1;
            result[neighbor_index].predecessor = current_index;

            if (needed_end_nodes <= 0)
              return true;

            unprocessed.push_back(neighbor_index);
          }
        }
      }
    }
    return false;
  }

}

#endif /* BIPARTITE_GRAPH_BFS_HPP_ */
