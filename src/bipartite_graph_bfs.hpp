
//          Copyright Matthias Walter 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BIPARTITE_GRAPH_BFS_HPP_
#define BIPARTITE_GRAPH_BFS_HPP_

#include "../config.h"
#include <vector>
#include <list>

namespace tu {

  struct bipartite_graph_bfs_node
  {
    int distance;
    size_t predecessor;

    inline bool is_reachable () const
    {
      return distance >= 0;
    }
  };

  struct bipartite_graph_dimensions
  {
    typedef size_t index_type;

    bipartite_graph_dimensions (size_t height, size_t width) :
      _height (height), _width (width)
    {

    }

    inline size_t height () const
    {
      return _height;
    }

    inline size_t width () const
    {
      return _width;
    }

    inline size_t size () const
    {
      return _width + _height;
    }

    inline bool is_row (index_type index) const
    {
      return index < _height;
    }

    inline bool is_column (index_type index) const
    {
      return index >= _height;
    }

    inline size_t index_to_row (index_type index) const
    {
      assert (is_row (index));
      return index;
    }

    inline size_t index_to_column (index_type index) const
    {
      assert (is_column (index));
      return index - _height;
    }

    inline index_type row_to_index (size_t row) const
    {
      return row;
    }

    inline index_type column_to_index (size_t column) const
    {
      return column + _height;
    }

    inline std::pair <size_t, size_t> indexes_to_coordinates (index_type first, index_type second) const
    {
      if (first < _height)
        return std::make_pair (first, second - _height);
      else
        return std::make_pair (second, first - _height);
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
  inline bool bipartite_graph_bfs (const MatrixType& matrix, const bipartite_graph_dimensions& dimensions,
      const StartNodesContainerType& start_nodes, const EndNodesContainerType& end_nodes, bool reach_all,
      std::vector <bipartite_graph_bfs_node>& result)
  {
    typedef MatrixType matrix_type;
    typedef std::vector <bipartite_graph_bfs_node> node_vector_type;
    typedef std::list <size_t> node_list_type;

    const size_t size = dimensions.size ();
    const size_t height = dimensions.height ();

    //    std::cout << "Looking at BFS(M) where M is\n";
    //    for (size_t y = 0; y < height; y++)
    //    {
    //        for (size_t x = 0; x < width + 1; x++)
    //        {
    //            std::cout << matrix (y, x) << " ";
    //        }
    //        std::cout << std::endl;
    //    }

    assert (height <= matrix.size1 ());
    assert (dimensions.width () <= matrix.size2 ());

    // Initialize datastructure
    int needed_end_nodes = (reach_all ? end_nodes.size () : 1);
    //    std::cout << "needed end nodes: " << needed_end_nodes << std::endl;
    node_list_type unprocessed;
    result.resize (size);
    for (typename node_vector_type::iterator iter = result.begin (); iter != result.end (); ++iter)
    {
      iter->distance = -1;
      iter->predecessor = 0;
    }
    for (typename StartNodesContainerType::const_iterator iter = start_nodes.begin (); iter != start_nodes.end (); ++iter)
    {
      result[*iter].distance = 0;
      result[*iter].predecessor = *iter;
      unprocessed.push_back (*iter);
    }
    for (typename EndNodesContainerType::const_iterator iter = end_nodes.begin (); iter != end_nodes.end (); ++iter)
    {
      if (result[*iter].distance == 0)
        --needed_end_nodes;
      result[*iter].distance = -2;
    }

    if (needed_end_nodes <= 0)
      return true;

    //    std::cout << "At real start, we still need to find " << needed_end_nodes << " end nodes" << std::endl;

    while (!unprocessed.empty ())
    {
      size_t current_index = unprocessed.front ();
      int current_distance = result[current_index].distance;
      unprocessed.pop_front ();

      //        std::cout << "Looking at (neighbors of) " << current_index << std::endl;

      if (current_index < height)
      {
        // row -> columns
        const size_t row = current_index;
        for (size_t neighbor_index = height; neighbor_index < size; ++neighbor_index)
        {
          const size_t column = neighbor_index - height;
          if (!result[neighbor_index].is_reachable () && matrix (row, column))
          {
            //                    std::cout << "  " << neighbor_index << std::endl;

            if (result[neighbor_index].distance == -2)
              --needed_end_nodes;

            result[neighbor_index].distance = current_distance + 1;
            result[neighbor_index].predecessor = current_index;

            if (needed_end_nodes <= 0)
              return true;

            //                    std::cout << "Pushing back " << neighbor_index << std::endl;
            unprocessed.push_back (neighbor_index);
          }
        }
      }
      else
      {
        // column -> rows
        const size_t column = current_index - height;
        for (size_t neighbor_index = 0; neighbor_index < height; ++neighbor_index)
        {
          const size_t row = neighbor_index;
          if (!result[neighbor_index].is_reachable () && matrix (row, column))
          {
            if (result[neighbor_index].distance == -2)
              --needed_end_nodes;

            result[neighbor_index].distance = current_distance + 1;
            result[neighbor_index].predecessor = current_index;

            if (needed_end_nodes <= 0)
              return true;

            //                    std::cout << "Pushing back " << neighbor_index << std::endl;
            unprocessed.push_back (neighbor_index);
          }
        }
      }
    }
    return false;
  }

}

#endif /* BIPARTITE_GRAPH_BFS_HPP_ */
