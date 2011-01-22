/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef GEN_NETWORK_HPP_
#define GEN_NETWORK_HPP_

#include "gen_generic.hpp"
#include "matrix.hpp"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/two_bit_color_map.hpp>

template <typename OutputIterator, typename Graph, typename IndexMap, typename Vertex>
OutputIterator find_path(Graph graph, IndexMap index_map, Vertex s, Vertex t, OutputIterator result)
{
  if (s == t)
  {
    *result++ = s;
    return result;
  }

  typedef boost::graph_traits <Graph> graph_traits_t;
  typedef typename graph_traits_t::vertex_descriptor vertex_t;

  /// Define a color map for DFS
  typedef boost::two_bit_color_map <> color_map_t;
  color_map_t color_map(boost::num_vertices(graph));

  /// Declare predecessor map
  typedef std::vector <vertex_t> predecessors_t;
  typedef boost::iterator_property_map <typename predecessors_t::iterator, IndexMap> predecessor_map_t;
  typedef boost::predecessor_recorder <predecessor_map_t, IndexMap> predecessor_recorder_t;

  predecessors_t predecessors(boost::num_vertices(graph), graph_traits_t::null_vertex());
  predecessor_map_t predecessor_map(predecessors.begin(), index_map);

  boost::depth_first_visit(graph, s, boost::make_dfs_visitor(record_predecessors(predecessor_map, boost::on_tree_edge())), color_map);

  vertex_t current_vertex = t;
  while (current_vertex != s)
  {
    vertex_t next_vertex = boost::get(predecessor_map, current_vertex);
    if (next_vertex == graph_traits_t::null_vertex())
      return result;

    *result++ = current_vertex;

    current_vertex = next_vertex;
  }

  *result++ = s;
  return result;
}

class network_matrix_generator: public matrix_generator
{
public:
  network_matrix_generator(size_t height, size_t width, tu::log_level level) :
    matrix_generator("network", height, width, level)
  {

  }

  virtual ~network_matrix_generator()
  {

  }

  template <typename MatrixType>
  void generate(MatrixType& matrix)
  {
    size_t nodes = matrix.size1() + 1;

    if (_level != tu::LOG_QUIET)
      std::cerr << "Creating a spanning tree with " << nodes << " nodes..." << std::flush;

    /// Create a spanning tree
    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> tree_graph_t;
    typedef boost::graph_traits <tree_graph_t> tree_traits_t;

    tree_graph_t tree_graph(nodes);
    std::vector <tree_traits_t::vertex_descriptor> used_vertices;

    tree_traits_t::vertex_iterator vertex_iter, vertex_beyond;
    for (boost::tie(vertex_iter, vertex_beyond) = boost::vertices(tree_graph); vertex_iter != vertex_beyond; ++vertex_iter)
    {
      if (!used_vertices.empty())
      {
        boost::uniform_int <int> dist(0, used_vertices.size() - 1);
        tree_traits_t::vertex_descriptor other = used_vertices[dist(_rng)];

        boost::add_edge(*vertex_iter, other, tree_graph);
      }
      used_vertices.push_back(*vertex_iter);
    }

    if (_level != tu::LOG_QUIET)
      std::cerr << " done.\nAdding edges and filling matrix..." << std::flush;

    for (size_t column = 0; column < _matrix.size2(); ++column)
    {
      /// Choose an edge not in the tree
      boost::uniform_int <int> dist(0, nodes - 1);
      tree_traits_t::vertex_descriptor u, v;
      do
      {
        u = boost::vertex(dist(_rng), tree_graph);
        v = boost::vertex(dist(_rng), tree_graph);
      }
      while (u == v || boost::edge(u, v, tree_graph).second);

      std::set <tree_traits_t::vertex_descriptor> path_vertices;
      find_path(tree_graph, boost::get(boost::vertex_index, tree_graph), u, v, std::inserter(path_vertices, path_vertices.end()));

      size_t row = 0;
      tree_traits_t::edge_iterator edge_iter, edge_beyond;
      for (boost::tie(edge_iter, edge_beyond) = boost::edges(tree_graph); edge_iter != edge_beyond; ++edge_iter)
      {
        u = boost::source(*edge_iter, tree_graph);
        v = boost::target(*edge_iter, tree_graph);

        _matrix(row, column) = (path_vertices.find(u) != path_vertices.end() && path_vertices.find(v) != path_vertices.end()) ? 1 : 0;

        ++row;
      }
    }

    if (_level != tu::LOG_QUIET)
      std::cerr << " done.\nCorrecting the signs..." << std::flush;
    sign();
    if (_level != tu::LOG_QUIET)
      std::cerr << " done." << std::endl;
  }

  virtual void generate()
  {
    if (_height <= _width)
    {
      generate(_matrix);
    }
    else
    {
      tu::matrix_transposed <tu::integer_matrix> transposed(_matrix);
      generate(transposed);
    }
  }

};
#endif
