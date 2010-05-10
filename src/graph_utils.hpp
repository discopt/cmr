/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef GRAPH_UTILS_HPP_
#define GRAPH_UTILS_HPP_

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace tu {

  namespace util {

    /**
     * Reconnects an edge in a graph, retaining the edge property, by replacing the target vertex.
     *
     * @param graph A given graph
     * @param edge An edge in the graph
     * @param old_vertex One of the vertices of the edge
     * @param new_vertex Another vertex in the graph, the edge should have as target
     */

    template <typename Graph>
    void reconnect_edge_target (Graph& graph, typename boost::graph_traits <Graph>::edge_descriptor edge,
        typename boost::graph_traits <Graph>::vertex_descriptor old_vertex, typename boost::graph_traits <Graph>::vertex_descriptor new_vertex)
    {
      if (old_vertex != new_vertex)
      {
        boost::add_edge(boost::source(edge, graph), new_vertex, *((typename Graph::edge_property_type*) edge.get_property()), graph);
        boost::remove_edge(edge, graph);
      }
    }

    /**
     * Reconnects an edge in a graph, retaining the edge property, by replacing the target vertex.
     *
     * @param graph A given graph
     * @param edge An edge in the graph
     * @param old_vertex One of the vertices of the edge
     * @param new_vertex Another vertex in the graph, the edge should have as source
     */

    template <typename Graph>
    void reconnect_edge_source (Graph& graph, typename boost::graph_traits <Graph>::edge_descriptor edge,
        typename boost::graph_traits <Graph>::vertex_descriptor old_vertex, typename boost::graph_traits <Graph>::vertex_descriptor new_vertex)
    {
      if (old_vertex != new_vertex)
      {
        boost::add_edge(boost::target(edge, graph), new_vertex, *((typename Graph::edge_property_type*) edge.get_property()), graph);
        boost::remove_edge(edge, graph);
      }
    }

    /**
     * Reconnects an edge in a graph, retaining the edge property, by replacing a given old vertex by a new one.
     *
     * @param graph The given graph
     * @param edge An edge in the graph
     * @param old_vertex One of the vertices of the edge
     * @param new_vertex Another vertex in the graph, which should replace the old vertex
     */

    template <typename Graph>
    void reconnect_edge (Graph& graph, typename boost::graph_traits <Graph>::edge_descriptor edge,
        typename boost::graph_traits <Graph>::vertex_descriptor old_vertex, typename boost::graph_traits <Graph>::vertex_descriptor new_vertex)
    {
      if (boost::source(edge, graph) == old_vertex)
        reconnect_edge_source(graph, edge, old_vertex, new_vertex);
      else
        reconnect_edge_target(graph, edge, old_vertex, new_vertex);
    }

    namespace detail {

      /**
       * BFS / DFS event filter to write vertices to an OutputIterator
       */

      template <typename OutputIterator, typename EventTag>
      struct vertex_writer
      {
        typedef EventTag event_filter;

        /**
         * Creates the filter
         *
         * @param iterator Iterator to write vertices to
         */

        vertex_writer (OutputIterator iterator) :
          _iterator(iterator)
        {

        }

        /**
         * Writes the given vertex into the OutputIterator
         *
         * @param vertex Event-vertex
         * @param graph
         */

        template <typename Vertex, typename Graph>
        void operator() (Vertex vertex, const Graph& graph)
        {
          *_iterator++ = vertex;
        }

        /**
         * @return The current iterator
         */

        OutputIterator iterator () const
        {
          return _iterator;
        }

      private:
        OutputIterator _iterator;
      };

      /**
       * Creates a vertex_write object.
       *
       * @param iterator The OutputIterator to write vertices to
       * @param tag Event tag for the event filter
       * @return A constructed vertex_write
       */

      template <typename OutputIterator, typename EventTag>
      inline vertex_writer <OutputIterator, EventTag> write_vertex (OutputIterator iterator, EventTag tag)
      {
        return vertex_writer <OutputIterator, EventTag> (iterator);
      }
    }

    /**
     * Tests, if a given graph is a path and stores the vertices in order into a given sequence.
     *
     * @param graph A given graph
     * @param index_map An index map for this graph
     * @param path The resulting sequence of nodes in order
     * @return true if and only if this graph is a path
     */

    template <typename Graph, typename IndexMap, typename VertexSequence>
    bool is_path (const Graph& graph, const IndexMap index_map, VertexSequence& path)
    {
      typedef boost::graph_traits <Graph> traits;
      typename traits::vertex_iterator vertex_iter, vertex_end;
      typename traits::vertex_descriptor start = traits::null_vertex();

      size_t path_length = 0;
      for (boost::tie(vertex_iter, vertex_end) = boost::vertices(graph); vertex_iter != vertex_end; ++vertex_iter)
      {
        size_t degree = boost::out_degree(*vertex_iter, graph);
        if (degree > 2)
          return false;
        else if (degree == 1)
          start = *vertex_iter;
        path_length += degree;
      }
      if (start == traits::null_vertex())
        return false;

      path_length = path_length / 2 + 1;

      size_t size_before = path.size();

      /// Call breadth first search
      boost::breadth_first_search(graph, start, boost::vertex_index_map(index_map).visitor(boost::make_bfs_visitor(detail::write_vertex(
          std::back_inserter(path), boost::on_discover_vertex()))));

      return path.size() - size_before == path_length;
    }

    /**
     * Tests, if a given graph is a path.
     *
     * @param graph A given graph
     * @param index_map An index map for this graph
     * @return true if and only if this graph is a path
     */

    template <typename Graph, typename IndexMap>
    bool is_path (const Graph& graph, const IndexMap index_map)
    {
      std::vector <typename boost::graph_traits <Graph>::vertex_descriptor> path;
      path.reserve(boost::num_vertices(graph));
      return is_path(graph, index_map, path);
    }

    /**
     * Tests, if a given graph is a path. The graph must have an index map attached.
     *
     * @param graph A given graph
     * @return true if and only if this graph is a path
     */

    template <typename Graph>
    bool is_path (const Graph& g)
    {
      return is_path(g, boost::get(boost::vertex_index, g));
    }

    /**
     * Tests, if a given graph is a star and stores the central vertex.
     *
     * @param graph A given graph
     * @param index_map An index map for this graph
     * @param vertex Returns the central vertex
     * @return true if and only if the graph is a star
     */

    template <typename Graph, typename IndexMap, typename Vertex>
    bool find_star_vertex (const Graph& graph, const IndexMap index_map, Vertex& vertex)
    {
      typedef boost::graph_traits <Graph> traits;
      typename traits::vertex_iterator vertex_iter, vertex_end;
      typename traits::edge_iterator edge_iter, edge_end;

      size_t edges = 0;
      for (boost::tie(edge_iter, edge_end) = boost::edges(graph); edge_iter != edge_end; ++edge_iter)
      {
        ++edges;
      }

      for (boost::tie(vertex_iter, vertex_end) = boost::vertices(graph); vertex_iter != vertex_end; ++vertex_iter)
      {
        if (boost::out_degree(*vertex_iter, graph) == edges)
        {
          vertex = *vertex_iter;
          return true;
        }
      }

      return false;
    }

    /**
     * Tests, if a given graph is a star and stores the central vertex. The graph must have an index map attached.
     *
     * @param graph A given graph
     * @param vertex Returns the central vertex
     * @return true if and only if the graph is a star
     */

    template <typename Graph, typename Vertex>
    bool find_star_vertex (const Graph& graph, Vertex& vertex)
    {
      return find_star_vertex(graph, boost::get(boost::vertex_index, graph), vertex);
    }

    /**
     * Tests, if a given graph is a star.
     *
     * @param graph A given graph
     * @param index_map An index map for this graph
     * @return true if and only if the graph is a star
     */

    template <typename Graph, typename IndexMap>
    bool is_star (const Graph& graph, const IndexMap index_map)
    {
      typename boost::graph_traits <Graph>::vertex_descriptor vertex;
      return find_star_vertex(graph, index_map, vertex);
    }

    /**
     * Tests, if a given graph is a star and stores the central vertex. The graph must have an index map attached.
     *
     * @param graph A given graph
     * @return true if and only if the graph is a star
     */

    template <typename Graph>
    bool is_star (const Graph& graph)
    {
      return is_star(graph, boost::get(boost::vertex_index, graph));
    }

    /**
     * Collects the non-isolated vertices of a given graph.
     *
     * @param graph A given graph
     * @param vertex_set A set to be filled with all vertices that are in an edge
     */

    template <typename Graph, typename VertexSet>
    void used_vertices (const Graph& graph, VertexSet& vertex_set)
    {
      typedef boost::graph_traits <Graph> traits;
      typename traits::edge_iterator edge_iter, edge_end;

      for (boost::tie(edge_iter, edge_end) = boost::edges(graph); edge_iter != edge_end; ++edge_iter)
      {
        vertex_set.insert(boost::source(*edge_iter, graph));
        vertex_set.insert(boost::target(*edge_iter, graph));
      }
    }

  }
}

#endif /* GRAPH_UTILS_HPP_ */
