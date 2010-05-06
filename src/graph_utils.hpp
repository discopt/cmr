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

    template <typename Graph>
    void reconnect_edge (Graph& g, typename boost::graph_traits <Graph>::edge_descriptor edge,
        typename boost::graph_traits <Graph>::vertex_descriptor old_vertex, typename boost::graph_traits <Graph>::vertex_descriptor new_vertex)
    {
      if (boost::source(edge, g) == old_vertex)
        reconnect_edge_source(g, edge, old_vertex, new_vertex);
      else
        reconnect_edge_target(g, edge, old_vertex, new_vertex);
    }

    namespace detail {

      template <typename OutputIterator, typename EventTag>
      struct vertex_writer
      {
        typedef EventTag event_filter;

        vertex_writer (OutputIterator iterator) :
          iterator_(iterator)
        {

        }

        template <typename Vertex, typename Graph>
        void operator() (Vertex vertex, const Graph& g)
        {
          *iterator_ = vertex;
        }

        OutputIterator iterator () const
        {
          return iterator_;
        }

      private:
        OutputIterator iterator_;
      };

      template <typename OutputIterator, typename EventTag>
      inline vertex_writer <OutputIterator, EventTag> write_vertex (OutputIterator iterator, EventTag tag)
      {
        return vertex_writer <OutputIterator, EventTag> (iterator);
      }
    }

    template <typename Graph, typename IndexMap, typename VertexSequence>
    bool is_path (const Graph& g, const IndexMap index_map, VertexSequence& path)
    {
      typedef boost::graph_traits <Graph> traits;
      typename traits::vertex_iterator vertex_iter, vertex_end;
      typename traits::vertex_descriptor start = traits::null_vertex();

      size_t path_length = 0;
      for (boost::tie(vertex_iter, vertex_end) = boost::vertices(g); vertex_iter != vertex_end; ++vertex_iter)
      {
        size_t degree = boost::out_degree(*vertex_iter, g);
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

      boost::breadth_first_search(g, start, boost::vertex_index_map(index_map).visitor(boost::make_bfs_visitor(detail::write_vertex(
          std::back_inserter(path), boost::on_discover_vertex()))));

      return path.size() - size_before == path_length;
    }

    template <typename Graph, typename IndexMap>
    bool is_path (const Graph& g, const IndexMap index_map)
    {
      std::vector <typename boost::graph_traits <Graph>::vertex_descriptor> path;
      path.reserve(boost::num_vertices(g));
      return is_path(g, index_map, path);
    }

    template <typename Graph>
    bool is_path (const Graph& g)
    {
      return is_path(g, boost::get(boost::vertex_index, g));
    }

    template <typename Graph, typename IndexMap, typename Vertex>
    bool find_star_vertex (const Graph& g, const IndexMap index_map, Vertex& vertex)
    {
      typedef boost::graph_traits <Graph> traits;
      typename traits::vertex_iterator vertex_iter, vertex_end;
      typename traits::edge_iterator edge_iter, edge_end;

      size_t edges = 0;
      for (boost::tie(edge_iter, edge_end) = boost::edges(g); edge_iter != edge_end; ++edge_iter)
      {
        ++edges;
      }

      for (boost::tie(vertex_iter, vertex_end) = boost::vertices(g); vertex_iter != vertex_end; ++vertex_iter)
      {
        if (boost::out_degree(*vertex_iter, g) == edges)
        {
          vertex = *vertex_iter;
          return true;
        }
      }

      return false;
    }

    template <typename Graph, typename Vertex>
    bool find_star_vertex (const Graph& g, Vertex& vertex)
    {
      return find_star_vertex(g, boost::get(boost::vertex_index, g), vertex);
    }

    template <typename Graph, typename IndexMap>
    bool is_star (const Graph& g, const IndexMap index_map)
    {
      typename boost::graph_traits <Graph>::vertex_descriptor vertex;
      return find_star_vertex(g, index_map, vertex);
    }

    template <typename Graph>
    bool is_star (const Graph& g)
    {
      return is_star(g, boost::get(boost::vertex_index, g));
    }

    template <typename Graph, typename VertexSet>
    void used_vertices (const Graph& g, VertexSet& vertex_set)
    {
      typedef boost::graph_traits <Graph> traits;
      typename traits::edge_iterator edge_iter, edge_end;

      for (boost::tie(edge_iter, edge_end) = boost::edges(g); edge_iter != edge_end; ++edge_iter)
      {
        vertex_set.insert(boost::source(*edge_iter, g));
        vertex_set.insert(boost::target(*edge_iter, g));
      }
    }

  }
}

#endif /* GRAPH_UTILS_HPP_ */
