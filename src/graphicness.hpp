/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef GRAPHICNESS_HPP_
#define GRAPHICNESS_HPP_

#include "boost/graph/bipartite.hpp"
#include <boost/graph/detail/set_adaptor.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/property_map/property_map.hpp>
#include <vector>

#include "matroid_graph.hpp"
#include "graph_utils.hpp"

#include <boost/graph/adjacency_list_io.hpp>

namespace tu
{

  /**
   * Finds a matroid element which is parallel to a row element in minor of a given matroid or
   * a unit vector in a column element in the minor.
   *
   * @param matroid The given matroid
   * @param matrix Representation matrix of the given matroid
   * @param minor_height Number of rows of the minor
   * @param minor_width Number of columns of the minor
   * @param row The element must be parallel to this row (along the minor)
   * @return The matroid element that was found to be parallel / unit vector
   */

  template <typename MatroidType, typename MatrixType>
  int find_parallel_to_row(const MatroidType& matroid, const MatrixType& matrix, const size_t minor_height, const size_t minor_width,
      const size_t row)
  {
    /// Look for unit vector
    size_t last = 0;
    size_t count = 0;
    for (size_t c = 0; c < minor_width; ++c)
    {
      if (matrix(row, c) != 0)
      {
        last = c;
        ++count;
      }
    }
    if (count == 1)
      return matroid.name2(last);

    /// Check for parallel
    for (size_t r = 0; r < minor_height; ++r)
    {
      bool same = true;
      for (size_t c = 0; c < minor_width; ++c)
      {
        if (matrix(r, c) != matrix(row, c))
        {
          same = false;
          break;
        }
      }
      if (same)
      {
        return matroid.name1(r);
      }
    }

    throw std::logic_error("find_parallel_to_row did not find any parallel / unit vector!");
  }

  /**
   * A filter for a graph which excludes some edges and a given vertex
   */

  template <typename Graph>
  struct articulation_edge_filter
  {
    typedef typename boost::graph_traits <Graph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits <Graph>::edge_descriptor edge_descriptor;
    typedef std::set <edge_descriptor> edge_set;

    /**
     * Constructs the filter.
     */

    articulation_edge_filter() :
      graph_(NULL), articulation_vertex_(NULL), evil_edges_(NULL)
    {
    }

    /**
     * Constructs the filter
     *
     * @param graph The original graph
     * @param articulation_vertex The excluded vertex
     * @param evil_edges The set of excluded edges
     */

    articulation_edge_filter(const Graph* graph, const vertex_descriptor* articulation_vertex, const edge_set* evil_edges) :
      graph_(graph), articulation_vertex_(articulation_vertex), evil_edges_(evil_edges)
    {
    }

    /**
     * Filter operator
     *
     * @param e The edge to be considered
     * @return true if and only if this edge is included in the filtered graph
     */

    template <typename Edge>
    bool operator()(const Edge& e) const
    {
      if (evil_edges_->find(e) != evil_edges_->end())
        return false;

      if (boost::source(e, *graph_) == *articulation_vertex_)
        return false;

      if (boost::target(e, *graph_) == *articulation_vertex_)
        return false;

      return true;
    }

  private:
    const Graph* graph_;
    const vertex_descriptor* articulation_vertex_;
    const edge_set* evil_edges_;
  };

  /**
   * Constructs an articulation edge filter for a given graph.
   *
   * @param graph The given graph
   * @param articulation_vertex Vertex to be excluded
   * @param evil_edges Edges to be excluded
   * @return The constructed filter
   */

  template <typename Graph, typename Vertex, typename EdgeSet>
  inline struct articulation_edge_filter <Graph> make_articulation_edge_filter(const Graph* graph, const Vertex* articulation_vertex,
      const EdgeSet* evil_edges)
  {
    return articulation_edge_filter <Graph> (graph, articulation_vertex, evil_edges);
  }

  /**
   * Tries to extend a 3-connected graph of a graphic minor or detects that
   * the extended minor is not graphic.
   *
   * @param graph The graph of the graphic minor
   * @param matroid The matroid
   * @param matrix Representation matrix of the matroid
   * @param minor_height Number of row elements of the minor
   * @param minor_width Number of column elements of the minor
   * @param extension_type Type of the extension
   * @return true if and only if the extension was possible.
   */

  template <typename MatroidType, typename MatrixType>
  bool extend_graph(matroid_graph& graph, const MatroidType& matroid, const MatrixType& matrix, const size_t minor_height, const size_t minor_width,
      const nested_minor_sequence::extension_type extension_type)
  {
    typedef boost::graph_traits <matroid_graph> traits;

    matroid_element_map element_map = boost::get(edge_matroid_element, graph);
    std::map <int, traits::edge_descriptor> reverse_element_map;

    typename traits::edge_iterator edge_iter, edge_end;
    for (boost::tie(edge_iter, edge_end) = boost::edges(graph); edge_iter != edge_end; ++edge_iter)
    {
      reverse_element_map[element_map[*edge_iter]] = *edge_iter;
    }

    /// Distinct extension cases
    switch (extension_type)
    {
    case nested_minor_sequence::TWO_ROWS_ONE_COLUMN:
    {
      int first_edge_element = find_parallel_to_row(matroid, matrix, minor_height, minor_width, minor_height);
      int second_edge_element = find_parallel_to_row(matroid, matrix, minor_height, minor_width, minor_height + 1);

      /// Find vertices of corresponding edges

      traits::vertex_descriptor first_vertex1 = boost::source(reverse_element_map[first_edge_element], graph);
      traits::vertex_descriptor first_vertex2 = boost::target(reverse_element_map[first_edge_element], graph);
      traits::vertex_descriptor second_vertex1 = boost::source(reverse_element_map[second_edge_element], graph);
      traits::vertex_descriptor second_vertex2 = boost::target(reverse_element_map[second_edge_element], graph);

      if (first_vertex1 == second_vertex2)
        std::swap(second_vertex1, second_vertex2);
      if (first_vertex2 == second_vertex1)
        std::swap(first_vertex1, first_vertex2);
      if (first_vertex2 == second_vertex2)
      {
        std::swap(first_vertex1, first_vertex2);
        std::swap(second_vertex1, second_vertex2);
      }
      if (first_vertex1 != second_vertex1)
        return false;

      /// Remove old edges

      boost::remove_edge(reverse_element_map[first_edge_element], graph);
      boost::remove_edge(reverse_element_map[second_edge_element], graph);

      /// Create new vertices

      traits::vertex_descriptor first_breaker = boost::add_vertex(graph);
      traits::vertex_descriptor second_breaker = boost::add_vertex(graph);

      /// Create new edges

      boost::add_edge(first_vertex2, first_breaker, first_edge_element, graph);
      boost::add_edge(first_breaker, first_vertex1, matroid.name1(minor_height), graph);
      boost::add_edge(second_vertex2, second_breaker, second_edge_element, graph);
      boost::add_edge(second_breaker, second_vertex1, matroid.name1(minor_height + 1), graph);
      boost::add_edge(first_breaker, second_breaker, matroid.name2(minor_width), graph);

      return true;
    }
    break;
    case nested_minor_sequence::ONE_ROW_TWO_COLUMNS:
    {
      int first_edge_element = find_parallel_to_row(make_transposed_matroid(matroid), make_transposed_matrix(matrix), minor_width, minor_height,
          minor_width);
      int second_edge_element = find_parallel_to_row(make_transposed_matroid(matroid), make_transposed_matrix(matrix), minor_width, minor_height,
          minor_width + 1);

      /// Find vertices of corresponding edges

      traits::vertex_descriptor first_vertex1 = boost::source(reverse_element_map[first_edge_element], graph);
      traits::vertex_descriptor first_vertex2 = boost::target(reverse_element_map[first_edge_element], graph);
      traits::vertex_descriptor second_vertex1 = boost::source(reverse_element_map[second_edge_element], graph);
      traits::vertex_descriptor second_vertex2 = boost::target(reverse_element_map[second_edge_element], graph);

      if (first_vertex1 == second_vertex2)
        std::swap(second_vertex1, second_vertex2);
      if (first_vertex2 == second_vertex1)
        std::swap(first_vertex1, first_vertex2);
      if (first_vertex2 == second_vertex2)
      {
        std::swap(first_vertex1, first_vertex2);
        std::swap(second_vertex1, second_vertex2);
      }
      if (first_vertex1 != second_vertex1)
        return false;

      /// Create new vertex

      traits::vertex_descriptor additional_vertex = boost::add_vertex(graph);

      /// Create new edges

      boost::add_edge(first_vertex2, additional_vertex, matroid.name2(minor_width), graph);
      boost::add_edge(second_vertex2, additional_vertex, matroid.name2(minor_width + 1), graph);
      boost::add_edge(first_vertex1, additional_vertex, matroid.name1(minor_height), graph);

      return true;
    }
    break;
    case nested_minor_sequence::ONE_ROW_ONE_COLUMN:
    {
      int row_edge_element = find_parallel_to_row(matroid, matrix, minor_height, minor_width, minor_height);
      int column_edge_element = find_parallel_to_row(make_transposed_matroid(matroid), make_transposed_matrix(matrix), minor_width, minor_height,
          minor_width);

      /// Find vertices of corresponding edges

      traits::vertex_descriptor row_vertex1 = boost::source(reverse_element_map[row_edge_element], graph);
      traits::vertex_descriptor row_vertex2 = boost::target(reverse_element_map[row_edge_element], graph);
      traits::vertex_descriptor column_vertex1 = boost::source(reverse_element_map[column_edge_element], graph);
      traits::vertex_descriptor column_vertex2 = boost::target(reverse_element_map[column_edge_element], graph);

      if (row_vertex1 == column_vertex2)
        std::swap(column_vertex1, column_vertex2);
      if (row_vertex2 == column_vertex1)
        std::swap(row_vertex1, row_vertex2);
      if (row_vertex2 == column_vertex2)
      {
        std::swap(row_vertex1, row_vertex2);
        std::swap(column_vertex1, column_vertex2);
      }
      if (row_vertex1 != column_vertex1)
        return false;

      /// Remove old edges

      boost::remove_edge(reverse_element_map[row_edge_element], graph);

      /// Create new vertex

      traits::vertex_descriptor row_breaker = boost::add_vertex(graph);

      /// Create new edges

      boost::add_edge(row_vertex2, row_breaker, row_edge_element, graph);
      boost::add_edge(row_breaker, row_vertex1, matroid.name1(minor_height), graph);
      boost::add_edge(row_breaker, column_vertex2, matroid.name2(minor_width), graph);

      return true;
    }
    break;
    case nested_minor_sequence::ONE_COLUMN:
    {
      typedef std::set <traits::edge_descriptor> edge_set_t;
      typedef boost::filtered_graph <matroid_graph, boost::is_in_subset <edge_set_t> > edge_subset_graph_t;

      edge_set_t edge_set;

      for (size_t r = 0; r < minor_height; ++r)
      {
        if (matrix(r, minor_width) == 1)
          edge_set.insert(reverse_element_map[matroid.name1(r)]);
      }

      /// Check if edges form a path
      std::vector <traits::vertex_descriptor> path;

      edge_subset_graph_t edge_subset_graph(graph, boost::is_in_subset <edge_set_t>(edge_set));
      if (!tu::util::is_path(edge_subset_graph, boost::get(boost::vertex_index, edge_subset_graph), path))
      {
        return false;
      }

      boost::add_edge(path[0], path[path.size() - 1], matroid.name2(minor_width), graph);
      return true;
    }
    break;
    case nested_minor_sequence::ONE_ROW:
    {
      typedef traits::edge_descriptor edge_t;
      typedef traits::vertex_descriptor vertex_t;
      typedef std::set <vertex_t> vertex_set;
      typedef std::set <edge_t> edge_set;
      typedef std::vector <vertex_t> vertex_vector_t;
      typedef std::vector <edge_t> edge_vector_t;

      typedef boost::vec_adj_list_vertex_id_map <boost::no_property, unsigned int> IndexMap;
      IndexMap index_map = boost::get(boost::vertex_index, graph);

      /// Collect all 1-edges

      edge_set one_edges;
      for (size_t c = 0; c < minor_width; ++c)
      {
        if (matrix(minor_height, c) != 0)
          one_edges.insert(reverse_element_map[matroid.name2(c)]);
      }

      traits::vertex_descriptor the_vertex = traits::null_vertex();
      if (tu::util::find_star_vertex(boost::make_filtered_graph(graph, boost::is_in_subset <edge_set>(one_edges)), the_vertex))
      {
        /// Create the new vertex and connect it with the_vertex.

        traits::vertex_descriptor new_vertex = boost::add_vertex(graph);
        boost::add_edge(the_vertex, new_vertex, matroid.name1(minor_height), graph);

        /// Reconnect the 1-edges

        for (edge_set::const_iterator edge_iter = one_edges.begin(); edge_iter != one_edges.end(); ++edge_iter)
        {
          traits::edge_descriptor edge = *edge_iter;
          traits::vertex_descriptor other_vertex = boost::source(edge, graph);
          if (other_vertex == the_vertex)
            other_vertex = boost::target(edge, graph);
          int name = element_map[edge];
          boost::remove_edge(the_vertex, other_vertex, graph);
          boost::add_edge(new_vertex, other_vertex, name, graph);
        }

        return true;
      }

      /// Count for each vertex the number of paths that use it

      std::vector <size_t> common_vertex_count(boost::num_vertices(graph), 0);
      for (size_t c = 0; c < minor_width; ++c)
      {
        if (matrix(minor_height, c) == 0)
          continue;

        edge_set edges;
        for (size_t r = 0; r < minor_height; ++r)
        {
          if (matrix(r, c) == 1)
            edges.insert(reverse_element_map[matroid.name1(r)]);
        }

        vertex_set vertices;
        tu::util::used_vertices(boost::make_filtered_graph(graph, boost::is_in_subset <edge_set>(edges)), vertices);

        for (typename vertex_set::const_iterator iter = vertices.begin(); iter != vertices.end(); ++iter)
        {
          common_vertex_count[boost::get(index_map, *iter)]++;
        }
      }

      /// Find articulation points for the graph without 1-edges

      vertex_vector_t articulation_points;
      boost::articulation_points(boost::make_filtered_graph(graph, boost::is_not_in_subset <edge_set>(one_edges)), std::back_inserter(
          articulation_points));

      the_vertex = traits::null_vertex();
      for (vertex_vector_t::const_iterator iter = articulation_points.begin(); iter != articulation_points.end(); ++iter)
      {
        if (common_vertex_count[boost::get(index_map, *iter)] == one_edges.size())
        {
          if (the_vertex == traits::null_vertex())
            the_vertex = *iter;
          else
            return false;
        }
      }

      if (the_vertex == traits::null_vertex())
        return false;

      /// Filter the unique articulation point and the one-edges and check the remaining graph

      vertex_set articulation_set;
      articulation_set.insert(the_vertex);
      std::vector <size_t> component_vector(boost::num_vertices(graph));

      size_t num_components = boost::connected_components(boost::make_filtered_graph(graph, make_articulation_edge_filter(&graph, &the_vertex,
          &one_edges), boost::is_not_in_subset <vertex_set>(articulation_set)),
          boost::make_iterator_property_map(component_vector.begin(), index_map));

      /// We should really have articulation point + 2 further components
      assert (num_components >= 2);

      typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> component_graph_t;
      component_graph_t component_graph(num_components);

      for (typename edge_set::const_iterator iter = one_edges.begin(); iter != one_edges.end(); ++iter)
      {
        typename boost::graph_traits <component_graph_t>::vertex_descriptor source, target;
        source = boost::source(*iter, graph);
        target = boost::target(*iter, graph);
        if (source == the_vertex || target == the_vertex)
          continue;

        size_t source_component = component_vector[boost::get(index_map, source)];
        size_t target_component = component_vector[boost::get(index_map, target)];

        if (source_component == target_component)
        {
          /// There cannot be a 1-edge inside one component.
          return false;
        }

        boost::add_edge(boost::vertex(source_component, component_graph), boost::vertex(target_component, component_graph), component_graph);
      }

      boost::one_bit_color_map <boost::vec_adj_list_vertex_id_map <boost::no_property, unsigned int> > bipartition(num_components, boost::get(
          boost::vertex_index, component_graph));

      if (!boost::is_bipartite(component_graph, boost::get(boost::vertex_index, component_graph), bipartition))
      {
        return false;
      }

      for (size_t i = 0; i < component_vector.size(); ++i)
      {
        if (boost::vertex(i, graph) == the_vertex)
          continue;
      }

      vertex_t new_vertex = boost::add_vertex(graph);

      edge_vector_t reconnect_edges;
      typename traits::out_edge_iterator out_edge_iter, out_edge_end;
      for (boost::tie(out_edge_iter, out_edge_end) = boost::incident_edges(the_vertex, graph); out_edge_iter != out_edge_end; ++out_edge_iter)
      {
        vertex_t incident_vertex = boost::target(*out_edge_iter, graph);

        bool reconnect = boost::get(bipartition, boost::vertex(component_vector[boost::get(index_map, incident_vertex)], component_graph))
            != boost::one_bit_white;

        if (one_edges.find(*out_edge_iter) != one_edges.end())
          reconnect = !reconnect;

        if (reconnect)
          reconnect_edges.push_back(*out_edge_iter);
      }

      for (typename edge_vector_t::const_iterator iter = reconnect_edges.begin(); iter != reconnect_edges.end(); ++iter)
      {
        util::reconnect_edge(graph, *iter, the_vertex, new_vertex);
      }

      boost::add_edge(the_vertex, new_vertex, matroid.name1(minor_height), graph);

      return true;

    }
    break;
    default:
      throw std::logic_error("Unknown extension in graphicness test.");
    }
  }

  /**
   * Either constructs a graph whose forest matroid is the given matroid or detects
   * that the given matroid is non-graphic.
   *
   * @param matroid The given matroid
   * @param matrix Representation matrix of the given matroid
   * @param nested_minors Sequence of nested minors
   * @return The constructed graph or NULL if the matroid is not graphic.
   */

  template <typename MatroidType, typename MatrixType, typename NestedMinorSequenceType>
  matroid_graph* construct_matroid_graph(const MatroidType& matroid, const MatrixType& matrix, const NestedMinorSequenceType& nested_minors)
  {
    typedef boost::graph_traits <matroid_graph>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits <matroid_graph>::edge_descriptor edge_descriptor;

    /// Initialize W3

    matroid_graph* graph = new matroid_graph(4);

    vertex_descriptor center_vertex = boost::vertex(0, *graph);
    vertex_descriptor border_vertex1 = boost::vertex(1, *graph);
    vertex_descriptor border_vertex2 = boost::vertex(2, *graph);
    vertex_descriptor border_vertex3 = boost::vertex(3, *graph);

    edge_descriptor spoke1 = boost::add_edge(center_vertex, border_vertex1, matroid.name1(0), *graph).first;
    edge_descriptor spoke2 = boost::add_edge(center_vertex, border_vertex2, matroid.name1(1), *graph).first;
    edge_descriptor spoke3 = boost::add_edge(center_vertex, border_vertex3, matroid.name2(2), *graph).first;
    edge_descriptor rim12 = boost::add_edge(border_vertex1, border_vertex2, matroid.name2(0), *graph).first;
    edge_descriptor rim13 = boost::add_edge(border_vertex1, border_vertex3, matroid.name2(1), *graph).first;
    edge_descriptor rim23 = boost::add_edge(border_vertex2, border_vertex3, matroid.name1(2), *graph).first;

    size_t minor_height = 3;
    size_t minor_width = 3;

    for (size_t i = 0; i < nested_minors.size(); ++i)
    {
      if (!extend_graph(*graph, matroid, matrix, minor_height, minor_width, nested_minors.get_extension(i)))
      {
        delete graph;
        return NULL;
      }

      minor_height += nested_minors.get_extension_height(i);
      minor_width += nested_minors.get_extension_width(i);
    }

    return graph;
  }

}

#endif /* GRAPHICNESS_HPP_ */
