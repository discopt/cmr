
//          Copyright Matthias Walter 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef GRAPHICNESS_HPP_
#define GRAPHICNESS_HPP_

#include <boost/graph/detail/set_adaptor.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/property_map/property_map.hpp>
#include "boost/graph/bipartite.hpp"
#include <vector>

#include "matroid_graph.hpp"
#include "graph_utils.hpp"

#include <boost/graph/adjacency_list_io.hpp>

namespace tu {

  template <typename MatroidType, typename MatrixType>
  int find_parallel_to_row (const MatroidType& matroid, const MatrixType& matrix, const size_t minor_height, const size_t minor_width,
      const size_t row)
  {
    /// Look for unit vector
    size_t last = 0;
    size_t count = 0;
    for (size_t c = 0; c < minor_width; ++c)
    {
      if (matrix (row, c) != 0)
      {
        last = c;
        ++count;
      }
    }
    if (count == 1)
      return matroid.name2 (last);

    /// Check for parallel
    for (size_t r = 0; r < minor_height; ++r)
    {
      bool same = true;
      for (size_t c = 0; c < minor_width; ++c)
      {
        if (matrix (r, c) != matrix (row, c))
        {
          same = false;
          break;
        }
      }
      if (same)
      {
        return matroid.name1 (r);
      }
    }

    throw std::logic_error ("find_parallel_to_row did not find any parallel / unit vector!");
  }

  template <typename Graph>
  struct articulation_edge_filter
  {
    typedef typename boost::graph_traits <Graph>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits <Graph>::edge_descriptor edge_descriptor;
    typedef std::set <edge_descriptor> edge_set;

    articulation_edge_filter () :
      graph_ (NULL), articulation_vertex_ (NULL), evil_edges_ (NULL)
    {
    }

    articulation_edge_filter (const Graph* graph, const vertex_descriptor* articulation_vertex, const edge_set* evil_edges) :
      graph_ (graph), articulation_vertex_ (articulation_vertex), evil_edges_ (evil_edges)
    {
    }

    template <typename Edge>
    bool operator() (const Edge& e) const
    {
      if (evil_edges_->find (e) != evil_edges_->end ())
        return false;

      if (boost::source (e, *graph_) == *articulation_vertex_)
        return false;

      if (boost::target (e, *graph_) == *articulation_vertex_)
        return false;

      return true;
    }

  private:
    const Graph* graph_;
    const vertex_descriptor* articulation_vertex_;
    const edge_set* evil_edges_;
  };

  template <typename Graph, typename Vertex, typename EdgeSet>
  inline struct articulation_edge_filter <Graph> make_articulation_edge_filter (const Graph* graph, const Vertex* articulation_vertex,
      const EdgeSet* evil_edges)
  {
    return articulation_edge_filter <Graph> (graph, articulation_vertex, evil_edges);
  }

  template <typename MatroidType, typename MatrixType>
  bool extend_graph (matroid_graph& graph, const MatroidType& matroid, const MatrixType& matrix, const size_t minor_height, const size_t minor_width,
      const nested_minor_sequence::extension_type extension_type)
  {
//    matroid_print (matroid, matrix);
//    std::cout << "Extending the following graph with an extension of ";
//    switch (extension_type)
//    {
//    case nested_minor_sequence::ONE_ROW:
//      std::cout << "1-row-type";
//    break;
//    case nested_minor_sequence::ONE_ROW_TWO_COLUMNS:
//      std::cout << "1-row-2-columns-type";
//    break;
//    case nested_minor_sequence::ONE_ROW_ONE_COLUMN:
//      std::cout << "1-row-1-column-type";
//    break;
//    case nested_minor_sequence::TWO_ROWS_ONE_COLUMN:
//      std::cout << "2-rows-1-column-type";
//    break;
//    case nested_minor_sequence::ONE_COLUMN:
//      std::cout << "1-column-type";
//    break;
//    default:
//      std::cout << "<unknown>";
//    break;
//    }
//    std::cout << " and current minor: " << minor_height << " x " << minor_width << std::endl;
//    std::cout << graph << std::endl;

    typedef boost::graph_traits <matroid_graph> traits;

    matroid_element_map element_map = boost::get (edge_matroid_element, graph);
    std::map <int, traits::edge_descriptor> reverse_element_map;

    typename traits::edge_iterator edge_iter, edge_end;
    for (boost::tie (edge_iter, edge_end) = boost::edges (graph); edge_iter != edge_end; ++edge_iter)
    {
      reverse_element_map[element_map[*edge_iter]] = *edge_iter;
    }

    switch (extension_type)
    {
    case nested_minor_sequence::TWO_ROWS_ONE_COLUMN:
    {
      int first_edge_element = find_parallel_to_row (matroid, matrix, minor_height, minor_width, minor_height);
      int second_edge_element = find_parallel_to_row (matroid, matrix, minor_height, minor_width, minor_height + 1);

      /// Find vertices of corresponding edges

      traits::vertex_descriptor first_vertex1 = boost::source (reverse_element_map[first_edge_element], graph);
      traits::vertex_descriptor first_vertex2 = boost::target (reverse_element_map[first_edge_element], graph);
      traits::vertex_descriptor second_vertex1 = boost::source (reverse_element_map[second_edge_element], graph);
      traits::vertex_descriptor second_vertex2 = boost::target (reverse_element_map[second_edge_element], graph);

      if (first_vertex1 == second_vertex2)
        std::swap (second_vertex1, second_vertex2);
      if (first_vertex2 == second_vertex1)
        std::swap (first_vertex1, first_vertex2);
      if (first_vertex2 == second_vertex2)
      {
        std::swap (first_vertex1, first_vertex2);
        std::swap (second_vertex1, second_vertex2);
      }
      if (first_vertex1 != second_vertex1)
        return false;

      /// Remove old edges

      boost::remove_edge (reverse_element_map[first_edge_element], graph);
      boost::remove_edge (reverse_element_map[second_edge_element], graph);

      /// Create new vertices

      traits::vertex_descriptor first_breaker = boost::add_vertex (graph);
      traits::vertex_descriptor second_breaker = boost::add_vertex (graph);

      /// Create new edges

      boost::add_edge (first_vertex2, first_breaker, first_edge_element, graph);
      boost::add_edge (first_breaker, first_vertex1, matroid.name1 (minor_height), graph);
      boost::add_edge (second_vertex2, second_breaker, second_edge_element, graph);
      boost::add_edge (second_breaker, second_vertex1, matroid.name1 (minor_height + 1), graph);
      boost::add_edge (first_breaker, second_breaker, matroid.name2 (minor_width), graph);

      return true;
    }
    break;
    case nested_minor_sequence::ONE_ROW_TWO_COLUMNS:
    {
      int first_edge_element = find_parallel_to_row (view_matroid_transposed (matroid), view_matrix_transposed (matrix), minor_width, minor_height,
          minor_width);
      int second_edge_element = find_parallel_to_row (view_matroid_transposed (matroid), view_matrix_transposed (matrix), minor_width, minor_height,
          minor_width + 1);

      /// Find vertices of corresponding edges

      traits::vertex_descriptor first_vertex1 = boost::source (reverse_element_map[first_edge_element], graph);
      traits::vertex_descriptor first_vertex2 = boost::target (reverse_element_map[first_edge_element], graph);
      traits::vertex_descriptor second_vertex1 = boost::source (reverse_element_map[second_edge_element], graph);
      traits::vertex_descriptor second_vertex2 = boost::target (reverse_element_map[second_edge_element], graph);

      if (first_vertex1 == second_vertex2)
        std::swap (second_vertex1, second_vertex2);
      if (first_vertex2 == second_vertex1)
        std::swap (first_vertex1, first_vertex2);
      if (first_vertex2 == second_vertex2)
      {
        std::swap (first_vertex1, first_vertex2);
        std::swap (second_vertex1, second_vertex2);
      }
      if (first_vertex1 != second_vertex1)
        return false;

      /// Create new vertex

      traits::vertex_descriptor additional_vertex = boost::add_vertex (graph);

      /// Create new edges

      boost::add_edge (first_vertex2, additional_vertex, matroid.name2 (minor_width), graph);
      boost::add_edge (second_vertex2, additional_vertex, matroid.name2 (minor_width + 1), graph);
      boost::add_edge (first_vertex1, additional_vertex, matroid.name1 (minor_height), graph);

      return true;
    }
    break;
    case nested_minor_sequence::ONE_ROW_ONE_COLUMN:
    {
      int row_edge_element = find_parallel_to_row (matroid, matrix, minor_height, minor_width, minor_height);
      int column_edge_element = find_parallel_to_row (view_matroid_transposed (matroid), view_matrix_transposed (matrix), minor_width, minor_height,
          minor_width);

      /// Find vertices of corresponding edges

      traits::vertex_descriptor row_vertex1 = boost::source (reverse_element_map[row_edge_element], graph);
      traits::vertex_descriptor row_vertex2 = boost::target (reverse_element_map[row_edge_element], graph);
      traits::vertex_descriptor column_vertex1 = boost::source (reverse_element_map[column_edge_element], graph);
      traits::vertex_descriptor column_vertex2 = boost::target (reverse_element_map[column_edge_element], graph);

      if (row_vertex1 == column_vertex2)
        std::swap (column_vertex1, column_vertex2);
      if (row_vertex2 == column_vertex1)
        std::swap (row_vertex1, row_vertex2);
      if (row_vertex2 == column_vertex2)
      {
        std::swap (row_vertex1, row_vertex2);
        std::swap (column_vertex1, column_vertex2);
      }
      if (row_vertex1 != column_vertex1)
        return false;

      /// Remove old edges

      boost::remove_edge (reverse_element_map[row_edge_element], graph);

      /// Create new vertex

      traits::vertex_descriptor row_breaker = boost::add_vertex (graph);

      /// Create new edges

      boost::add_edge (row_vertex2, row_breaker, row_edge_element, graph);
      boost::add_edge (row_breaker, row_vertex1, matroid.name1 (minor_height), graph);
      boost::add_edge (row_breaker, column_vertex2, matroid.name2 (minor_width), graph);

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
        if (matrix (r, minor_width) == 1)
          edge_set.insert (reverse_element_map[matroid.name1 (r)]);
      }

      /// Check if edges form a path
      std::vector <traits::vertex_descriptor> path;

      edge_subset_graph_t edge_subset_graph (graph, boost::is_in_subset <edge_set_t> (edge_set));
      if (!tu::util::is_path (edge_subset_graph, boost::get (boost::vertex_index, edge_subset_graph), path))
      {
        //        std::cout << "1-edges did not form a path" << std::endl;
        return false;
      }

      //      std::cout << "1-edges did form a path:";
      //      for (size_t i = 0; i < path.size (); ++i)
      //        std::cout << " " << path[i] << std::endl;

      boost::add_edge (path[0], path[path.size () - 1], matroid.name2 (minor_width), graph);
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
      IndexMap index_map = boost::get (boost::vertex_index, graph);

      /// Collect all 1-edges

      edge_set one_edges;
      for (size_t c = 0; c < minor_width; ++c)
      {
        if (matrix (minor_height, c) != 0)
          one_edges.insert (reverse_element_map[matroid.name2 (c)]);
      }

      traits::vertex_descriptor the_vertex = traits::null_vertex ();
      if (tu::util::find_star_vertex (boost::make_filtered_graph (graph, boost::is_in_subset <edge_set> (one_edges)), the_vertex))
      {
        //        std::cout << "More simple case on vertex " << the_vertex << std::endl;

        /// Create the new vertex and connect it with the_vertex.

        traits::vertex_descriptor new_vertex = boost::add_vertex (graph);
        boost::add_edge (the_vertex, new_vertex, matroid.name1 (minor_height), graph);

        /// Reconnect the 1-edges

        for (edge_set::const_iterator edge_iter = one_edges.begin (); edge_iter != one_edges.end (); ++edge_iter)
        {
          traits::edge_descriptor edge = *edge_iter;
          traits::vertex_descriptor other_vertex = boost::source (edge, graph);
          if (other_vertex == the_vertex)
            other_vertex = boost::target (edge, graph);
          int name = element_map[edge];
          boost::remove_edge (the_vertex, other_vertex, graph);
          boost::add_edge (new_vertex, other_vertex, name, graph);
        }

        return true;
      }

      /// Count for each vertex the number of paths that use it

      std::vector <size_t> common_vertex_count (boost::num_vertices (graph), 0);
      for (size_t c = 0; c < minor_width; ++c)
      {
        if (matrix (minor_height, c) == 0)
          continue;

        edge_set edges;
        for (size_t r = 0; r < minor_height; ++r)
        {
          if (matrix (r, c) == 1)
            edges.insert (reverse_element_map[matroid.name1 (r)]);
        }

        vertex_set vertices;
        tu::util::used_vertices (boost::make_filtered_graph (graph, boost::is_in_subset <edge_set> (edges)), vertices);

        for (typename vertex_set::const_iterator iter = vertices.begin (); iter != vertices.end (); ++iter)
        {
          common_vertex_count[boost::get (index_map, *iter)]++;
        }
      }

      /// Find articulation points for the graph without 1-edges

//      std::cout << "1-edges: ";
//      std::copy (one_edges.begin (), one_edges.end (), std::ostream_iterator <edge_t> (std::cout, " "));
//      std::cout << std::endl;

      vertex_vector_t articulation_points;
      boost::articulation_points (boost::make_filtered_graph (graph, boost::is_not_in_subset <edge_set> (one_edges)), std::back_inserter (
          articulation_points));

//      std::cout << "APs: ";
//      std::copy (articulation_points.begin (), articulation_points.end (), std::ostream_iterator <vertex_t> (std::cout, " "));
//      std::cout << std::endl;

      the_vertex = traits::null_vertex ();
      for (vertex_vector_t::const_iterator iter = articulation_points.begin (); iter != articulation_points.end (); ++iter)
      {
        if (common_vertex_count[boost::get (index_map, *iter)] == one_edges.size ())
        {
//          std::cout << "candidate vertex (AP and common to 1-edges): " << *iter << std::endl;
          if (the_vertex == traits::null_vertex ())
            the_vertex = *iter;
          else
          {
            //                        std::cout << "Found 2 articulation points" << std::endl;
            return false;
          }
        }
      }

      if (the_vertex == traits::null_vertex ())
        return false;

      /// Filter the unique articulation point and the one-edges and check the remaining graph

      vertex_set articulation_set;
      articulation_set.insert (the_vertex);
      std::vector <size_t> component_vector (boost::num_vertices (graph));

      size_t num_components = boost::connected_components (boost::make_filtered_graph (graph, make_articulation_edge_filter (&graph, &the_vertex,
          &one_edges), boost::is_not_in_subset <vertex_set> (articulation_set)), boost::make_iterator_property_map (component_vector.begin (),
          index_map));

//      for (size_t i = 0; i < component_vector.size (); ++i)
//      {
//        std::cout << "compo[" << i << "] = " << component_vector[i] << std::endl;
//      }

      /// We should really have articulation point + 2 further components
      assert (num_components >= 2);

      typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> component_graph_t;
      component_graph_t component_graph (num_components);

      for (typename edge_set::const_iterator iter = one_edges.begin (); iter != one_edges.end (); ++iter)
      {
        typename boost::graph_traits <component_graph_t>::vertex_descriptor source, target;
        source = boost::source (*iter, graph);
        target = boost::target (*iter, graph);
        if (source == the_vertex || target == the_vertex)
          continue;

        size_t source_component = component_vector[boost::get (index_map, source)];
        size_t target_component = component_vector[boost::get (index_map, target)];

        if (source_component == target_component)
        {
          /// There cannot be a 1-edge inside one component.
          return false;
        }

        boost::add_edge (boost::vertex (source_component, component_graph), boost::vertex (target_component, component_graph), component_graph);
      }

      boost::one_bit_color_map <boost::vec_adj_list_vertex_id_map <boost::no_property, unsigned int> > bipartition (num_components, boost::get (
          boost::vertex_index, component_graph));

      if (!boost::is_bipartite (component_graph, boost::get (boost::vertex_index, component_graph), bipartition))
      {
        //        std::cout << "Component graph is not bipartite!" << std::endl;

        return false;
      }

      for (size_t i = 0; i < component_vector.size (); ++i)
      {
        if (boost::vertex (i, graph) == the_vertex)
          continue;

//        std::cout << "Vertex " << i << " is in component " << component_vector[i] << " and colored ";
//        if (boost::get (bipartition, boost::vertex (component_vector[i], component_graph)) == boost::one_bit_white)
//          std::cout << "white";
//        else
//          std::cout << "black";
//        std::cout << std::endl;
      }

      vertex_t new_vertex = boost::add_vertex (graph);

//      std::cout << "the_vertex = " << the_vertex << ", new_vertex = " << new_vertex << std::endl;

      edge_vector_t reconnect_edges;
      typename traits::out_edge_iterator out_edge_iter, out_edge_end;
      for (boost::tie (out_edge_iter, out_edge_end) = boost::incident_edges (the_vertex, graph); out_edge_iter != out_edge_end; ++out_edge_iter)
      {
        //        std::cout << "\nChecking edge " << *out_edge_iter << std::endl;

        vertex_t incident_vertex = boost::target (*out_edge_iter, graph);

        //        std::cout << "incident vertex = " << incident_vertex << std::endl;

        bool reconnect = boost::get (bipartition, boost::vertex (component_vector[boost::get (index_map, incident_vertex)], component_graph))
            != boost::one_bit_white;

        //        std::cout << "is " << (reconnect ? "white" : "black") << std::endl;

        if (one_edges.find (*out_edge_iter) != one_edges.end ())
        {
          //          std::cout << "edge is 1-edge" << std::endl;
          reconnect = !reconnect;
        }

        if (reconnect)
        {
          //          std::cout << "Adding to reconnect list" << std::endl;
          reconnect_edges.push_back (*out_edge_iter);
        }
      }

      for (typename edge_vector_t::const_iterator iter = reconnect_edges.begin (); iter != reconnect_edges.end (); ++iter)
      {
        util::reconnect_edge (graph, *iter, the_vertex, new_vertex);
      }

      boost::add_edge (the_vertex, new_vertex, matroid.name1 (minor_height), graph);

//      std::cout << "Resulting graph:\n" << graph << std::endl;

      return true;

    }
    break;
    default:
      throw std::logic_error ("Unknown extension in graphicness test.");
    }
  }

  template <typename MatroidType, typename MatrixType, typename NestedMinorSequenceType>
  matroid_graph* construct_matroid_graph (const MatroidType& matroid, const MatrixType& matrix, const NestedMinorSequenceType& nested_minors)
  {
    //    std::cout << "  [GRAPHICNESS TEST]\n";
    //    matroid_print (matroid, matrix);

    typedef boost::graph_traits <matroid_graph>::vertex_descriptor vertex_descriptor;
    typedef boost::graph_traits <matroid_graph>::edge_descriptor edge_descriptor;

    /// Initialize W3

    matroid_graph* g = new matroid_graph (4);

    vertex_descriptor center_vertex = boost::vertex (0, *g);
    vertex_descriptor border_vertex1 = boost::vertex (1, *g);
    vertex_descriptor border_vertex2 = boost::vertex (2, *g);
    vertex_descriptor border_vertex3 = boost::vertex (3, *g);

    edge_descriptor spoke1 = boost::add_edge (center_vertex, border_vertex1, matroid.name1 (0), *g).first;
    edge_descriptor spoke2 = boost::add_edge (center_vertex, border_vertex2, matroid.name1 (1), *g).first;
    edge_descriptor spoke3 = boost::add_edge (center_vertex, border_vertex3, matroid.name2 (2), *g).first;
    edge_descriptor rim12 = boost::add_edge (border_vertex1, border_vertex2, matroid.name2 (0), *g).first;
    edge_descriptor rim13 = boost::add_edge (border_vertex1, border_vertex3, matroid.name2 (1), *g).first;
    edge_descriptor rim23 = boost::add_edge (border_vertex2, border_vertex3, matroid.name1 (2), *g).first;

    size_t minor_height = 3;
    size_t minor_width = 3;

    for (size_t i = 0; i < nested_minors.size (); ++i)
    {
      if (!extend_graph (*g, matroid, matrix, minor_height, minor_width, nested_minors.get_extension (i)))
      {
        //        std::cout << "Test failed: not (co)graphic :-(" << std::endl;

        delete g;
        return NULL;
      }

      minor_height += nested_minors.get_extension_height (i);
      minor_width += nested_minors.get_extension_width (i);

      //        std::cout << "Currently having minor " << minor_height << " x " << minor_width << ":\n" << *g << std::endl;
    }

    return g;
  }

  inline void unittest_graphicness ()
  {
    //    std::cout << "Testing graphicness routine" << std::endl;

    matroid_graph graph;
    boost::graph_traits <matroid_graph>::vertex_descriptor V[11];

    for (int i = 0; i < 11; i++)
      V[i] = boost::add_vertex (graph);

    boost::add_edge (V[0], V[1], -1, graph);
    boost::add_edge (V[1], V[2], -2, graph);
    boost::add_edge (V[1], V[3], -3, graph);
    boost::add_edge (V[0], V[4], -4, graph);
    boost::add_edge (V[4], V[5], -5, graph);
    boost::add_edge (V[0], V[6], -6, graph);
    boost::add_edge (V[6], V[7], -7, graph);
    boost::add_edge (V[0], V[8], -8, graph);
    boost::add_edge (V[8], V[10], -9, graph);
    boost::add_edge (V[10], V[9], -10, graph);
    boost::add_edge (V[2], V[3], 1, graph);
    boost::add_edge (V[2], V[7], 2, graph);
    boost::add_edge (V[3], V[9], 3, graph);
    boost::add_edge (V[0], V[10], 4, graph);
    boost::add_edge (V[0], V[5], 5, graph);
    boost::add_edge (V[5], V[10], 6, graph);
    boost::add_edge (V[8], V[9], 7, graph);
    boost::add_edge (V[0], V[9], 8, graph);
    boost::add_edge (V[2], V[0], 9, graph);

    integer_matrix matrix (11, 9);
    for (int i = 0; i < 11 * 9; i++)
      matrix (i % 11, i / 11) = 0;
    matrix (1, 0) = 1;
    matrix (2, 0) = 1;

    matrix (0, 1) = 1;
    matrix (1, 1) = 1;
    matrix (5, 1) = 1;
    matrix (6, 1) = 1;
    matrix (10, 1) = 1;

    matrix (0, 2) = 1;
    matrix (2, 2) = 1;
    matrix (7, 2) = 1;
    matrix (8, 2) = 1;
    matrix (9, 2) = 1;
    matrix (10, 2) = 1;

    matrix (7, 3) = 1;
    matrix (8, 3) = 1;
    matrix (10, 3) = 1;

    matrix (3, 4) = 1;
    matrix (4, 4) = 1;
    matrix (10, 4) = 1;

    matrix (3, 5) = 1;
    matrix (4, 5) = 1;
    matrix (7, 5) = 1;
    matrix (8, 5) = 1;
    matrix (10, 5) = 1;

    matrix (8, 6) = 1;
    matrix (9, 6) = 1;

    matrix (7, 7) = 1;
    matrix (8, 7) = 1;
    matrix (9, 7) = 1;

    matrix (0, 8) = 1;
    matrix (1, 8) = 1;
    matrix (10, 8) = 1;

    integer_matroid matroid (11, 9);

    //    matroid_print (matroid, matrix);
    //    std::cout << graph << std::endl;

    extend_graph (graph, matroid, matrix, 10, 9, (nested_minor_sequence::extension_type) -2);
  }

}

#endif /* GRAPHICNESS_HPP_ */
