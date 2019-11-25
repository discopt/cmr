#pragma once

#include <iostream>

#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace tu
{

  /**
   * Matroid element property
   */

  enum edge_matroid_element_t
  {
    edge_matroid_element
  };

  typedef boost::property <tu::edge_matroid_element_t, int> matroid_element_property;
  typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, matroid_element_property> matroid_graph;

} /* namespace tu */

namespace boost
{

  /**
   * Matroid element property tag
   */

  template <>
  struct property_kind <tu::edge_matroid_element_t>
  {
    typedef edge_property_tag type;
  };

} /* namespace boost */
 

namespace tu
{

  typedef boost::property_map <matroid_graph, edge_matroid_element_t>::const_type const_matroid_element_map;
  typedef boost::property_map <matroid_graph, edge_matroid_element_t>::type matroid_element_map;

  /**
   * Output operator for a matroid graph
   *
   * @param stream Output stream
   * @param graph The given matroid graph
   * @return The output stream after writing to it
   */

  std::ostream& operator<<(std::ostream& stream, const tu::matroid_graph& graph);

} /* namespace tu */
