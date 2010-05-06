/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#include "matroid_graph.hpp"
#include <map>
#include <string>
#include <iomanip>

namespace tu {

  std::ostream& operator<< (std::ostream& stream, const tu::matroid_graph& graph)
  {
    typedef boost::graph_traits <tu::matroid_graph> traits;

    std::map <size_t, std::string> names;

    traits::vertex_iterator vertex_iter, vertex_end;
    std::string base = "";
    size_t index = boost::num_vertices(graph);
    if (index > 0)
      index--;
    while (index > 0)
    {
      index /= 26;
      base += ' ';
    }

    index = 0;
    for (boost::tie(vertex_iter, vertex_end) = boost::vertices(graph); vertex_iter != vertex_end; ++vertex_iter, ++index)
    {
      std::string name = base;
      size_t j = index;
      for (size_t i = 0; i < base.length(); i++)
      {
        name[name.length() - 1 - i] = (char) ('A' + j % 26);
        j /= 26;
      }
      names[(size_t) (*vertex_iter)] = name;

      stream << "Vertex " << name << " has descriptor " << std::hex << (int) *vertex_iter << std::dec << '\n';
    }
    stream << '\n';

    tu::const_matroid_element_map element_map = boost::get(tu::edge_matroid_element, graph);
    traits::edge_iterator edge_iter, edge_end;
    for (boost::tie(edge_iter, edge_end) = boost::edges(graph); edge_iter != edge_end; ++edge_iter)
    {
      traits::vertex_descriptor s = boost::source(*edge_iter, graph);
      traits::vertex_descriptor t = boost::target(*edge_iter, graph);

      stream << names[(size_t) (s)] << " - " << names[(size_t) (t)] << " via";
      stream << std::setw(6) << element_map[*edge_iter] << '\n';
    }

    return stream;
  }

}
