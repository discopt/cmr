/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#include "matroid_decomposition.hpp"

namespace tu {

  /**
   * Constructs the tree node
   *
   * @param elements Set of matroid elements
   * @param extra_elements Set of extra elements, i.e. those that pivots were made upon
   */

  decomposed_matroid::decomposed_matroid (const matroid_element_set& elements, const matroid_element_set& extra_elements) :
    _elements(elements)
  {
    std::set_difference(extra_elements.begin(), extra_elements.end(), elements.begin(), elements.end(), std::inserter(_extra_elements,
        _extra_elements.end()));
  }

  /**
   * Destructor
   */

  decomposed_matroid::~decomposed_matroid ()
  {

  }

  /**
   * Constructs the node
   *
   * @param graph A corresponding graph or NULL if not graphic
   * @param cograph A corresponding cograph or NULL if not cographic
   * @param is_R10 Whether this matroid is isomorphic to R10
   * @param elements Set of elements
   * @param extra_elements Set of extra elements
   */

  decomposed_matroid_leaf::decomposed_matroid_leaf (matroid_graph* graph, matroid_graph* cograph, bool is_R10, const matroid_element_set& elements,
      const matroid_element_set& extra_elements) :
    decomposed_matroid(elements, extra_elements), _graph(graph), _cograph(cograph), _is_R10(is_R10)
  {

  }

  /**
   * Destructor
   */

  decomposed_matroid_leaf::~decomposed_matroid_leaf ()
  {
    if (_graph)
      delete _graph;
    if (_cograph)
      delete _cograph;
  }

  /**
   * Constructs the node.
   *
   * @param first First component
   * @param second Second component
   * @param type Type of the separation
   * @param elements Set of elements
   * @param extra_elements Set of extra elements
   * @return
   */

  decomposed_matroid_separator::decomposed_matroid_separator (decomposed_matroid* first, decomposed_matroid* second, int type,
      const matroid_element_set& elements, const matroid_element_set& extra_elements) :
    decomposed_matroid(elements, extra_elements), _first(first), _second(second), _type(type)
  {

  }

  /**
   * Destructor
   */

  decomposed_matroid_separator::~decomposed_matroid_separator ()
  {
    delete _first;
    delete _second;
  }

}
