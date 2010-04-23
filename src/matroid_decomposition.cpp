/*
 * matroid_decomposition.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: xammy
 */

#include "../config.h"
#include "matroid_decomposition.hpp"

namespace tu {

  decomposed_matroid::decomposed_matroid (const matroid_element_set& elements, const matroid_element_set& extra_elements) :
    _elements (elements)
  {
    std::set_difference (extra_elements.begin (), extra_elements.end (), elements.begin (), elements.end (), std::inserter (_extra_elements,
        _extra_elements.end ()));
  }

  decomposed_matroid::~decomposed_matroid ()
  {

  }

  decomposed_matroid_leaf::decomposed_matroid_leaf (matroid_graph* graph, matroid_graph* cograph, bool is_R10, const matroid_element_set& elements,
      const matroid_element_set& extra_elements) :
    decomposed_matroid (elements, extra_elements), _graph (graph), _cograph (cograph), _is_R10 (is_R10)
  {

  }

  decomposed_matroid_leaf::~decomposed_matroid_leaf ()
  {
    if (_graph)
      delete _graph;
    if (_cograph)
      delete _cograph;
  }

  decomposed_matroid_separator::decomposed_matroid_separator (decomposed_matroid* first, decomposed_matroid* second, int type,
      const matroid_element_set& elements, const matroid_element_set& extra_elements) :
    decomposed_matroid (elements, extra_elements), _first (first), _second (second), _type (type)
  {

  }

  decomposed_matroid_separator::~decomposed_matroid_separator ()
  {
    delete _first;
    delete _second;
  }

}
