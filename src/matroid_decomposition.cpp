/*
 * matroid_decomposition.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: xammy
 */

#include "../config.h"
#include "total_unimodularity.hpp"
#include "matroid_decomposition.hpp"
#include "matroid.hpp"

#include "algorithm.hpp"

namespace tu {

  decomposed_matroid::decomposed_matroid (const std::set <int>& elements) :
    _elements (elements)
  {

  }

  decomposed_matroid::~decomposed_matroid ()
  {

  }

  decomposed_matroid_leaf::decomposed_matroid_leaf (matroid_graph* graph, matroid_graph* cograph, bool is_R10,
      const std::set <int>& elements) :
    decomposed_matroid (elements), _graph (graph), _cograph (cograph), _is_R10 (is_R10)
  {

  }

  decomposed_matroid_leaf::~decomposed_matroid_leaf ()
  {
    if (_graph)
      delete _graph;
    if (_cograph)
      delete _cograph;
  }

  decomposed_matroid_separator::decomposed_matroid_separator (decomposed_matroid* first, decomposed_matroid* second,
      int type, const std::set <int>& elements) :
    decomposed_matroid (elements), _first (first), _second (second), _type (type)
  {

  }

  decomposed_matroid_separator::~decomposed_matroid_separator ()
  {
    delete _first;
    delete _second;
  }

  /// Returns a decomposition of a given binary matroid into a k-sum-decomposition (k=1,2,3)
  /// in graphic, cographic, R10 and maybe irregular components. 

  decomposed_matroid* decompose_binary_matroid (const boost::numeric::ublas::matrix <int>& input_matrix)
  {
    if (!is_zero_one_matrix (input_matrix))
      return NULL;

    integer_matrix matrix (input_matrix);
    integer_matroid matroid (matrix.size1 (), matrix.size2 ());

    return decompose_binary_matroid (matroid, matrix, true).second;
  }

}
