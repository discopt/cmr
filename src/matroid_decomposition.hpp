
//          Copyright Matthias Walter 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef MATROID_DECOMPOSITION_HPP_
#define MATROID_DECOMPOSITION_HPP_

#include <set>

#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/graph/properties.hpp>

#include "matroid_graph.hpp"

namespace tu {

  typedef std::set <int> matroid_element_set;

  class decomposed_matroid
  {
  public:
    explicit decomposed_matroid (const matroid_element_set& elements, const matroid_element_set& extra_elements);
    virtual ~decomposed_matroid ();

    virtual bool is_leaf () const = 0;
    virtual bool is_graphic () const = 0;
    virtual bool is_cographic () const = 0;
    virtual bool is_regular () const = 0;

    bool is_network () const
    {
      return is_graphic () || is_cographic ();
    }

    bool is_planar () const
    {
      return is_graphic () && is_cographic ();
    }

    inline const matroid_element_set& elements () const
    {
      return _elements;
    }

    inline const matroid_element_set& extra_elements () const
    {
      return _extra_elements;
    }

  private:
    matroid_element_set _elements;
    matroid_element_set _extra_elements;
  };

  /// Class for 3-connected component

  class decomposed_matroid_leaf: public decomposed_matroid
  {
  public:

    explicit decomposed_matroid_leaf (matroid_graph* graph, matroid_graph* cograph, bool is_R10, const std::set <int>& elements,
        const matroid_element_set& extra_elements);
    virtual ~decomposed_matroid_leaf ();

    virtual bool is_leaf () const
    {
      return true;
    }

    virtual const matroid_graph* graph () const
    {
      return _graph;
    }

    virtual const matroid_graph* cograph () const
    {
      return _cograph;
    }

    virtual bool is_R10 () const
    {
      return _is_R10;
    }

    virtual bool is_graphic () const
    {
      return graph () != NULL;
    }

    virtual bool is_cographic () const
    {
      return cograph () != NULL;
    }

    virtual bool is_regular () const
    {
      return is_network () || is_R10 ();
    }

  protected:
    matroid_graph* _graph;
    matroid_graph* _cograph;
    bool _is_R10;
  };

  /// Class for separator

  class decomposed_matroid_separator: public decomposed_matroid
  {
  public:
    explicit decomposed_matroid_separator (decomposed_matroid* first, decomposed_matroid* second, int type, const std::set <int>& elements,
        const matroid_element_set& extra_elements);
    virtual ~decomposed_matroid_separator ();

    static const int ONE_SEPARATION = 1;
    static const int TWO_SEPARATION = 2;
    static const int THREE_SEPARATION = 3;

    inline int separation_type () const
    {
      return _type;
    }

    inline decomposed_matroid* first () const
    {
      return _first;
    }

    inline decomposed_matroid* second () const
    {
      return _second;
    }

    virtual bool is_leaf () const
    {
      return false;
    }

    virtual bool is_graphic () const
    {
      return _first->is_graphic () && _second->is_graphic ();
    }

    virtual bool is_cographic () const
    {
      return _first->is_cographic () && _second->is_cographic ();
    }

    virtual bool is_regular () const
    {
      return _first->is_regular () && _second->is_regular ();
    }

  protected:
    decomposed_matroid* _first;
    decomposed_matroid* _second;
    int _type;
  };

}

#endif /* MATROID_DECOMPOSITION_HPP_ */
