/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef MATROID_DECOMPOSITION_HPP_
#define MATROID_DECOMPOSITION_HPP_

#include <set>

#include <boost/graph/adjacency_list.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/graph/properties.hpp>

#include "matroid_graph.hpp"

namespace unimod
{

  /**
   * A set of matroid elements. Here, the original row elements are numbered -1, ..., -h
   * and the original column elements +1, ..., +w
   */

  typedef std::set <int> matroid_element_set;

  /**
   * Abstract node in a decomposition tree
   */

  class decomposed_matroid
  {
  public:
    /**
     * Constructs the tree node
     *
     * @param elements Set of matroid elements
     * @param extra_elements Set of extra elements, i.e. those that pivots were made upon
     */

    explicit decomposed_matroid(const matroid_element_set& elements, const matroid_element_set& extra_elements);

    /**
     * Destructor
     */

    virtual ~decomposed_matroid();

    /**
     * @return true if and only if this node is a leaf
     */

    virtual bool is_leaf() const = 0;

    /**
     * @return true if and only if the matroid is graphic
     */

    virtual bool is_regular() const = 0;

    /**
     * @return A read-only reference to the elements
     */

    inline const matroid_element_set& elements() const
    {
      return _elements;
    }

    /**
     * @return A read-only reference to the extra elements
     */

    inline const matroid_element_set& extra_elements() const
    {
      return _extra_elements;
    }

  private:
    matroid_element_set _elements;
    matroid_element_set _extra_elements;
  };

  /**
   * Leaf node in a decomposition tree
   */

  class decomposed_matroid_leaf: public decomposed_matroid
  {
  public:

    /**
     * Constructs the node
     *
     * @param graph A corresponding graph or NULL if not graphic
     * @param cograph A corresponding cograph or NULL if not cographic
     * @param is_R10 Whether this matroid is isomorphic to R10
     * @param elements Set of elements
     * @param extra_elements Set of extra elements
     */

    explicit decomposed_matroid_leaf(matroid_graph* graph, matroid_graph* cograph, bool is_R10, const std::set <int>& elements,
        const matroid_element_set& extra_elements);

    /**
     * Destructor
     */

    virtual ~decomposed_matroid_leaf();

    /**
     * @return true
     */

    virtual bool is_leaf() const
    {
      return true;
    }

    /**
     * @return A corresponding graph or NULL if the matroid is not graphic
     */

    const matroid_graph* graph() const
    {
      return _graph;
    }

    /**
     * @returnA corresponding cograph or NULL if the matroid is not cographic
     */

    const matroid_graph* cograph() const
    {
      return _cograph;
    }

    /**
     * @return true if and only if the matroid is isomorphic to R10
     */

    bool is_R10() const
    {
      return _is_R10;
    }

    /**
     * @return true if and only if the matroid is graphic
     */

    bool is_graphic() const
    {
      return graph() != NULL;
    }

    /**
     * @return true if and only if the matroid is cographic
     */

    bool is_cographic() const
    {
      return cograph() != NULL;
    }

    /**
     * @return true if and only if the matroid is regular
     */

    virtual bool is_regular() const
    {
      return is_graphic() || is_cographic() || is_R10();
    }

  protected:
    matroid_graph* _graph;
    matroid_graph* _cograph;
    bool _is_R10;
  };

  /**
   * Separator node in a decomposition tree
   */

  class decomposed_matroid_separator: public decomposed_matroid
  {
  public:

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

    explicit decomposed_matroid_separator(decomposed_matroid* first, decomposed_matroid* second, int type, const std::set <int>& elements,
        const matroid_element_set& extra_elements);

    /**
     * Destructor
     */

    virtual ~decomposed_matroid_separator();

    enum { ONE_SEPARATION = 1 };
    enum { TWO_SEPARATION = 2 };
    enum { THREE_SEPARATION = 3 };

    /**
     * @return Type of the separation
     */

    inline int separation_type() const
    {
      return _type;
    }

    /**
     * @return First component
     */

    inline decomposed_matroid* first() const
    {
      return _first;
    }

    /**
     * @return Second component
     */

    inline decomposed_matroid* second() const
    {
      return _second;
    }

    /**
     * @return false
     */

    virtual bool is_leaf() const
    {
      return false;
    }

    /**
     * @return true if and only both components are regular
     */

    virtual bool is_regular() const
    {
      return _first->is_regular() && _second->is_regular();
    }

  protected:
    decomposed_matroid* _first;
    decomposed_matroid* _second;
    int _type;
  };

}

#endif /* MATROID_DECOMPOSITION_HPP_ */
