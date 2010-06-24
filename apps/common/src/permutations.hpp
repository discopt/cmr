/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef PERMUTATION_HPP_
#define PERMUTATION_HPP_

#include <utility>
#include <iostream>
#include <vector>
#include <map>
#include <exception>
#include <cassert>

namespace tu {

  /**
   * Exception to indicate that the size of a permutation could not be be reduced.
   */

  class permutation_shrink_exception: std::exception
  {
  public:
    permutation_shrink_exception ()
    {

    }

    virtual ~permutation_shrink_exception () throw ()
    {

    }

    virtual const char* what () const throw ()
    {
      return "Cannot shrink permutation";
    }
  };

  /**
   * A permutation which maps integers from 0 to size-1 to the same range
   * in any possible way. This implementation stores the image of each
   * of the values in a vector.
   */

  class permutation
  {
  public:
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef size_t value_type;
    typedef std::vector <value_type> data_type;

  protected:
    data_type _data;

  public:

    /**
     * Constructs a permutation of a given size,
     * initializing to the identity.
     *
     * @param size Optional size of the permutation
     */

    permutation (size_type size = 0)
    {
      reset(size);
    }

    /**
     * Copy constructor
     *
     * @param other Another permutation
     */

    permutation (const permutation& other)
    {
      _data.resize(other.size());
      for (size_type i = 0; i < _data.size(); ++i)
        _data[i] = other._data[i];
    }

    /**
     * Destructor
     */

    virtual ~permutation ()
    {

    }

    /**
     * Resizes the permutation and resets it to identity.
     *
     * @param new_size New size
     */

    void reset (size_t new_size)
    {
      _data.resize(new_size);
      for (size_type i = 0; i < new_size; ++i)
      {
        _data[i] = i;
      }
    }

    /**
     * Resets the permutation to identity.
     */

    inline void reset ()
    {
      reset(_data.size());
    }

    /**
     * @return The current size of the permutation
     */

    inline size_type size () const
    {
      return _data.size();
    }

    /**
     * The image of a specific integer in range 0 .. size-1.
     *
     * @param index Given integer
     * @return The integers image
     */

    inline value_type operator() (value_type index) const
    {
      return get(index);
    }

    /**
     * The image of a specific integer in range 0 .. size-1.
     * @param index Given integer
     * @return The integers image
     */

    inline value_type get (value_type index) const
    {
      assert (index < _data.size());

      return _data[index];
    }

    /**
     * Swaps the images of two integers in range 0 .. size-1.
     *
     * @param a First integer
     * @param b Second integer
     */

    inline void swap (value_type a, value_type b)
    {
      std::swap(_data[a], _data[b]);
    }

    /**
     * Swaps the preimages of two integers in range 0 .. size-1.
     *
     * @param a First integer
     * @param b Second integer
     */

    void rswap (value_type a, value_type b)
    {
      value_type tmp, pa = a, pb = b;
      while ((tmp = get(pa)) != a)
        pa = tmp;
      while ((tmp = get(pb)) != b)
        pb = tmp;
      swap(pa, pb);
    }

    /**
     * Makes this permutation its own inverse.
     */

    void revert ()
    {
      /// Create a temporary copy
      value_type* temp = new value_type[size()];
      for (size_type i = 0; i < size(); ++i)
        temp[i] = _data[i];

      for (size_type i = 0; i < size(); ++i)
        _data[temp[i]] = i;
      delete[] temp;
    }

    /**
     * @return The inverse permutation of this permutation
     */

    permutation reverse () const
    {
      permutation result(size());
      for (size_type i = 0; i < size(); ++i)
        result._data[get(i)] = i;
      return result;
    }

    /**
     * Assignment operator.
     *
     * @param other Another permutation
     * @return A reference to this permutation
     */

    permutation& operator= (const permutation& other)
    {
      _data.resize(other.size());
      for (size_type i = 0; i < _data.size(); ++i)
        _data[i] = other._data[i];
      return *this;
    }

    /**
     * Calculates the product of this and another permutation.
     * The result is equivalent to applying the second permutation
     * first and this one afterwards.
     *
     * @param rhs Right hand side permutation
     * @return The resulting permutation
     */

    permutation operator* (const permutation& rhs) const
    {
      assert (size() == rhs.size());

      permutation result(size());
      for (size_type i = 0; i < size(); ++i)
        result._set(i, _data[rhs(i)]);

      return result;
    }

    /**
     * Resizes the permutation, retaining the contents.
     * This may fail if it cannot be shrunken this way and
     * will throw a permutation_shrink_exception in that case.
     *
     * @param new_size New size of the permutation
     */

    void resize (size_type new_size)
    {
      size_type old_size = size();
      for (size_type i = new_size; i < old_size; ++i)
      {
        if (_data[i] < new_size)
          throw permutation_shrink_exception();
      }

      _data.resize(new_size);
      for (size_type i = old_size; i < new_size; i++)
        _data[i] = i;
    }

    /**
     * Grows the permutation by a given number of elements.
     *
     * @param by Number of elements to increase the size by
     */

    inline void grow (difference_type by)
    {
      resize(size() + by);
    }

    /**
     * Shrinks the permutation by a given number of elements.
     *
     * @param by Number of elements to decrease the size by
     */

    inline void shrink (difference_type by)
    {
      resize(size() - by);
    }

  protected:

    /// Class to enumerate permutations

    friend class permutation_enumerator;

    /// Class to setup a permutation which can be put before a vector.

    template <class Less>
    friend void sort (permutation& permutation, size_t first, size_t beyond, Less& less);

    /**
     * Sets a specific value without checking validity.
     *
     * @param index Index of the entry
     * @param value New value
     */

    void _set (value_type index, value_type value)
    {
      _data[index] = value;
    }

    /**
     * @return A reference to the data vector
     */

    inline data_type& get_data ()
    {
      return _data;
    }

  };

  /**
   * Enumerates all permutations of a given size with
   * optional constraints. You must derive from it
   * and override the visitor() method.
   */

  class permutation_enumerator
  {
  public:
    typedef size_t size_type;
    typedef std::vector <permutation::value_type> memberlist_type;
    typedef std::pair <memberlist_type*, permutation*> groupinfo_type;
    typedef std::map <int, groupinfo_type> state_type;

  private:
    size_type _size;
    state_type _state;
    permutation _permutation;

  public:

    /**
     * Constructs an enumerator without constraints.
     * Has O(size!) running time.
     *
     * @param size Size of the permutations
     */

    permutation_enumerator (permutation::size_type size) :
      _permutation(size)
    {
      if (size)
      {
        memberlist_type* memberlist = new memberlist_type();
        memberlist->resize(size);
        for (permutation::size_type i = 0; i < size; i++)
          (*memberlist)[i] = i;
        permutation* perm = new permutation(size);
        _state[0] = groupinfo_type(memberlist, perm);
      }
    }

    /**
     * Constructs an enumerator with grouping constraints:
     * Only those permutations are enumerated where k and p(k) are
     * in the same group for all i.
     *
     * @param groups Vector of groups
     */

    permutation_enumerator (const std::vector <permutation::value_type>& groups) :
      _permutation(groups.size())
    {
      _size = groups.size();
      for (size_t i = 0; i < groups.size(); ++i)
      {
        state_type::iterator iter = _state.find(groups[i]);
        if (iter == _state.end())
        {
          _state[groups[i]] = groupinfo_type(new memberlist_type(), NULL);
          _state[groups[i]].first->push_back(i);
        }
        else
        {
          iter->second.first->push_back(i);
        }
      }
      for (state_type::iterator iter = _state.begin(); iter != _state.end(); ++iter)
      {
        groupinfo_type& info = iter->second;
        info.second = new permutation(info.first->size());
      }
    }

    /**
     * Destructor.
     */

    virtual ~permutation_enumerator ()
    {
      for (state_type::iterator iter = _state.begin(); iter != _state.end(); ++iter)
      {
        delete iter->second.first;
        delete iter->second.second;
      }
    }

    /**
     * @return true if and only if the size of each permutation is zero
     */

    inline bool empty ()
    {
      return _size == 0;
    }

    /**
     * Starts the enumeration which is aborted if any
     * visitor returns false.
     *
     * @return true if and only if all visitors returned true
     */

    virtual bool enumerate ()
    {
      if (empty())
        return true;

      return enumerate_group(_state.begin(), 0);
    }

  protected:

    /**
     * Method which is called for each enumerated permutation.
     * Needs to be overridden.
     *
     * @param perm Current permutation
     * @return false to abort enumeration
     */

    virtual bool visitor (const permutation& perm) = 0;

  private:

    /**
     * Enumerates possible values at a given index in a given group.
     *
     * @param state Current state and group
     * @param index Current index to enumerate values
     * @return false to abort enumeration
     */

    bool enumerate_group (state_type::iterator state, permutation::size_type index)
    {
      /// Nothing to enumerate - just call the visitor.
      if (state == _state.end())
      {
        for (state = _state.begin(); state != _state.end(); ++state)
        {
          const memberlist_type& memberlist = *(state->second.first);
          const permutation& perm = *(state->second.second);
          for (permutation::size_type i = 0; i < perm.size(); i++)
          {
            _permutation._set(memberlist[i], memberlist[perm(i)]);
          }
        }

        return visitor(_permutation);
      }

      /// Jump to next group
      if (index >= state->second.first->size())
      {
        return enumerate_group(++state, 0);
      }

      /// Enumerate the current groups permutation.
      permutation& p = *(state->second.second);
      for (permutation::size_type i = 0; i < p.size(); i++)
      {
        bool found = false;
        for (permutation::size_type j = 0; j < index; j++)
        {
          if (p(j) == i)
          {
            found = true;
            break;
          }
        }
        if (!found)
        {
          p._set(index, i);
          if (!enumerate_group(state, index + 1))
            return false;
        }
      }

      return true;
    }
  };

  /**
   * Compares two permutations for equality.
   *
   * @param p First permutation
   * @param q Second permutation
   * @return true if and only if the permutations are the same
   */

  inline bool operator== (const permutation& p, const permutation& q)
  {
    if (p.size() != q.size())
      return false;

    for (size_t i = 0; i < p.size(); ++i)
    {
      if (p(i) != q(i))
        return false;
    }
    return true;
  }

  /**
   * Output operator for permutations.
   *
   * @param stream A given output stream
   * @param p A given permutation
   * @return The stream after writing the permutation
   */

  inline std::ostream& operator<< (std::ostream& stream, const permutation& p)
  {
    if (p.size() > 0)
    {
      stream << p(0);
      for (size_t i = 1; i < p.size(); i++)
        stream << ' ' << p(i);
    }
    return stream;
  }

  /**
   * Arranges part of a given permutation such that it can be used to
   * view a vector in a sorted way.
   *
   * @param permutation The given permutation to be changed
   * @param first First index of the part
   * @param beyond Beyond index of the part
   * @param less Functor to compare two elements
   */

  template <class Less>
  inline void sort (permutation& permutation, size_t first, size_t beyond, Less& less)
  {
    permutation::data_type& data = permutation.get_data();
    std::sort(data.begin() + first, data.begin() + beyond, less);
  }

  /**
   * Arranges a given permutation such that it can be used to
   * view a vector in a sorted way.
   *
   * @param permutation The given permutation to be changed
   * @param less Functor to compare two elements
   */

  template <class Less>
  inline void sort (permutation& permutation, Less& less)
  {
    sort(permutation, 0, permutation.size(), less);
  }

}

#endif /* PERMUTATION_HPP_ */
