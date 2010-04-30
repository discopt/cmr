
//          Copyright Matthias Walter 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef PERMUTATION_HPP_
#define PERMUTATION_HPP_

#include "../config.h"
#include <utility>
#include <iostream>
#include <vector>
#include <map>
#include <exception>
#include <cassert>

namespace tu {

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

  class permutation
  {
  public:
    typedef size_t size_type;
    typedef ptrdiff_t difference_type;
    typedef size_t value_type;
    typedef std::vector <value_type> data_type;

  protected:
    data_type data_;

  public:
    permutation (size_type size = 0)
    {
      reset (size);
    }

    permutation (const permutation& other)
    {
      data_.resize (other.size ());
      for (size_type i = 0; i < data_.size (); ++i)
        data_[i] = other.data_[i];
    }

    virtual ~permutation ()
    {

    }

    void reset (size_t new_size)
    {
      data_.resize (new_size);
      for (size_type i = 0; i < new_size; ++i)
      {
        data_[i] = i;
      }
    }

    inline void reset ()
    {
      reset (data_.size ());
    }

    inline size_type size () const
    {
      return data_.size ();
    }

    // Operations

    inline value_type operator() (value_type index) const
    {
      assert (index < data_.size());
      return get (index);
    }

    inline value_type get (value_type index) const
    {
      return data_[index];
    }

    inline void swap (value_type a, value_type b)
    {
      std::swap (data_[a], data_[b]);
    }

    void rswap (value_type a, value_type b)
    {
      value_type tmp, pa = a, pb = b;
      while ((tmp = get (pa)) != a)
        pa = tmp;
      while ((tmp = get (pb)) != b)
        pb = tmp;
      swap (pa, pb);
    }

    void revert ()
    {
      // temporary copy
      value_type* temp = new value_type[size ()];
      for (size_type i = 0; i < size (); ++i)
        temp[i] = data_[i];

      for (size_type i = 0; i < size (); ++i)
        data_[temp[i]] = i;
      delete[] temp;
    }

    permutation reverse () const
    {
      permutation result (size ());
      for (size_type i = 0; i < size (); ++i)
        result.data_[get (i)] = i;
      return result;
    }

    permutation& operator= (const permutation& other)
    {
      data_.resize (other.size ());
      for (size_type i = 0; i < data_.size (); ++i)
        data_[i] = other.data_[i];
      return *this;
    }

    permutation operator* (const permutation& rhs) const
    {
      assert (size() == rhs.size());

      permutation result (size ());
      for (size_type i = 0; i < size (); ++i)
        result._set (i, data_[rhs (i)]);

      return result;
    }

    // Resizing stuff

    void resize (size_type new_size)
    {
      size_type old_size = size ();
      for (size_type i = new_size; i < old_size; ++i)
      {
        if (data_[i] < new_size)
          throw permutation_shrink_exception ();
      }

      data_.resize (new_size);
      for (size_type i = old_size; i < new_size; i++)
        data_[i] = i;
    }

    inline void grow (difference_type by)
    {
      resize (size () + by);
    }

    inline void shrink (difference_type by)
    {
      resize (size () - by);
    }

  protected:

    friend class permutation_enumerator;

    template <class Less>
    friend void sort (permutation& permutation, size_t first, size_t beyond, Less& less);

    void _set (value_type index, value_type value)
    {
      data_[index] = value;
    }

    inline data_type& get_data ()
    {
      return data_;
    }

  };

  class permutation_enumerator
  {
  public:
    typedef size_t size_type;
    typedef std::vector <permutation::value_type> groups_type;
    typedef std::vector <permutation::value_type> memberlist_type;
    typedef std::pair <memberlist_type*, permutation*> groupinfo_type;
    typedef std::map <int, groupinfo_type> state_type;

  private:
    size_type _size;
    state_type _state;
    permutation _permutation;

  public:
    permutation_enumerator (permutation::size_type size) :
      _permutation (size)
    {
      if (size)
      {
        memberlist_type* memberlist = new memberlist_type ();
        memberlist->resize (size);
        for (permutation::size_type i = 0; i < size; i++)
          (*memberlist)[i] = i;
        permutation* perm = new permutation (size);
        _state[0] = groupinfo_type (memberlist, perm);
      }
    }

    permutation_enumerator (const groups_type& groups) :
      _permutation (groups.size ())
    {
      _size = groups.size ();
      for (size_t i = 0; i < groups.size (); ++i)
      {
        state_type::iterator iter = _state.find (groups[i]);
        if (iter == _state.end ())
        {
          _state[groups[i]] = groupinfo_type (new memberlist_type (), NULL);
          _state[groups[i]].first->push_back (i);
        }
        else
        {
          iter->second.first->push_back (i);
        }
      }
      for (state_type::iterator iter = _state.begin (); iter != _state.end (); ++iter)
      {
        groupinfo_type& info = iter->second;
        info.second = new permutation (info.first->size ());
      }
    }

    virtual ~permutation_enumerator ()
    {
      for (state_type::iterator iter = _state.begin (); iter != _state.end (); ++iter)
      {
        delete iter->second.first;
        delete iter->second.second;
      }
    }

    inline bool empty ()
    {
      return _size == 0;
    }

    virtual bool enumerate ()
    {
      if (empty ())
        return true;

      return enumerate_group (_state.begin (), 0);
    }

  protected:

    virtual bool visitor (const permutation& perm) = 0;

  private:

    bool enumerate_group (state_type::iterator state, permutation::size_type index)
    {
      if (state == _state.end ())
      {
        for (state = _state.begin (); state != _state.end (); ++state)
        {
          const memberlist_type& memberlist = *(state->second.first);
          const permutation& perm = *(state->second.second);
          for (permutation::size_type i = 0; i < perm.size (); i++)
          {
            _permutation._set (memberlist[i], memberlist[perm (i)]);
          }
        }

        return visitor (_permutation);
      }

      if (index >= state->second.first->size ())
      {
        return enumerate_group (++state, 0);
      }

      permutation& p = *(state->second.second);

      for (permutation::size_type i = 0; i < p.size (); i++)
      {
        bool found = false;
        for (permutation::size_type j = 0; j < index; j++)
        {
          if (p (j) == i)
          {
            found = true;
            break;
          }
        }
        if (!found)
        {
          p._set (index, i);
          if (!enumerate_group (state, index + 1))
            return false;
        }
      }

      return true;
    }
  };

  inline bool operator== (const permutation& p, const permutation& q)
  {
    if (p.size () != q.size ())
      return false;

    for (size_t i = 0; i < p.size (); ++i)
    {
      if (p (i) != q (i))
        return false;
    }
    return true;
  }

  inline std::ostream& operator<< (std::ostream& stream, const permutation& p)
  {
    if (p.size () > 0)
    {
      stream << p (0);
      for (size_t i = 1; i < p.size (); i++)
        stream << ' ' << p (i);
    }
    return stream;
  }

  template <class Less>
  inline void sort (permutation& permutation, size_t first, size_t beyond, Less& less)
  {
    permutation::data_type& data = permutation.get_data ();
    std::sort (data.begin () + first, data.begin () + beyond, less);
  }

  template <class Less>
  inline void sort (permutation& permutation, Less& less)
  {
    sort (permutation, 0, permutation.size (), less);
  }

}

#endif /* PERMUTATION_HPP_ */
