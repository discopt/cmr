
//          Copyright Matthias Walter 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef BINARYLINEARVECTORSPACE_H_
#define BINARYLINEARVECTORSPACE_H_

#include "../config.h"
#include "permutations.hpp"
#include <vector>

namespace tu {

  class binary_linear_space
  {
  public:

    typedef std::vector <bool> vector_type;
    typedef std::vector <vector_type> data_type;

  public:
    explicit binary_linear_space (size_t size, int reserve = 1)
    {
      _size = size;
      _dimension = 0;
      _vectors.reserve (reserve);
      _vectors.push_back (vector_type (size, false));
    }

    binary_linear_space (const binary_linear_space& other)
    {
      _size = other._size;
      _dimension = other._dimension;
      _vectors.resize (other._vectors.size ());
      data_type::const_iterator other_iter = other._vectors.begin ();
      for (data_type::iterator iter = _vectors.begin (); iter != _vectors.end (); ++iter)
      {
        *iter = *other_iter;
        ++other_iter;
      }
    }

    ~binary_linear_space ()
    {

    }

    const vector_type& operator[] (size_t index) const
    {
      return _vectors[index];
    }

    inline size_t size () const
    {
      return _size;
    }

    inline size_t vectors () const
    {
      return _vectors.size ();
    }

    inline size_t dimension () const
    {
      return _dimension;
    }

    inline bool insert_checked (const std::vector <bool>& vector)
    {
      if (is_spanned (vector))
        return false;

      insert (vector);
      return true;
    }

    void insert (const vector_type& vector)
    {
      size_t begin = 1 << _dimension;
      _dimension++;
      size_t end = 1 << _dimension;
      _vectors.resize (end, vector);
      for (size_t current = begin + 1; current < end; current++)
      {
        for (size_t i = 0; i < _size; ++i)
        {
          if (_vectors[current - begin][i])
          {
            _vectors[current][i] = !_vectors[current][i];
          }
        }
      }
    }

    int get_span_index (const std::vector <bool>& vector) const
    {
      for (size_t j = 0; j < _vectors.size (); ++j)
      {
        const vector_type& current_vector = _vectors[j];
        bool same = true;
        for (size_t i = 0; i < _size; ++i)
        {
          if ((current_vector[i] && !vector[i]) || (!current_vector[i] && vector[i]))
          {
            same = false;
            break;
          }
        }
        if (same)
          return j;
      }
      return -1;
    }

    inline bool is_spanned (const std::vector <bool>& vector) const
    {
      return get_span_index (vector) >= 0;
    }

  private:
    size_t _size;
    size_t _dimension;
    data_type _vectors;
  };

  inline std::ostream& operator<< (std::ostream& stream, const binary_linear_space& space)
  {
    for (size_t i = 0; i < space.vectors (); ++i)
    {
      for (size_t j = 0; j < space.size (); ++j)
      {
        stream << ' ' << (space[i][j] ? '1' : '0');
      }
      stream << '\n';
    }
    return stream;
  }

}

#endif /* BINARYLINEARVECTORSPACE_H_ */
