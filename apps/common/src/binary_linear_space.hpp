/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef BINARY_LINEAR_SPACE_HPP_
#define BINARY_LINEAR_SPACE_HPP_

#include <vector>

#include "permutations.hpp"

namespace tu {

  /**
   * Models a binary linear vector space by storing all vectors contained in the vector space to quickly
   * test for membership. Space consumption is O(length * 2 ^ dimension).
   */

  class binary_linear_space
  {
  public:

    typedef std::vector <bool> vector_type;

    /**
     * Contructs a binary linear vector space containing
     * vectors of the given length.
     *
     * @param length Length of each vector in the space.
     */

    explicit binary_linear_space (size_t length)
    {
      _length = length;
      _dimension = 0;
      _vectors.push_back(vector_type(length, false));
    }

    /**
     * Copy constructor which copies all vectors.
     *
     * @param other Instance of another binary linear vector space.
     */

    binary_linear_space (const binary_linear_space& other)
    {
      _length = other._length;
      _dimension = other._dimension;
      _vectors.clear();
      _vectors.reserve(other.vectors());
      std::copy(other._vectors.begin(), other._vectors.end(), std::back_inserter(_vectors));
    }

    /**
     * Read-only access to the vectors in this vector space.
     *
     * @param index Index in range [0, vectors() )
     */

    const vector_type& operator[] (size_t index) const
    {
      return _vectors[index];
    }

    /**
     * @return Length of each of the vectors
     */

    inline size_t length () const
    {
      return _length;
    }

    /**
     * @return Number of vectors in this vector space
     */

    inline size_t vectors () const
    {
      return _vectors.size();
    }

    /**
     * @return Dimension of this vector space
     */

    inline size_t dimension () const
    {
      return _dimension;
    }

    /**
     * Inserts a given vector into the vector space if it is not contained already.
     * Its size() must match the result of length().
     *
     * @param vector A given vector of appropriate size
     * @return true if and only if the vector was inserted
     */

    inline bool insert_checked (const vector_type& vector)
    {
      assert (vector.size() == length());

      if (is_spanned(vector))
        return false;

      insert(vector);
      return true;
    }

    /**
     * Inserts a given vector into the vector space without checking whether it exists.
     * Its size() must match the result of length().
     *
     * @param vector A given vector of appropriate size
     */

    void insert (const vector_type& vector)
    {
      assert (vector.size() == length());

      size_t begin = 1 << _dimension;
      ++_dimension;
      size_t end = 1 << _dimension;
      _vectors.resize(end, vector);
      for (size_t current = begin + 1; current < end; current++)
        vector_sum(vector, _vectors[current - begin], _vectors[current]);
    }

    /**
     * Searches a given vector in this vector space.
     *
     * @param vector Given vector of appropriate size
     * @return Index of the given vector or -1 if it is not spanned.
     */

    int find (const vector_type& vector) const
    {
      assert (vector.size() == length());

      for (size_t j = 0; j < _vectors.size(); ++j)
      {
        const vector_type& current_vector = _vectors[j];
        bool same = true;
        for (size_t i = 0; i < _length; ++i)
        {
          if (current_vector[i] ? !vector[i] : vector[i])
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

    /**
     * Tests whether a given vector is spanned by this vector space.
     *
     * @param vector A given vector of appropriate size
     * @return
     */

    inline bool is_spanned (const vector_type& vector) const
    {
      assert (vector.size() == length());

      return find(vector) >= 0;
    }

  protected:

    /**
     * Calculates the sum of two binary vectors and stores
     * it in a third one. All must have the same size.
     *
     * @param first First summand vector
     * @param second Second summand vector
     * @param result Vector containing the result
     */

    static void vector_sum (const vector_type& first, const vector_type& second, vector_type& result)
    {
      assert (first.size() == second.size());
      assert (first.size() == result.size());

      for (size_t i = 0; i < first.size(); ++i)
      {
        result[i] = first[i] ? !second[i] : second[i];
      }
    }

  private:
    size_t _length;
    size_t _dimension;
    std::vector <vector_type> _vectors;
  };

  /**
   * Streams a binary linear vector space.
   *
   * @param stream Output stream
   * @param space The given vector space
   * @return The output stream after processing
   */

  inline std::ostream& operator<< (std::ostream& stream, const binary_linear_space& space)
  {
    for (size_t i = 0; i < space.vectors(); ++i)
    {
      for (size_t j = 0; j < space.length(); ++j)
      {
        stream << ' ' << (space[i][j] ? '1' : '0');
      }
      stream << '\n';
    }
    return stream;
  }

}

#endif /* BINARYLINEARVECTORSPACE_H_ */
