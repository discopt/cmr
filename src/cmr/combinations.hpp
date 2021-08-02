#pragma once

#include <vector>

namespace tu
{

  /**
   * Exception to indicate that we are at the last combination already.
   */

  class last_combination_exception: std::exception
  {
  public:
    last_combination_exception()
    {

    }

    virtual ~last_combination_exception() throw ()
    {

    }

    virtual const char* what() const throw ()
    {
      return "Already last combination";
    }
  };

  class combination
  {
  protected:
    std::vector <size_t> _data;
    size_t _size;

  public:
    /// Init with (1, ..., 1, 0, ..., 0)

    combination(size_t n, size_t k) :
      _data(k), _size(n)
    {
      assert(k <= n);
      for (size_t i = 0; i < k; ++i)
      {
        _data[i] = i;
      }
    }

    /// Where is the index's 1 entry?

    size_t operator[](size_t index) const
    {
      return _data[index];
    }

    size_t n() const
    {
      return _size;
    }

    size_t k() const
    {
      return _data.size();
    }

    /**
     * Is this the last, i.e. (0, ..., 0, 1, ..., 1) ?
     * Time: O(1)
     */

    bool is_last() const
    {
      return _data.empty() || _data[0] == _size - _data.size();
    }

    /// Compute the next combination.
    /// Time: O(k)

    void next()
    {
      if (is_last())
        throw last_combination_exception();

      // Find end of first consecutive block.
      size_t last_consecutive_one = 0;
      while (last_consecutive_one + 1 < _data.size())
      {
        if (_data[last_consecutive_one + 1] != _data[last_consecutive_one] + 1)
          break;
        last_consecutive_one++;
      }

      // Move end one entry further.
      _data[last_consecutive_one]++;

      // Move all others of block to beginning.
      for (size_t i = 0; i < last_consecutive_one; ++i)
      {
        _data[i] = i;
      }
    }
  };

  /**
   * Output operator for combinations.
   *
   * @param stream A given output stream
   * @param p A given permutation
   * @return The stream after writing the combination
   */

  inline std::ostream& operator<<(std::ostream& stream, const combination& c)
  {
    size_t next = (c.k() > 0 && c[0] == 0) ? 1 : 0;
    stream << next;

    for (size_t i = 1; i < c.n(); ++i)
    {
      if ((next < c.k()) && c[next] == i)
      {
        stream << " 1";
        next++;
      }
      else
      {
        stream << " 0";
      }
    }
    return stream;
  }

} /* namespace tu */
