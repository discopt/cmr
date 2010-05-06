/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef NESTED_MINOR_SEQUENCE_HPP_
#define NESTED_MINOR_SEQUENCE_HPP_

#include "../config.h"
#include <vector>

namespace tu {

  class nested_minor_sequence
  {
  public:
    enum extension_type
    {
      FIRST_EXTENSION_TYPE = -2,
      ONE_ROW = -2,
      ONE_ROW_TWO_COLUMNS = -1,
      ONE_ROW_ONE_COLUMN = 0,
      TWO_ROWS_ONE_COLUMN = 1,
      ONE_COLUMN = 2,
      BEYOND_EXTENSION_TYPE = 3
    };

    nested_minor_sequence ();

    virtual ~nested_minor_sequence ();

    void push (extension_type);

    static size_t get_extension_height (extension_type ext);
    static size_t get_extension_width (extension_type ext);

    extension_type get_extension (size_t index) const
    {
      return extensions_[index];
    }

    size_t get_extension_height (size_t index) const
    {
      return get_extension_height(get_extension(index));
    }

    size_t get_extension_width (size_t index) const
    {
      return get_extension_width(get_extension(index));
    }

    size_t height () const
    {
      return height_;
    }

    size_t width () const
    {
      return width_;
    }

    size_t size () const
    {
      return extensions_.size();
    }

  private:

    std::vector <extension_type> extensions_;
    size_t height_;
    size_t width_;
  };

  class nested_minor_sequence_transposed
  {
  public:
    nested_minor_sequence_transposed (nested_minor_sequence& sequence) :
      sequence_(sequence)
    {

    }

    nested_minor_sequence_transposed (const nested_minor_sequence_transposed& other) :
      sequence_(other.sequence_)
    {

    }

    virtual ~nested_minor_sequence_transposed ()
    {

    }

    void push (nested_minor_sequence::extension_type type)
    {
      sequence_.push(nested_minor_sequence::extension_type(-int(type)));
    }

    nested_minor_sequence::extension_type get_extension (size_t index) const
    {
      return nested_minor_sequence::extension_type(-int(sequence_.get_extension(index)));
    }

    size_t get_extension_height (size_t index) const
    {
      return sequence_.get_extension_width(index);
    }

    size_t get_extension_width (size_t index) const
    {
      return sequence_.get_extension_height(index);
    }

    size_t height () const
    {
      return sequence_.width();
    }

    size_t width () const
    {
      return sequence_.height();
    }

    size_t size () const
    {
      return sequence_.size();
    }

  private:
    nested_minor_sequence& sequence_;
  };

  inline nested_minor_sequence_transposed view_nested_minor_sequence_transposed (nested_minor_sequence& sequence)
  {
    return nested_minor_sequence_transposed(sequence);
  }

}

#endif /* NESTED_MINOR_SEQUENCE_HPP_ */
