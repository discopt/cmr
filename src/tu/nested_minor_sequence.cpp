/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#include "nested_minor_sequence.hpp"

namespace tu
{

  /**
   * Constructs a sequence which only consists of a W3-minor.
   */

  nested_minor_sequence::nested_minor_sequence() :
    height_(3), width_(3)
  {

  }

  /**
   * Destructor
   */

  nested_minor_sequence::~nested_minor_sequence()
  {

  }

  /**
   * Returns the number of rows in an extension of a given type.
   *
   * @param ext Extension type
   * @return Number of rows
   */

  size_t nested_minor_sequence::get_extension_height(extension_type type)
  {
    static int ext_height[] =
    { 1, 1, 1, 2, 0 };

    return ext_height[type - FIRST_EXTENSION_TYPE];
  }

  /**
   * Returns the number of columns in an extension of a given type.
   *
   * @param ext Extension type
   * @return Number of columns
   */

  size_t nested_minor_sequence::get_extension_width(extension_type type)
  {
    static int ext_width[] =
    { 0, 2, 1, 1, 1 };

    return ext_width[type - FIRST_EXTENSION_TYPE];
  }

  /**
   * Augments the sequence with a minor of the specified extension type.
   *
   * @param type The extension type of the biggest minor
   */

  void nested_minor_sequence::push(extension_type type)
  {
    extensions_.push_back(type);
    width_ += get_extension_width(type);
    height_ += get_extension_height(type);
  }

}
