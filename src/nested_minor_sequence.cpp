/*
 * nested_minor_sequence.cpp
 *
 *  Created on: Jan 15, 2010
 *      Author: xammy
 */

#include "nested_minor_sequence.hpp"

#include "../config.h"

namespace tu {

  nested_minor_sequence::nested_minor_sequence () :
    height_ (3), width_ (3)
  {

  }

  nested_minor_sequence::~nested_minor_sequence ()
  {

  }

  size_t nested_minor_sequence::get_extension_height (extension_type ext)
  {
    static int ext_height[] = { 1, 1, 1, 2, 0 };

    return ext_height[ext - FIRST_EXTENSION_TYPE];
  }

  size_t nested_minor_sequence::get_extension_width (extension_type ext)
  {
    static int ext_width[] = { 0, 2, 1, 1, 1 };

    return ext_width[ext - FIRST_EXTENSION_TYPE];
  }

  void nested_minor_sequence::push (extension_type ext)
  {
    extensions_.push_back (ext);
    width_ += get_extension_width (ext);
    height_ += get_extension_height (ext);
  }

}
