//          Copyright Matthias Walter 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include "logger.hpp"

namespace tu {

  logger::logger (level_t level) :
    _indent (0), _level (level)
  {
    _line = new std::stringstream ();
  }

  logger::~logger ()
  {
    delete _line;
  }

  std::ostream& operator<< (std::ostream& stream, logger& log)
  {
    stream << "\r";
    for (size_t i = 0; i < log._indent; ++i)
      stream << ' ';
    stream << log.line ().str () << std::flush;

    return stream;
  }
}
