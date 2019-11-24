/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#include "../config.h"
#include "logger.hpp"

namespace unimod
{

  /**
   * Creates a logger object.
   *
   * @param level The log level to work with
   */

  logger::logger(log_level level) :
    _indent(0), _level(level)
  {
    _line = new std::stringstream();
  }

  /**
   * Destructor
   */

  logger::~logger()
  {
    delete _line;
  }

  /**
   * Streams a line of a logger object and flushes the output stream.
   *
   * @param Output stream
   * @param Logger object
   * @return Output stream
   */

  std::ostream& operator<<(std::ostream& stream, logger& log)
  {
    stream << "\r";
    for (size_t i = 0; i < log._indent; ++i)
      stream << ' ';
    stream << log.line().str() << std::flush;

    return stream;
  }
}
