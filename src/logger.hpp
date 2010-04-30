//          Copyright Matthias Walter 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef LOGGER_HPP_
#define LOGGER_HPP_

#include <sstream>
#include <iostream>

#include "total_unimodularity.hpp"

namespace tu {

  class logger
  {
  public:

  public:
    logger (log_level level);
    virtual ~logger ();

    inline log_level level () const
    {
      return _level;
    }

    inline bool is_quiet () const
    {
      return _level == LOG_QUIET;
    }

    inline bool is_verbose () const
    {
      return _level == LOG_VERBOSE;
    }

    inline bool is_updating () const
    {
      return _level == LOG_UPDATING;
    }

    inline void indent (size_t amount = 2)
    {
      _indent += amount;
    }

    inline void unindent (size_t amount = 2)
    {
      _indent -= amount;
    }

    inline void clear ()
    {
      delete _line;
      _line = new std::stringstream ();
    }

    inline size_t size () const
    {
      return line ().str ().size ();
    }

    inline void erase (size_t position)
    {
      std::string data = line ().str ();
      data.erase (position);
      clear ();
      line () << data;
    }

    inline std::stringstream& line ()
    {
      return *_line;
    }

    inline const std::stringstream& line () const
    {
      return *_line;
    }

    friend std::ostream& operator<< (std::ostream&, logger&);

  private:
    size_t _indent;
    log_level _level;
    std::stringstream* _line;

  };

  std::ostream& operator<< (std::ostream&, logger&);
}

#endif /* LOGGER_HPP_ */
