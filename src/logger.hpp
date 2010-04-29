//          Copyright Matthias Walter 2010.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef LOGGER_HPP_
#define LOGGER_HPP_

#include <sstream>
#include <iostream>

namespace tu {

  class logger
  {
  public:
    enum level_t
    {
      QUIET, VERBOSE
    };

  public:
    logger (level_t level);
    virtual ~logger ();

    inline level_t level () const
    {
      return _level;
    }

    inline bool is_verbose () const
    {
      return _level == VERBOSE;
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
    level_t _level;
    std::stringstream* _line;

  };

  std::ostream& operator<< (std::ostream&, logger&);
}

#endif /* LOGGER_HPP_ */
