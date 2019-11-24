/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef LOGGER_HPP_
#define LOGGER_HPP_

#include <sstream>
#include <iostream>

#include <tu/common.hpp>

namespace unimod
{

  /**
   * This class manages the logging behaviour of the library. There are three level:
   * - quiet: prints nothing
   * - verbose: prints a line for each important step
   * - updating: updates the current line incrementally
   */

  class logger
  {
  public:

    /**
     * Creates a logger object.
     *
     * @param level The log level to work with
     */

    logger(log_level level);

    /**
     * Destructor
     */

    virtual ~logger();

    /**
     * @return The current log level
     */

    inline log_level level() const
    {
      return _level;
    }

    /**
     * @return true if and only if log level is quiet
     */

    inline bool is_quiet() const
    {
      return _level == LOG_QUIET;
    }

    /**
     * @return true if and only if log level is verbose
     */

    inline bool is_verbose() const
    {
      return _level == LOG_VERBOSE;
    }

    /**
     * @return true if and only if log level is updating
     */

    inline bool is_progressive() const
    {
      return _level == LOG_PROGRESSIVE;
    }

    /**
     * Increases the indent of current and further lines.
     *
     * @param amount Number of spaces to increase
     */

    inline void indent(size_t amount = 1)
    {
      _indent += amount;
    }

    /**
     * Decreases the indent of current and further lines.
     *
     * @param amount Number of space to decrease
     */

    inline void unindent(size_t amount = 1)
    {
      _indent -= amount;
    }

    /**
     * Clears the current line.
     */

    inline void clear()
    {
      delete _line;
      _line = new std::stringstream();
    }

    /**
     * @return Number of characters in the current line
     */

    inline size_t size() const
    {
      return line().str().size();
    }

    /**
     * Erases a portion from the current line.
     *
     * @param position Position from which to erase
     */

    inline void erase(size_t position)
    {
      std::string data = line().str();
      data.erase(position);
      clear();
      line() << data;
    }

    /**
     * @return The string stream holding the current line
     */

    inline std::stringstream& line()
    {
      return *_line;
    }

    /**
     * @return The string stream holding the current line
     */

    inline const std::stringstream& line() const
    {
      return *_line;
    }

    /// Stream operator

    friend std::ostream& operator<<(std::ostream&, logger&);

  private:
    size_t _indent;
    log_level _level;
    std::stringstream* _line;

  };

  /**
   * Streams a line of a logger object and flushes the output stream.
   *
   * @param Output stream
   * @param Logger object
   * @return Output stream
   */

  std::ostream& operator<<(std::ostream&, logger&);
}

#endif /* LOGGER_HPP_ */
