/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef GEN_RANDOM_HPP_
#define GEN_RANDOM_HPP_

#include "gen_generic.hpp"
#include "matrix.hpp"

class random_matrix_generator: public matrix_generator
{
public:
  random_matrix_generator(size_t height, size_t width, unimod::log_level level) :
    matrix_generator("random", height, width, level)
  {

  }

  virtual ~random_matrix_generator()
  {

  }

  virtual void generate()
  {
    log_generate_start();
    boost::uniform_int <int> dist(-1, 1);

    for (size_t row = 0; row < _height; ++row)
    {
      for (size_t column = 0; column < _width; ++column)
      {
        _matrix(row, column) = dist(_rng);
      }
    }
    log_generate_end();
  }

  virtual bool do_pivot(size_t row, size_t column)
  {
    unimod::matrix_ternary_pivot(_matrix, row, column);
    return true;
  }
};

#endif
