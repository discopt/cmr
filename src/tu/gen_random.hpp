/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef GEN_RANDOM_HPP_
#define GEN_RANDOM_HPP_

#include <boost/random/uniform_real.hpp>
#include "gen_generic.hpp"
#include "matrix.hpp"

class random_matrix_generator: public matrix_generator
{
private:
  double _nonzero_probability;

public:
  random_matrix_generator(size_t height, size_t width, double nonzero_probability, tu::log_level level) :
    matrix_generator("random", height, width, level), _nonzero_probability(nonzero_probability)
  {

  }

  virtual ~random_matrix_generator()
  {

  }

  virtual void generate()
  {
    log_generate_start();
    boost::uniform_real<double> dist;

    for (size_t row = 0; row < _height; ++row)
    {
      for (size_t column = 0; column < _width; ++column)
      {
        _matrix(row, column) = dist(_rng) > _nonzero_probability ? 0 : 1;
      }
    }
    log_generate_end();
    sign();
  }

  virtual bool do_pivot(size_t row, size_t column)
  {
    tu::matrix_ternary_pivot(_matrix, row, column);
    return true;
  }
};

#endif
