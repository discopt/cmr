#pragma once

#include "gen_generic.hpp"
#include "matrix.hpp"

class cycle_violator_matrix_generator: public matrix_generator
{
private:
  size_t _violator_size;

public:
  cycle_violator_matrix_generator(size_t size, size_t violator_size, tu::log_level level) :
    matrix_generator("cycle-violator", size, size, level)
  {
    _violator_size = violator_size;
  }

  cycle_violator_matrix_generator(size_t size, tu::log_level level) :
    matrix_generator("cycle-violator", size, size, level)
  {
    _violator_size = 0;
  }

  virtual ~cycle_violator_matrix_generator()
  {

  }

  virtual void generate()
  {
    log_generate_start();
    for (size_t row = 0; row < _height; ++row)
    {
      for (size_t column = 0; column < _width; ++column)
      {
        _matrix(row, column) = (row == column || row == (column + 1) % _height || (row < 2 && column >= 2 && column < _width - 1)) ? 1 : 0;
      }
    }
    log_generate_end();

    /// Hide violator
    if (_level != tu::LOG_QUIET)
      std::cerr << "Hiding the violator..." << std::flush;

    /// Random size
    if (_violator_size == 0)
    {
      boost::uniform_int <int> dist(3, _height);
      _violator_size = dist(_rng);
    }

    for (size_t i = 0; i < _height - _violator_size; ++i)
      tu::matrix_binary_pivot(_matrix, 2 + i, 2 + i);

    if (_level != tu::LOG_QUIET)
      std::cerr << " done. (size is " << _violator_size << " x " << _violator_size << ")" << std::endl;
  }
};
