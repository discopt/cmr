/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef GEN_GENERIC_HPP_
#define GEN_GENERIC_HPP_

#include <boost/random.hpp>
#include "total_unimodularity.hpp"
#include "matrix_transposed.hpp"
#include "matrix_reorder.hpp"
#include "permutations.hpp"

class matrix_generator
{
protected:
  size_t _height;
  size_t _width;
  unimod::integer_matrix _matrix;
  boost::mt19937 _rng;
  unimod::log_level _level;
  const char* _name;

public:
  matrix_generator(const char* name, size_t height, size_t width, unimod::log_level level) :
    _height(height), _width(width), _matrix(height, width), _rng(time(NULL)), _level(level), _name(name)
  {
  }

  virtual ~matrix_generator()
  {

  }

  virtual void generate() = 0;

  void log_generate_start()
  {
    if (_level != unimod::LOG_QUIET)
      std::cerr << "Generating " << _height << " x " << _width << " " << _name << " matrix..." << std::flush;
  }

  void log_generate_end()
  {
    if (_level != unimod::LOG_QUIET)
      std::cerr << " done." << std::endl;
  }

  virtual bool do_pivot(size_t row, size_t column)
  {
    return true;
  }

  template <typename MatrixType>
  void permute_rows(MatrixType& matrix)
  {
    unimod::permutation perm(matrix.size1(), _rng);
    unimod::matrix_apply_row_permutation(matrix, perm);
  }

  void permute_rows()
  {
    permute_rows(_matrix);
  }

  void permute_columns()
  {
    unimod::matrix_transposed <unimod::integer_matrix> transposed(_matrix);
    permute_rows(transposed);
  }

  virtual void randomize()
  {
    permute_rows();
    permute_columns();
  }

  virtual void sign()
  {
    unimod::sign_matrix(_matrix);
  }

  virtual void print()
  {
    std::cout << _height << " " << _width << "\n\n";
    for (size_t row = 0; row < _height; ++row)
    {
      for (size_t column = 0; column < _width; ++column)
      {
        std::cout << std::setw(3) << _matrix(row, column);
      }
      std::cout << "\n";
    }
    std::cout << std::flush;
  }
};

#endif
