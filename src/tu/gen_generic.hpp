#pragma once

#include <boost/random.hpp>
#include <tu/total_unimodularity.hpp>
#include "matrix_transposed.hpp"
#include "matrix_reorder.hpp"
#include "permutations.hpp"
#include "matrix_transposed_permuted.hpp"

class matrix_generator
{
protected:
  size_t _height;
  size_t _width;
  tu::integer_matrix _matrix;
  boost::mt19937 _rng;
  tu::log_level _level;
  const char* _name;

public:
  matrix_generator(const char* name, size_t height, size_t width, tu::log_level level) :
    _height(height), _width(width), _matrix(height, width), _rng(time(NULL)), _level(level), _name(name)
  {
  }

  virtual ~matrix_generator()
  {

  }

  virtual void generate() = 0;

  void log_generate_start()
  {
    if (_level != tu::LOG_QUIET)
      std::cerr << "Generating " << _height << " x " << _width << " " << _name << " matrix..." << std::flush;
  }

  void log_generate_end()
  {
    if (_level != tu::LOG_QUIET)
      std::cerr << " done." << std::endl;
  }

  virtual bool do_pivot(size_t row, size_t column)
  {
    return true;
  }

  template <typename MatrixType>
  void permute_rows(MatrixType& matrix)
  {
    tu::permutation perm(matrix.size1(), _rng);
    tu::matrix_apply_row_permutation(matrix, perm);
  }

  void permute_rows()
  {
    permute_rows(_matrix);
  }

  void permute_columns()
  {
    tu::matrix_transposed <tu::integer_matrix> transposed(_matrix);
    permute_rows(transposed);
  }

  virtual void randomize()
  {
    permute_rows();
    permute_columns();
  }

  virtual void sign()
  {
    tu::sign_matrix(_matrix);
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
