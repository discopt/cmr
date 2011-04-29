/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#include "../config.h"
#include <fstream>
#include <iomanip>

#include "total_unimodularity.hpp"
#include "matroid_decomposition.hpp"
#include "unimodularity.hpp"
#include "smith_normal_form.hpp"

void print_violator(const unimod::integer_matrix& matrix, const unimod::submatrix_indices& violator)
{
  typedef boost::numeric::ublas::matrix_indirect <const unimod::integer_matrix, unimod::submatrix_indices::indirect_array_type> indirect_matrix_t;

  const indirect_matrix_t indirect_matrix(matrix, violator.rows, violator.columns);

  for (size_t row = 0; row < indirect_matrix.size1(); ++row)
  {
    for (size_t column = 0; column < indirect_matrix.size2(); ++column)
    {
      std::cout << std::setw(4) << indirect_matrix(row, column);
    }
    std::cout << '\n';
  }
  std::cout << "\nRow indices in range [0," << (matrix.size1() - 1) << "]:\n";
  for (size_t row = 0; row < violator.rows.size(); ++row)
    std::cout << (row == 0 ? "" : " ") << violator.rows[row];
  std::cout << "\n\nColumn indices in range [0," << (matrix.size2() - 1) << "]:\n";
  for (size_t column = 0; column < violator.columns.size(); ++column)
    std::cout << (column == 0 ? "" : " ") << violator.columns[column] << ' ';
  std::cout << std::endl;
}

int run_decomposition(const std::string& file_name, bool show_certificates, unimod::log_level level)
{
  /// Open the file

  std::ifstream file(file_name.c_str());
  if (!file.good())
  {
    std::cout << "Error: cannot open file \"" << file_name << "\"." << std::endl;
    return EXIT_FAILURE;
  }

  /// Read height and width

  size_t height, width;
  file >> height >> width;
  if (!file.good())
  {
    std::cout << "Error: cannot read height and width from input file." << std::endl;
    return EXIT_FAILURE;
  }

  /// Read matrix entries

  unimod::integer_matrix matrix(height, width);
  for (size_t row = 0; row < height; ++row)
  {
    for (size_t column = 0; column < width; ++column)
    {
      if (!file.good())
      {
        std::cout << "Error: cannot read matrix data." << std::endl;
      }
      int value;
      file >> value;
      matrix(row, column) = value;
    }
  }

  file.close();

  if (show_certificates)
  {
    unimod::submatrix_indices violator;
    unimod::decomposed_matroid* decomposition;

    if (unimod::is_totally_unimodular(matrix, decomposition, violator, level))
    {
      std::cout << "\nThe " << matrix.size1() << " x " << matrix.size2() << " matrix is totally unimodular.\n" << std::endl;
    }
    else
    {
      std::cout << "\nThe " << matrix.size1() << " x " << matrix.size2() << " matrix is not totally unimodular." << std::endl;
      assert (violator.rows.size() == violator.columns.size());
      int det = unimod::submatrix_determinant(matrix, violator);
      std::cout << "\nThe violating submatrix (det = " << det << ") is " << violator.rows.size() << " x " << violator.columns.size() << ":\n\n"
          << std::flush;
      print_violator(matrix, violator);
    }
    delete decomposition;
  }
  else
  {
    if (unimod::is_totally_unimodular(matrix, level))
    {
      std::cout << "The " << matrix.size1() << " x " << matrix.size2() << " matrix is totally unimodular." << std::endl;
    }
    else
    {
      std::cout << "The " << matrix.size1() << " x " << matrix.size2() << " matrix is not totally unimodular." << std::endl;
    }
  }

  return EXIT_SUCCESS;
}

int run_column_enumeration(const std::string& file_name)
{
  /// Open the file

  std::ifstream file(file_name.c_str());
  if (!file.good())
  {
    std::cout << "Error: cannot open file \"" << file_name << "\"." << std::endl;
    return EXIT_FAILURE;
  }

  /// Read height and width

  size_t height, width;
  file >> height >> width;
  if (!file.good())
  {
    std::cout << "Error: cannot read height and width from input file." << std::endl;
    return EXIT_FAILURE;
  }

  /// Read matrix entries

  unimod::integer_matrix matrix(height, width);
  for (size_t row = 0; row < height; ++row)
  {
    for (size_t column = 0; column < width; ++column)
    {
      if (!file.good())
      {
        std::cout << "Error: cannot read matrix data." << std::endl;
      }
      int value;
      file >> value;
      matrix(row, column) = value;
    }
  }

  file.close();

  if (unimod::ghouila_houri_is_totally_unimodular(matrix))
  {
    std::cout << "The " << matrix.size1() << " x " << matrix.size2() << " matrix is totally unimodular." << std::endl;
  }
  else
  {
    std::cout << "The " << matrix.size1() << " x " << matrix.size2() << " matrix is not totally unimodular." << std::endl;
  }

  return EXIT_SUCCESS;
}

int run_submatrix(const std::string& file_name)
{
  /// Open the file

  std::ifstream file(file_name.c_str());
  if (!file.good())
  {
    std::cout << "Error: cannot open file \"" << file_name << "\"." << std::endl;
    return EXIT_FAILURE;
  }

  /// Read height and width

  size_t height, width;
  file >> height >> width;
  if (!file.good())
  {
    std::cout << "Error: cannot read height and width from input file." << std::endl;
    return EXIT_FAILURE;
  }

  /// Read matrix entries

  unimod::integer_matrix matrix(height, width);
  for (size_t row = 0; row < height; ++row)
  {
    for (size_t column = 0; column < width; ++column)
    {
      if (!file.good())
      {
        std::cout << "Error: cannot read matrix data." << std::endl;
      }
      int value;
      file >> value;
      matrix(row, column) = value;
    }
  }

  file.close();

  unimod::submatrix_indices violator_indices;
  if (unimod::determinant_is_totally_unimodular(matrix, violator_indices))
  {
    std::cout << "The " << matrix.size1() << " x " << matrix.size2() << " matrix is totally unimodular." << std::endl;
  }
  else
  {
    std::cout << "The " << matrix.size1() << " x " << matrix.size2() << " matrix is not totally unimodular with violator of size "
        << violator_indices.rows.size() << "." << std::endl;
  }

  return EXIT_SUCCESS;
}

bool extract_option(char c, char& algorithm, bool& certs, unimod::log_level& level, bool& help)
{
  if (c == 'D' || c == 'C' || c == 'S')
    algorithm = c;
  else if (c == 'h')
    help = true;
  else if (c == 'c')
    certs = true;
  else if (c == 'q')
    level = unimod::LOG_QUIET;
  else if (c == 'p')
    level = unimod::LOG_PROGRESSIVE;
  else if (c == 'v')
    level = unimod::LOG_VERBOSE;
  else
    return false;

  return true;
}

int main(int argc, char **argv)
{
  /// Possible parameters
  std::string matrix_file_name = "";
  bool certs = false;
  unimod::log_level level = unimod::LOG_PROGRESSIVE;
  bool help = false;
  char algorithm = 'm';

  bool options_done = false;
  for (int a = 1; a < argc; ++a)
  {
    const std::string current = argv[a];

    if (!options_done)
    {
      if (current == std::string("--"))
      {
        options_done = true;
        continue;
      }
      else if (current != "" && current[0] == '-')
      {
        for (size_t i = 1; i < current.size(); ++i)
        {
          if (!extract_option(current[i], algorithm, certs, level, help))
          {
            std::cerr << "Unknown option: -" << current[i] << "\nSee " << argv[0] << " -h for usage." << std::endl;
            return EXIT_FAILURE;
          }
        }
        continue;
      }
    }

    if (matrix_file_name != "")
    {
      std::cerr << "Matrix file was given twice!\nSee " << argv[0] << " -h for usage." << std::endl;
      return EXIT_FAILURE;
    }
    matrix_file_name = current;
  }

  if (help)
  {
    std::cerr << "Usage: " << argv[0] << " [OPTIONS] [--] MATRIX_FILE\n";
    std::cerr << "Options:\n";
    std::cerr << " -h Shows a help message.\n";
    std::cerr << " -D Test total unimodularity via matroid decomposition algorithm (default).\n";
    std::cerr << " -C Test total unimodularity via column enumeration algorithm (slow).\n";
    std::cerr << " -S Test total unimodularity via submatrix enumeration algorithm (very slow!).\n";
    std::cerr << " -c Calculates a violating submatrix if one exists. (only decomposition algorithm)\n";
    std::cerr << " -p Progressive logging (default, affects only decomposition algorithm).\n";
    std::cerr << " -v Verbose logging. (only decomposition algorithm)\n";
    std::cerr << " -q No logging at all. (only decomposition algorithm)\n";
    std::cerr << std::flush;
    return EXIT_SUCCESS;
  }

  if (matrix_file_name == "")
  {
    std::cerr << "No matrix file was given!\nSee " << argv[0] << " -h for usage." << std::endl;
    return EXIT_FAILURE;
  }

  if (algorithm != 'D' && certs)
  {
    std::cout << "Certificates are only available for decomposition algorithm!\nSee " << argv[0] << " -h for usage." << std::endl;
    return EXIT_FAILURE;
  }
  if (algorithm != 'D' && level != unimod::LOG_PROGRESSIVE)
  {
    std::cout << "Logging options only have an affect on decomposition algorithm!" << std::endl;
  }

  if (algorithm == 'D')
    return run_decomposition(matrix_file_name, certs, level);
  else if (algorithm == 'C')
    return run_column_enumeration(matrix_file_name);
  else if (algorithm == 'S')
    return run_submatrix(matrix_file_name);
  else
  {
    std::cerr << "Fatal error: Invalid algorithm selected." << std::endl;
    return EXIT_FAILURE;
  }
}
