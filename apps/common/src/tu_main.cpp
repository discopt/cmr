/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#include <fstream>
#include <iomanip>

#include "total_unimodularity.hpp"
#include "matroid_decomposition.hpp"
#include "config.h"

void print_matroid_graph (const tu::matroid_graph& graph, const std::string& indent = "")
{
  std::cout << boost::num_vertices(graph) << " nodes and " << boost::num_edges(graph) << " edges:";

  typedef boost::graph_traits <tu::matroid_graph> traits;
  traits::vertex_iterator vertex_iter, vertex_end;
  traits::out_edge_iterator edge_iter, edge_end;

  for (boost::tie(vertex_iter, vertex_end) = boost::vertices(graph); vertex_iter != vertex_end; ++vertex_iter)
  {
    std::cout << '\n' << indent << *vertex_iter << ':';
    for (boost::tie(edge_iter, edge_end) = boost::out_edges(*vertex_iter, graph); edge_iter != edge_end; ++edge_iter)
    {
      int matroid_element = boost::get(tu::edge_matroid_element, graph, *edge_iter);
      std::cout << ' ' << boost::target(*edge_iter, graph) << " (" << (matroid_element < 0 ? "row " : "column ") << matroid_element << ") ";
    }
  }
  std::cout << '\n';
}

void print_decomposition (const tu::decomposed_matroid* decomposition, std::string indent = "")
{
  if (decomposition->is_leaf())
  {
    tu::decomposed_matroid_leaf* leaf = (tu::decomposed_matroid_leaf*) (decomposition);

    if (leaf->is_R10())
    {
      std::cout << indent << "R10:";
      for (tu::matroid_element_set::const_iterator iter = leaf->elements().begin(); iter != leaf->elements().end(); ++iter)
        std::cout << " " << *iter;
      std::cout << "\n";
    }
    else if (leaf->is_graphic() && leaf->is_cographic())
    {
      std::cout << indent << "planar binary matroid.\n";
      std::cout << indent << "graph:\n" << indent << "{ ";
      print_matroid_graph(*leaf->graph(), indent + "  ");
      std::cout << indent << "}\n" << indent << "cograph:\n" << indent << "{ ";
      print_matroid_graph(*leaf->cograph(), indent + "  ");
      std::cout << indent << "}\n";
    }
    else if (leaf->is_graphic())
    {
      std::cout << indent << "graphic binary matroid.\n";
      std::cout << indent << "graph:\n" << indent << "{ ";
      print_matroid_graph(*leaf->graph(), indent + "  ");
      std::cout << indent << "}\n";
    }
    else if (leaf->is_cographic())
    {
      std::cout << indent << "cographic binary matroid.\n";
      std::cout << indent << "cograph:\n" << indent << "{ ";
      print_matroid_graph(*leaf->cograph(), indent + "  ");
      std::cout << indent << "}\n";
    }
    else
    {
      std::cout << indent << "irregular matroid.\n";
    }
  }
  else
  {
    tu::decomposed_matroid_separator* separator = (tu::decomposed_matroid_separator*) (decomposition);

    if (separator->separation_type() == tu::decomposed_matroid_separator::ONE_SEPARATION)
    {
      std::cout << indent << "1-separation:\n";

    }
    else if (separator->separation_type() == tu::decomposed_matroid_separator::TWO_SEPARATION)
    {
      std::cout << indent << "2-separation:\n";
    }
    else if (separator->separation_type() == tu::decomposed_matroid_separator::THREE_SEPARATION)
    {
      std::cout << indent << "3-separation:\n";
    }
    else
    {
      std::cout << indent << "invalid separation:\n";
    }
    std::cout << indent << "{\n";
    print_decomposition(separator->first(), indent + "  ");
    print_decomposition(separator->second(), indent + "  ");
    std::cout << indent << "}\n";
  }
}

void print_violator (const tu::integer_matrix& matrix, const tu::submatrix_indices& violator)
{
  typedef boost::numeric::ublas::matrix_indirect <const tu::integer_matrix, tu::submatrix_indices::indirect_array_type> indirect_matrix_t;

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

int run (const std::string& file_name, bool show_certificates, tu::log_level level)
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

  tu::integer_matrix matrix(height, width);
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
    tu::submatrix_indices violator;
    tu::decomposed_matroid* decomposition;

    if (tu::is_totally_unimodular(matrix, decomposition, violator, level))
    {
      std::cout << "\nThe " << matrix.size1() << " x " << matrix.size2() << " matrix is totally unimodular.\n" << std::endl;

      print_decomposition(decomposition);
    }
    else
    {
      std::cout << "\nThe " << matrix.size1() << " x " << matrix.size2() << " matrix is not totally unimodular." << std::endl;
      assert (violator.rows.size() == violator.columns.size());
      int det = tu::determinant_submatrix(matrix, violator);
      std::cout << "\nThe violating submatrix (det = " << det << ") is " << violator.rows.size() << " x " << violator.columns.size() << ":\n\n"
          << std::flush;
      print_violator(matrix, violator);
    }
    delete decomposition;
  }
  else
  {
    if (tu::is_totally_unimodular(matrix, level))
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

bool extract_option (char c, bool& certs, tu::log_level& level, bool& help)
{
  if (c == 'h')
    help = true;
  else if (c == 'c')
    certs = true;
  else if (c == 'q')
    level = tu::LOG_QUIET;
  else if (c == 'u')
    level = tu::LOG_UPDATING;
  else if (c == 'v')
    level = tu::LOG_VERBOSE;
  else
    return false;

  return true;
}

int main (int argc, char **argv)
{
  /// Possible parameters
  std::string matrix_file_name = "";
  bool certs = false;
  tu::log_level level = tu::LOG_UPDATING;
  bool help = false;

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
          if (!extract_option(current[i], certs, level, help))
          {
            std::cout << "Unknown option: -" << current[i] << "\nSee " << argv[0] << " -h for usage." << std::endl;
            return EXIT_FAILURE;
          }
        }
        continue;
      }
    }

    if (matrix_file_name != "")
    {
      std::cout << "Matrix file was given twice!\nSee " << argv[0] << " -h for usage." << std::endl;
      return EXIT_FAILURE;
    }
    matrix_file_name = current;
  }

  if (help)
  {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] [--] MATRIX_FILE\n";
    std::cout << "Options:\n";
    std::cout << " -h Shows a help message.\n";
    std::cout << " -c Prints certificates: A matroid decomposition if the matrix is totally unimodular and a violating submatrix otherwise.\n";
    std::cout << " -u Updating logging (default).\n";
    std::cout << " -v Verbose logging.\n";
    std::cout << " -q No logging at all.\n";
    std::cout << std::flush;
    return EXIT_SUCCESS;
  }

  if (matrix_file_name == "")
  {
    std::cout << "No matrix file was given!\nSee " << argv[0] << " -h for usage." << std::endl;
    return EXIT_FAILURE;
  }

  return run(matrix_file_name, certs, level);
}
