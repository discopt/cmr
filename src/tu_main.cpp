#include <fstream>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/program_options/parsers.hpp>

#include "total_unimodularity.hpp"
#include "matroid_decomposition.hpp"

// TODO: Debug
#include <boost/numeric/ublas/io.hpp>

void print_matroid_graph (const tu::matroid_graph& graph, const std::string& indent = "")
{
  std::cout << boost::num_vertices (graph) << " nodes and " << boost::num_edges (graph) << " edges:";

  typedef boost::graph_traits <tu::matroid_graph> traits;
  traits::vertex_iterator vertex_iter, vertex_end;
  traits::out_edge_iterator edge_iter, edge_end;

  for (boost::tie (vertex_iter, vertex_end) = boost::vertices (graph); vertex_iter != vertex_end; ++vertex_iter)
  {
    std::cout << '\n' << indent << *vertex_iter << ':';
    for (boost::tie (edge_iter, edge_end) = boost::out_edges (*vertex_iter, graph); edge_iter != edge_end; ++edge_iter)
    {
      int matroid_element = boost::get (tu::edge_matroid_element, graph, *edge_iter);
      std::cout << ' ' << boost::target (*edge_iter, graph) << " (" << (matroid_element < 0 ? "row " : "column ") << (matroid_element - 1) << ") ";
    }
  }
  std::cout << '\n';
}

void print_decomposition (const tu::decomposed_matroid* decomposition, std::string indent = "")
{
  if (decomposition->is_leaf ())
  {
    tu::decomposed_matroid_leaf* leaf = (tu::decomposed_matroid_leaf*) (decomposition);

    if (leaf->is_R10 ())
      std::cout << indent << "R10\n";
    else if (leaf->is_planar ())
    {
      std::cout << indent << "planar binary matroid.\n";
      std::cout << indent << "graph:\n" << indent << "{ ";
      print_matroid_graph (*leaf->graph (), indent + "  ");
      std::cout << indent << "}\n" << indent << "cograph:\n" << indent << "{ ";
      print_matroid_graph (*leaf->cograph (), indent + "  ");
      std::cout << indent << "}\n";
    }
    else if (leaf->is_graphic ())
    {
      std::cout << indent << "graphic binary matroid.\n";
      std::cout << indent << "graph:\n" << indent << "{ ";
      print_matroid_graph (*leaf->graph (), indent + "  ");
      std::cout << indent << "}\n";
    }
    else if (leaf->is_cographic ())
    {
      std::cout << indent << "cographic binary matroid.\n";
      std::cout << indent << "cograph:\n" << indent << "{ ";
      print_matroid_graph (*leaf->cograph (), indent + "  ");
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

    if (separator->separation_type () == tu::decomposed_matroid_separator::ONE_SEPARATION)
    {
      std::cout << indent << "1-separation:\n";

    }
    else if (separator->separation_type () == tu::decomposed_matroid_separator::TWO_SEPARATION)
    {
      std::cout << indent << "2-separation:\n";
    }
    else if (separator->separation_type () == tu::decomposed_matroid_separator::THREE_SEPARATION)
    {
      std::cout << indent << "3-separation:\n";
    }
    else
    {
      std::cout << indent << "invalid separation:\n";
    }
    std::cout << indent << "{\n";
    print_decomposition (separator->first (), indent + "  ");
    print_decomposition (separator->second (), indent + "  ");
    std::cout << indent << "}\n";
  }
}

void print_violator (const tu::integer_matrix& matrix, const tu::submatrix_indices& violator)
{
  typedef boost::numeric::ublas::matrix_indirect <const tu::integer_matrix, tu::submatrix_indices::indirect_array_type> indirect_matrix_t;

  const indirect_matrix_t indirect_matrix (matrix, violator.rows, violator.columns);

  for (size_t row = 0; row < indirect_matrix.size1 (); ++row)
  {
    for (size_t column = 0; column < indirect_matrix.size2 (); ++column)
    {
      std::cout << std::setw (4) << indirect_matrix (row, column);
    }
    std::cout << '\n';
  }
  std::cout << "\n0-indexed rows:";
  for (size_t row = 0; row < violator.rows.size (); ++row)
    std::cout << ' ' << violator.rows[row];
  std::cout << "\n0-indexed columns:";
  for (size_t column = 0; column < violator.columns.size (); ++column)
    std::cout << ' ' << violator.columns[column];
  std::cout << std::endl;
}

namespace po = boost::program_options;

int run (const std::string& file_name, bool show_certificates)
{
  /// Open the file

  std::ifstream file (file_name.c_str ());
  if (!file.good ())
  {
    std::cout << "Error: cannot open file \"" << file_name << "\"." << std::endl;
    return EXIT_FAILURE;
  }

  /// Read height and width

  size_t height, width;
  file >> height >> width;
  if (!file.good ())
  {
    std::cout << "Error: cannot read height and width from input file." << std::endl;
    return EXIT_FAILURE;
  }

  /// Read matrix entries

  tu::integer_matrix matrix (height, width);
  for (size_t row = 0; row < height; ++row)
  {
    for (size_t column = 0; column < width; ++column)
    {
      if (!file.good ())
      {
        std::cout << "Error: cannot read matrix data." << std::endl;
      }
      int value;
      file >> value;
      matrix (row, column) = value;
    }
  }

  file.close ();

  if (show_certificates)
  {
    tu::submatrix_indices violator;
    tu::decomposed_matroid* decomposition;

    if (tu::is_totally_unimodular (matrix, decomposition, violator))
    {
      std::cout << "Matrix is totally unimodular.\n" << std::endl;

      print_decomposition (decomposition);
    }
    else
    {
      std::cout << "Matrix is not totally unimodular." << std::endl;
      int det = tu::determinant_submatrix (matrix, violator);
      std::cout << "\nThe violating submatrix (det = " << det << ") is " << violator.rows.size () << " x " << violator.columns.size () << ":\n\n";
      print_violator (matrix, violator);
    }
    delete decomposition;
  }
  else
  {
    if (tu::is_totally_unimodular (matrix))
    {
      std::cout << "Matrix is totally unimodular." << std::endl;
    }
    else
    {
      std::cout << "Matrix is not totally unimodular." << std::endl;
    }
  }

  return EXIT_SUCCESS;
}

int main (int argc, char **argv)
{
  po::options_description options_desc ("Allowed options");
  options_desc.add_options () ("help,h", "Shows a help message.") ("matrix,m", po::value <std::string> (),
      "Input matrix to test for total unimodularity.") ("certs,c",
      "Prints certificates: A matroid decomposition, if the matrix is totally unimodular and a violating submatrix otherwise.");
  po::positional_options_description positional_options_desc;
  positional_options_desc.add ("matrix", -1);
  po::variables_map variables;

  try
  {
    po::store (po::command_line_parser (argc, argv).options (options_desc).positional (positional_options_desc).run (), variables);
  }
  catch (po::error e)
  {
    std::cout << "Error: " << e.what () << ". Please use --help for usage." << std::endl;
    return EXIT_FAILURE;
  }
  po::notify (variables);

  if (variables.count ("help"))
  {
    std::cout << options_desc << "\n";
    return EXIT_SUCCESS;
  }

  if (!variables.count ("matrix"))
  {
    std::cout << "Error: you need to specify a matrix to test. Please use --help for usage." << std::endl;
    return EXIT_FAILURE;
  }

  return run (variables["matrix"].as <std::string> (), variables.count ("certs"));
}
