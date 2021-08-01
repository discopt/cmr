/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#include <fstream>
#include <iomanip>
#include <map>

#include <tu/config.h>

#include <boost/logic/tribool.hpp>

#include <cmr/total_unimodularity.hpp>
#include <cmr/matroid_decomposition.hpp>
#include <cmr/unimodularity.hpp>
#include <cmr/smith_normal_form.hpp>

template <typename Set, typename Element>
bool contains(const Set& set, const Element& element)
{
  return set.find(element) != set.end();
}

void print_matroid_graph(const tu::matroid_graph& graph, const std::string& indent = "")
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

void print_decomposition(const tu::decomposed_matroid* decomposition, std::string indent = "")
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

void print_violator(const tu::integer_matrix& matrix, const tu::submatrix_indices& violator)
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

void print_result(std::ostream& stream, const std::string& name, boost::logic::tribool result)
{
  stream << std::setw(30) << name << ": ";

  if (result)
    stream << "yes\n";
  else if (!result)
    stream << "no\n";
  else
    stream << "not determined\n";
}

bool test_total_unimodularity(tu::integer_matrix& matrix, bool show_certificates, tu::log_level level)
{
  bool result;

  if (show_certificates)
  {
    tu::submatrix_indices violator;
    tu::decomposed_matroid* decomposition;

    result = tu::is_totally_unimodular(matrix, decomposition, violator, level);

    if (result)
    {
      std::cout << "\nThe matrix is totally unimodular due to the following decomposition:\n" << std::endl;

      print_decomposition(decomposition);
    }
    else
    {
      int det = tu::submatrix_determinant(matrix, violator);
      assert (violator.rows.size() == violator.columns.size());
      std::cout << "\nThe matrix is not totally unimodular due to the following " << violator.rows.size() << " x " << violator.columns.size()
          << " submatrix with determinant " << det << "." << std::endl;
      print_violator(matrix, violator);
    }
    delete decomposition;
  }
  else
  {
    result = tu::is_totally_unimodular(matrix, level);
    std::cout << "The matrix is " << (result ? "" : "not ") << "totally unimodular." << std::endl;
  }

  return result;
}

bool test_regularity(tu::integer_matrix& matrix, tu::log_level level)
{
  bool result;

  tu::integer_matrix copy = matrix;
  tu::sign_matrix(copy);
  result = tu::is_totally_unimodular(copy, level);
  std::cout << "The underlying binary matroid is " << (result ? "" : "not ") << "regular." << std::endl;

  return result;
}

int run(const std::string& file_name, const std::set <char>& tests, bool show_certificates, tu::log_level level)
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
    std::cout << "Error: cannot read matrix size from input file." << std::endl;
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

  std::cerr << "Combinatorial Matrix Recognition version " << CMR_VERSION_MAJOR << "." << CMR_VERSION_MINOR
    << CMR_VERSION_PATCH << " by Matthias Walter and Klaus Truemper.\n";
  std::cerr << "See http://matthiaswalter.org/CMRtest/ for references and citation.\n\n" << std::flush;

  std::cout << "Read a " << matrix.size1() << " x " << matrix.size2() << " matrix.\n" << std::endl;

  std::map <char, boost::logic::tribool> results;
  for (size_t i = 0; i < 7; ++i)
    results["tUuMmCr"[i]] = boost::logic::indeterminate;
  size_t rank = 0;
  bool know_rank = false;
  unsigned int k = 0;

  if (contains(tests, 's'))
  {
    tu::sign_matrix(matrix);

    std::cout << "Signed version of input matrix is the following.\n\n";
    std::cout << matrix.size1() << " " << matrix.size2() << "\n";
    for (std::size_t r = 0; r < matrix.size1(); ++r)
    {
      for (std::size_t c = 0; c < matrix.size2(); ++c)
        std::cout << ' ' << std::setw(2) << matrix(r,c);
      std::cout << "\n";
    }
    std::cout << std::endl;
  }

  if (contains(tests, 't'))
  {
    /// Let's test for total unimodularity.

    results['t'] = test_total_unimodularity(matrix, show_certificates, level);
    k = 1;
  }

  if (results['t'])
  {
    for (size_t i = 0; i < 5; ++i)
    {
      results["uUmMr"[i]] = true;
    }
    if (contains(tests, 'C'))
    {
       /// Test for complement total unimodularity.

      std::size_t complementedRow, complementedColumn;
      results['C'] = tu::is_complement_total_unimodular(matrix, complementedRow, complementedColumn, tu::LOG_QUIET);
      if (!results['C'] && show_certificates)
      {
        std::cout << "The matrix obtained by complementing ";
        if (complementedRow < matrix.size1())
          std::cout << "row " << complementedRow;
        else
          std::cout << "no row";
        std::cout << " and complementing ";
        if (complementedColumn < matrix.size2())
          std::cout << "column " << complementedColumn;
        else
          std::cout << "no column";
        std::cout << " is not totally unimodular." << std::endl;
      }
    }
  }
  else
  {
    if (!results['t'])
      results['C'] = false;

    if (contains(tests, 'm') || contains(tests, 'M'))
    {
      if (level != tu::LOG_QUIET)
        std::cout << "Testing matrix for k-modularity... " << std::flush;
      results['m'] = tu::is_k_modular(matrix, rank, k, tu::LOG_PROGRESSIVE);
      std::cout << "The matrix is " << (results['m'] ? "" : "not ") << "k-modular.\n" << std::flush;

      results['u'] = (results['m'] && k == 1);
      if (results['m'])
        std::cout << "The matrix is " << (results['u'] ? "" : "not ") << "unimodular.\n" << std::flush;

      if (!results['m'])
        results['M'] = false;
      if (!results['u'])
      {
        results['U'] = false;
        results['t'] = false;
      }
      know_rank = true;
    }
    else if (contains(tests, 'u') || contains(tests, 'U'))
    {
      if (level != tu::LOG_QUIET)
        std::cout << "Testing matrix for unimodularity... " << std::flush;
      results['u'] = tu::is_unimodular(matrix, rank, tu::LOG_QUIET);
      std::cout << "The matrix is " << (results['u'] ? "" : "not ") << "unimodular.\n" << std::flush;

      if (results['u'])
        results['m'] = true;
      if (!results['u'])
      {
        results['U'] = false;
        results['t'] = false;
      }
      know_rank = true;
    }
    if (contains(tests, 'M') && boost::logic::indeterminate(results['M']))
    {
      if (level != tu::LOG_QUIET)
        std::cout << "Testing transpose of matrix for k-modularity... " << std::flush;
      tu::matrix_transposed <tu::integer_matrix> transposed(matrix);
      results['M'] = tu::is_k_modular(transposed, rank, k, tu::LOG_QUIET);
      std::cout << "The transpose is " << (results['M'] ? "" : "not ") << "k-modular.\n" << std::flush;

      results['U'] = (results['M'] && k == 1);
      if (results['M'])
        std::cout << "The transpose is " << (results['U'] ? "" : "not ") << "unimodular.\n" << std::flush;

      if (!results['U'])
        results['t'] = false;
      know_rank = true;
    }
    else if (contains(tests, 'U') && boost::logic::indeterminate(results['U']))
    {
      if (level != tu::LOG_QUIET)
        std::cout << "Testing transpose of matrix for unimodularity... " << std::flush;
      tu::matrix_transposed <tu::integer_matrix> transposed(matrix);
      results['U'] = tu::is_unimodular(transposed, rank, tu::LOG_QUIET);
      std::cout << "The transpose is " << (results['U'] ? "" : "not ") << "unimodular.\n" << std::flush;

      if (!results['U'])
        results['t'] = false;
      know_rank = true;
    }
  }

  /// Print a summary

  if (know_rank)
    std::cout << "\nSummary of rank " << rank << " matrix:\n\n";
  else
    std::cout << "\nSummary:\n\n";

  print_result(std::cout, "Totally unimodular", results['t']);
  print_result(std::cout, "Strongly unimodular", results['U']);
  print_result(std::cout, "Complement total unimodularity", results['C']);
  print_result(std::cout, "Unimodular", results['u']);
  print_result(std::cout, "Strongly k-modular", results['M']);
  print_result(std::cout, "k-modular", results['m']);
  if (results['m'])
    std::cout << "                  k = " << k << "\n";
  if (know_rank)
    print_result(std::cout, "Dantzig property", results['m'] && rank == matrix.size1());

  std::cout << std::flush;

  return EXIT_SUCCESS;
}

bool extract_option(char c, std::set <char>& tests, bool& certs, tu::log_level& level, bool& help)
{
  if (c == 't' || c == 'u' || c == 'm' || c == 'U' || c == 'M' || c == 's')
    tests.insert(c);
  else if (c == 'C')
  {
    tests.insert('t');
    tests.insert('C');
  }
  else if (c == 'a')
  {
    tests.insert('t');
    tests.insert('u');
    tests.insert('U');
    tests.insert('C');
    tests.insert('m');
    tests.insert('M');
  }
  else if (c == 'h')
    help = true;
  else if (c == 'c')
    certs = true;
  else if (c == 'q')
    level = tu::LOG_QUIET;
  else if (c == 'p')
    level = tu::LOG_PROGRESSIVE;
  else if (c == 'v')
    level = tu::LOG_VERBOSE;
  else
    return false;

  return true;
}

int main(int argc, char **argv)
{
  /// Possible parameters
  std::string matrix_file_name = "";
  bool certs = false;
  tu::log_level level = tu::LOG_PROGRESSIVE;
  bool help = false;
  std::set <char> tests;

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
          if (!extract_option(current[i], tests, certs, level, help))
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

  if (tests.empty())
  {
    tests.insert('t');
    tests.insert('u');
    tests.insert('U');
    tests.insert('C');
    tests.insert('m');
    tests.insert('M');
  }

  if (help)
  {
    std::cerr << "Usage: " << argv[0] << " [OPTIONS] [--] MATRIX_FILE\n";
    std::cerr << "Options:\n";
    std::cerr << " -h Shows a help message.\n";
    std::cerr << " -s Signs the matrix beforehand.\n";
    std::cerr << " -a Test for everything possible (default).\n";
    std::cerr << " -t Test for total unimodularity.\n";
    std::cerr << " -U Test for strong unimodularity.\n";
    std::cerr << " -C Test for complement total unimodularity.\n";
    std::cerr << " -u Test for unimodularity.\n";
    std::cerr << " -M Test for strong k-modularity.\n";
    std::cerr << " -m Test for k-modularity.\n";
    std::cerr << " -c Prints certificates: Try to find certificates for the results.\n";
    std::cerr << " -p Progressive logging (default).\n";
    std::cerr << " -v Verbose logging.\n";
    std::cerr << " -q No logging at all.\n";
    std::cerr << std::flush;
    return EXIT_SUCCESS;
  }

  if (matrix_file_name == "")
  {
    std::cerr << "No matrix file was given!\nSee " << argv[0] << " -h for usage." << std::endl;
    return EXIT_FAILURE;
  }

  return run(matrix_file_name, tests, certs, level);
}
