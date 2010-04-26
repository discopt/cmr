/*
 * total_unimodularity.cpp
 *
 *  Created on: Dec 20, 2009
 *      Author: xammy
 */

#include "../config.h"
#include "algorithm.hpp"
#include "total_unimodularity.hpp"
#include "matroid.hpp"
#include "violator_search.hpp"
#include "signing.hpp"

#include <boost/numeric/ublas/io.hpp>

namespace tu {

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
        std::cout << indent << "  matroid elements: ";
        std::copy (leaf->elements ().begin (), leaf->elements ().end (), std::ostream_iterator <int> (std::cout, " "));
        std::cout << "\n" << indent << "  extra elements: ";
        std::copy (leaf->extra_elements ().begin (), leaf->extra_elements ().end (), std::ostream_iterator <int> (std::cout, " "));
        std::cout << "\n";
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

  /// Returns a decomposition of a given binary matroid into a k-sum-decomposition (k=1,2,3)
  /// in graphic, cographic, R10 and maybe irregular components. 

  decomposed_matroid* decompose_binary_matroid (const boost::numeric::ublas::matrix <int>& input_matrix)
  {
    if (!is_zero_one_matrix (input_matrix))
      return NULL;

    integer_matrix matrix (input_matrix);
    integer_matroid matroid (matrix.size1 (), matrix.size2 ());

    return decompose_binary_matroid (matroid, matrix, matroid_element_set (), true).second;
  }

  /// Returns true, iff the given matrix is totally unimodular.

  bool is_totally_unimodular (const integer_matrix& input_matrix)
  {
    std::cout << "Checking for a -1/0/1 matrix" << std::endl;

    if (!is_zero_plus_minus_one_matrix (input_matrix))
    {
      return false;
    }

    std::cout << "Checking for signedness" << std::endl;

    if (!is_signed_matrix (input_matrix))
    {
      //        submatrix_indices submatrix;
      //        is_signed_matrix (input_matrix, submatrix);
      //
      //        for (size_t col = 0; col < submatrix.columns.size (); col++)
      //        {
      //            std::cout << "column " << col << " = " << submatrix.columns[col] << std::endl;
      //        }
      //        for (size_t row = 0; row < submatrix.rows.size (); row++)
      //        {
      //            std::cout << "row " << row << " = " << submatrix.rows[row] << std::endl;
      //        }
      //
      //        const boost::numeric::ublas::matrix_indirect<const boost::numeric::ublas::matrix<int>,
      //                submatrix_indices::indirect_array_type> indirect_matrix (input_matrix, submatrix.rows,
      //                submatrix.columns);
      //
      //        boost::numeric::ublas::matrix<float> matrix (indirect_matrix);
      //        
      //        std::cout << matrix << std::endl;

      return false;
    }

    std::cout << "Copying matrix and creating matroid" << std::endl;

    integer_matrix matrix (input_matrix);
    integer_matroid matroid (matrix.size1 (), matrix.size2 ());

    support_matrix (matrix);

    std::pair <bool, decomposed_matroid*> result = decompose_binary_matroid (matroid, matrix, matroid_element_set (), false);

    assert (result.second == NULL);

    return result.first;
  }

  /// Returns true, iff the given matrix is totally unimodular.
  /// decomposition points to a k-sum-decomposition (k=1,2,3) in graphic, cographic, R10 and maybe irregular components.

  bool is_totally_unimodular (const integer_matrix& input_matrix, decomposed_matroid*& decomposition)
  {
    if (!is_zero_plus_minus_one_matrix (input_matrix) || !is_signed_matrix (input_matrix))
    {
      decomposition = NULL;
      return false;
    }

    integer_matrix matrix (input_matrix);
    integer_matroid matroid (matrix.size1 (), matrix.size2 ());

    support_matrix (matrix);

    std::pair <bool, decomposed_matroid*> result = decompose_binary_matroid (matroid, matrix, matroid_element_set (), true);
    decomposition = result.second;

    return result.first;
  }

  /// Returns true, iff the given matrix is totally unimodular.
  /// decomposition points to a k-sum-decomposition (k=1,2,3) in graphic, cographic, R10 and maybe irregular components.

  bool is_totally_unimodular (const integer_matrix& input_matrix, decomposed_matroid*& decomposition, submatrix_indices& violator)
  {
    /// Check each entry 
    std::pair <unsigned int, unsigned int> entry;
    if (!is_zero_plus_minus_one_matrix (input_matrix, entry))
    {
      violator.rows = submatrix_indices::indirect_array_type (1);
      violator.rows[0] = entry.first;
      violator.columns = submatrix_indices::indirect_array_type (1);
      violator.columns[0] = entry.second;
      decomposition = NULL;
      return false;
    }

    /// Signing test
    if (!is_signed_matrix (input_matrix, violator))
    {
      decomposition = NULL;
      assert (violator.rows.size() == violator.columns.size());
      return false;
    }

    integer_matroid matroid (input_matrix.size1 (), input_matrix.size2 ());
    integer_matrix matrix (input_matrix);

    /// Remove sign from matrix
    support_matrix (matrix);

    /// Matroid decomposition
    bool is_tu;
    boost::tie (is_tu, decomposition) = decompose_binary_matroid (matroid, matrix, matroid_element_set (), true);
    if (is_tu)
      return true;

    matroid_element_set rows, columns, elements = detail::find_smallest_irregular_minor (decomposition);
    detail::split_elements (elements.begin (), elements.end (), std::inserter (rows, rows.end ()), std::inserter (columns, columns.end ()));

    detail::violator_strategy* strategy = new detail::single_violator_strategy (input_matrix, rows, columns);

    strategy->search ();
    strategy->create_matrix (violator);

    delete strategy;

    assert (violator.rows.size() == violator.columns.size());

    return false;
  }

  /// Returns true, iff the given matrix is totally unimodular.
  /// If this is not the case, violator describes a violating submatrix.  

  bool is_totally_unimodular (const integer_matrix& matrix, submatrix_indices& violator)
  {
    decomposed_matroid* decomposition;

    bool is_tu = is_totally_unimodular (matrix, decomposition, violator);

    delete decomposition;

    return is_tu;
  }

  /// Returns true, iff the given matrix is already a signed version of its support matrix.
  /// Running time: O(height * width * min (height, width) )

  bool is_signed_matrix (const boost::numeric::ublas::matrix <int>& matrix)
  {
    if (matrix.size2 () > matrix.size1 ())
    {
      const matrix_transposed <const boost::numeric::ublas::matrix <int> > transposed (matrix);
      return sign_matrix (transposed, NULL);
    }
    else
    {
      return sign_matrix (matrix, NULL);
    }
  }

  /// Returns true, iff the given matrix is already a signed version of its support matrix.
  /// If not, violator describes a violating submatrix.
  /// Running time: O(height * width * min (height, width) )

  bool is_signed_matrix (const boost::numeric::ublas::matrix <int>& matrix, submatrix_indices& violator)
  {
    if (matrix.size2 () > matrix.size1 ())
    {
      const matrix_transposed <const boost::numeric::ublas::matrix <int> > transposed (matrix);
      return sign_matrix (transposed, &violator);
    }
    else
    {
      return sign_matrix (matrix, &violator);
    }
  }

  /// Returns true, iff the given matrix is already a signed version of its support matrix.
  /// If not, the given matrix is changed to be such a signed version.
  /// Running time: O(height * width * min (height, width) )

  bool sign_matrix (boost::numeric::ublas::matrix <int>& matrix)
  {
    if (matrix.size2 () > matrix.size1 ())
    {
      matrix_transposed <boost::numeric::ublas::matrix <int> > transposed (matrix);
      return sign_matrix (transposed, NULL);
    }
    else
    {
      return sign_matrix (matrix, NULL);
    }
  }

  /// Drops all signs in the given matrix.

  void support_matrix (boost::numeric::ublas::matrix <int>& matrix)
  {
    for (size_t i = 0; i < matrix.size1 (); ++i)
    {
      for (size_t j = 0; j < matrix.size2 (); ++j)
      {
        const int value = matrix (i, j);
        if (value != 0)
        {
          matrix (i, j) = 1;
        }
      }
    }
  }

}
