/*
 * algorithm.hpp
 *
 *  Created on: Dec 20, 2009
 *      Author: xammy
 */

#ifndef ALGORITHM_HPP_
#define ALGORITHM_HPP_

#include "../config.h"
#include "total_unimodularity.hpp"
#include "matroid_decomposition.hpp"
#include "find_wheel_minor.hpp"
#include "matroid_separate.hpp"
#include "nested_minor_sequence.hpp"
#include "find_minor_sequence.hpp"
#include "graphicness.hpp"
#include "r10.hpp"
#include "enumeration.hpp"

namespace tu {

  template <typename MatroidType, typename MatrixType>
  std::pair <bool, decomposed_matroid*> decompose_binary_matroid (MatroidType& matroid, MatrixType& matrix, matroid_element_set extra_elements,
      bool construct_decomposition);

  template <typename MatroidType, typename MatrixType>
  std::pair <bool, decomposed_matroid*> decompose_minor_sequence (matroid_permuted <MatroidType>& permuted_matroid,
      matrix_permuted <MatrixType>& permuted_matrix, nested_minor_sequence& nested_minors, matroid_element_set extra_elements,
      bool construct_decomposition)
  {
    std::cout << "Searching for sequence of nested minors in " << permuted_matroid.size1 () << " x " << permuted_matroid.size2 () << " matroid... "
        << std::flush;

    //    {
    //      integer_matrix copy = permuted_matrix;
    //      sign_matrix (copy);
    //      bool copy_result = ghouila_houri_is_totally_unimodular (copy);
    //      std::cout << "  [The matroid is " << (copy_result ? "" : "NOT ") << "regular.]  ";
    //    }

    separation sep = find_minor_sequence (permuted_matroid, permuted_matrix, nested_minors, extra_elements);
    if (sep.is_valid ())
    {
      std::cout << "found a " << (sep.rank () + 1) << "-separation instead." << std::endl;

      //      std::cout << "find_minor_sequence returned a separation of rank " << sep.rank () << std::endl;
      //      std::cout << "split is at " << sep.split ().first << "," << sep.split ().second << std::endl;
      //      std::cout << "witness is at " << sep.witness ().first << "," << sep.witness ().second << std::endl;
      //      if (sep.has_special_swap ())
      //      {
      //        if (sep.has_special_row_swap ())
      //          std::cout << "special row swap at " << sep.get_special_swap_index () << std::endl;
      //        else
      //          std::cout << "special column swap at " << sep.get_special_swap_index () << std::endl;
      //      }

      //      matroid_print (permuted_matroid, permuted_matrix);

      //      {
      //        integer_matrix copy = permuted_matrix;
      //        sign_matrix (copy);
      //        bool copy_result = ghouila_houri_is_totally_unimodular (copy);
      //        std::cout << "The matroid is " << (copy_result ? "" : "NOT ") << "regular.\n";
      //      }

      integer_matroid upper_left_matroid;
      integer_matrix upper_left_matrix;
      integer_matroid lower_right_matroid;
      integer_matrix lower_right_matrix;

      sep.create_components (permuted_matroid, permuted_matrix, upper_left_matroid, upper_left_matrix, lower_right_matroid, lower_right_matrix);

      //      std::cout << "Decomposed into\n";
      //      matroid_print (upper_left_matroid, upper_left_matrix);
      //      std::cout << "\nand\n";
      //      matroid_print (lower_right_matroid, lower_right_matrix);
      //      std::cout << std::endl;

      //      {
      //        integer_matrix copy = upper_left_matrix;
      //        sign_matrix (copy);
      //        bool copy_result = ghouila_houri_is_totally_unimodular (copy);
      //        std::cout << "The upper left matroid is " << (copy_result ? "" : "NOT ") << "regular.\n";
      //      }
      //
      //      {
      //        integer_matrix copy = lower_right_matrix;
      //        sign_matrix (copy);
      //        bool copy_result = ghouila_houri_is_totally_unimodular (copy);
      //        std::cout << "The lower right matroid is " << (copy_result ? "" : "NOT ") << "regular.\n";
      //      }

      // TODO: if rank == 0, filter some extra_elements!
      //      matroid_element_set upper_left_elements;
      //      matroid_element_set lower_right_elements;

      matroid_permuted <integer_matroid> permuted_upper_left_matroid (upper_left_matroid);
      matrix_permuted <integer_matrix> permuted_upper_left_matrix (upper_left_matrix);

      //        std::cout << "separation successful.\n             Looking at part 1\n" << std::endl;
      //        matroid_print (upper_left_matroid, upper_left_matrix);

      if (sep.has_special_swap ())
      {
        //            std::cout << "separation requests a special swap!" << std::endl;

        if (sep.has_special_row_swap ())
        {
          //                std::cout << "swapping " << (upper_left_matroid.size1 () - 1) << " and "
          //                        << sep.get_special_swap_index () << std::endl;
          matroid_permute1 (permuted_upper_left_matroid, permuted_upper_left_matrix, upper_left_matroid.size1 () - 1, sep.get_special_swap_index ());
        }
        else
        {
          matroid_permute2 (permuted_upper_left_matroid, permuted_upper_left_matrix, upper_left_matroid.size2 () - 1, sep.get_special_swap_index ());
        }

        //        std::cout << "Upper left after special swap" << std::endl;
        //        matroid_print (permuted_upper_left_matroid, permuted_upper_left_matrix);
      }

      std::pair <bool, decomposed_matroid*> upper_left_result = decompose_minor_sequence (permuted_upper_left_matroid, permuted_upper_left_matrix,
          nested_minors, extra_elements, construct_decomposition);
      //        decompose_binary_matroid (upper_left_matroid, upper_left_matrix, construct_decomposition);

      if (!construct_decomposition && !upper_left_result.first)
        return std::pair <bool, decomposed_matroid*> (false, NULL);

      //      std::cout << "\n               Looking at part 2\n" << std::endl;
      //      matroid_print (lower_right_matroid, lower_right_matrix);

      std::pair <bool, decomposed_matroid*> lower_right_result = decompose_binary_matroid (lower_right_matroid, lower_right_matrix, extra_elements,
          construct_decomposition);

      if (construct_decomposition)
      {
        int type = sep.rank () == 0 ? decomposed_matroid_separator::ONE_SEPARATION : decomposed_matroid_separator::TWO_SEPARATION;

        //        std::cout << "\nSTATUS OF 1/2-separation = " << lower_right_result.first << " && " << upper_left_result.first
        //            << "\n\n";

        return std::pair <bool, decomposed_matroid*> (lower_right_result.first && upper_left_result.first, new decomposed_matroid_separator (
            upper_left_result.second, lower_right_result.second, type, matroid_elements (permuted_matroid), extra_elements));
      }
      else
        return std::pair <bool, decomposed_matroid*> (lower_right_result.first, NULL);
    }

    std::cout << "done (W3 and " << nested_minors.size () << " extensions)." << std::endl;

    //    for (size_t i = 0; i < nested_minors.size (); i++)
    //    {
    //      std::cout << "extension " << i << " is " << nested_minor_sequence::get_extension_height (
    //          nested_minors.get_extension (i)) << " x " << nested_minor_sequence::get_extension_width (
    //          nested_minors.get_extension (i)) << " of type " << int(nested_minors.get_extension (i)) << std::endl;
    //    }

    //    if (permuted_matroid.size2 () < 100)
    //      matroid_print (permuted_matroid, permuted_matrix);
    //    else
    //      std::cout << "Omitting matrix of size " << permuted_matroid.size1 () << "," << permuted_matroid.size2 ()
    //          << std::endl;

    std::cout << "Checking for graphicness... " << std::flush;

    matroid_graph* graph = construct_matroid_graph (permuted_matroid, permuted_matrix, nested_minors);

    std::cout << (graph == NULL ? "not graphic." : "graphic.") << std::endl;

    if (!construct_decomposition && graph)
    {
      //      std::cout << "Matroid is graphic, thus we omit cographicness check!" << std::endl;

      //      std::cout << *graph << std::endl;

      delete graph;
      return std::make_pair (true, (decomposed_matroid*) NULL);
    }

    std::cout << "Checking for cographicness... " << std::flush;

    matroid_graph* cograph = construct_matroid_graph (view_matroid_transposed (permuted_matroid), view_matrix_transposed (permuted_matrix),
        view_nested_minor_sequence_transposed (nested_minors));

    std::cout << (cograph == NULL ? "not cographic." : "cographic.") << std::endl;

    if (!construct_decomposition && cograph)
    {
      //      std::cout << "Matroid is cographic" << std::endl;
      delete cograph;
      return std::make_pair (true, (decomposed_matroid*) NULL);
    }

    if (construct_decomposition && (graph || cograph))
    {
      return std::make_pair (true, new decomposed_matroid_leaf (graph, cograph, false, matroid_elements (permuted_matroid), extra_elements));
    }

    std::cout << "Checking for being R10... " << std::flush;

    if (is_r10 (permuted_matrix))
    {
      std::cout << "successful." << std::endl;
      if (construct_decomposition)
      {
        return std::make_pair (true, new decomposed_matroid_leaf (NULL, NULL, true, matroid_elements (permuted_matroid), extra_elements));
      }
      else
        return std::make_pair (true, (decomposed_matroid*) NULL);
    }
    else
    {
      std::cout << "failed." << std::endl;
    }

    sep = enumerate_separations (permuted_matroid, permuted_matrix, nested_minors, extra_elements);
    if (sep.is_valid ())
    {
      //      std::cout << "Decomposing a (3|4)-separation. at split " << sep.split ().first << " x " << sep.split ().second << std::endl;
      //      matroid_print (permuted_matroid, permuted_matrix);

      //      {
      //        integer_matrix copy = permuted_matrix;
      //        sign_matrix (copy);
      //        bool copy_result = ghouila_houri_is_totally_unimodular (copy);
      //        std::cout << "The matroid is " << (copy_result ? "" : "NOT ") << "regular." << std::endl;
      //      }

      integer_matroid upper_left_matroid;
      integer_matrix upper_left_matrix;
      integer_matroid lower_right_matroid;
      integer_matrix lower_right_matrix;

      sep.create_components (permuted_matroid, permuted_matrix, upper_left_matroid, upper_left_matrix, lower_right_matroid, lower_right_matrix);

      //      {
      //        integer_matrix copy = upper_left_matrix;
      //        sign_matrix (copy);
      //        bool copy_result = ghouila_houri_is_totally_unimodular (copy);
      //        std::cout << "The upper left matroid is " << (copy_result ? "" : "NOT ") << "regular." << std::endl;
      //      }
      //
      //      {
      //        integer_matrix copy = lower_right_matrix;
      //        sign_matrix (copy);
      //        bool copy_result = ghouila_houri_is_totally_unimodular (copy);
      //        std::cout << "The lower right matroid is " << (copy_result ? "" : "NOT ") << "regular." << std::endl;
      //      }

      matroid_permuted <integer_matroid> permuted_upper_left_matroid (upper_left_matroid);
      matrix_permuted <integer_matrix> permuted_upper_left_matrix (upper_left_matrix);

      //      std::cout << "Decomposed into\n";
      //      matroid_print (permuted_upper_left_matroid, permuted_upper_left_matrix);
      //      std::cout << "\nand\n";
      //      matroid_print (lower_right_matroid, lower_right_matrix);
      //      std::cout << std::endl;

      std::pair <bool, decomposed_matroid*> upper_left_result = decompose_binary_matroid (permuted_upper_left_matroid, permuted_upper_left_matrix,
          extra_elements, construct_decomposition);

      if (!construct_decomposition && !upper_left_result.first)
        return std::pair <bool, decomposed_matroid*> (false, NULL);
      //
      //      //        std::cout << "\n               Looking at part 2\n" << std::endl;
      //      //        matroid_print (lower_right_matroid, lower_right_matrix);
      //
      std::pair <bool, decomposed_matroid*> lower_right_result = decompose_binary_matroid (lower_right_matroid, lower_right_matrix, extra_elements,
          construct_decomposition);

      if (construct_decomposition)
      {
        //        std::cout << "\nSTATUS OF 3-separation = " << lower_right_result.first << " && " << upper_left_result.first
        //            << "\n\n";

        return std::make_pair (lower_right_result.first && upper_left_result.first, (decomposed_matroid *) (new decomposed_matroid_separator (
            upper_left_result.second, lower_right_result.second, decomposed_matroid_separator::THREE_SEPARATION, matroid_elements (permuted_matroid),
            extra_elements)));
      }
      else
        return std::pair <bool, decomposed_matroid*> (lower_right_result.first, NULL);
    }

    //    std::cout << "irregular minor:\n";
    //    matroid_print (permuted_matroid, permuted_matrix);

    if (construct_decomposition)
    {
      return std::make_pair (false, new decomposed_matroid_leaf (NULL, NULL, false, matroid_elements (permuted_matroid), extra_elements));
    }
    else
      return std::make_pair (false, (decomposed_matroid*) NULL);
  }

  /**
   * Create the graph for 2 x w or h x 2 matrices.
   * 
   * @param matroid
   * @param matrix
   * @return
   */

  template <typename MatroidType, typename MatrixType>
  matroid_graph* construct_small_matroid_graph (MatroidType& matroid, MatrixType& matrix)
  {
    assert (matroid.size1() <= 2 || matroid.size2() <= 2);

    //    std::cout << "construct_small_matroid_graph (" << matroid.size1 () << ", " << matroid.size2 () << ")" << std::endl;

    matroid_graph* graph = new matroid_graph (matroid.size1 () + 1);

    if (matroid.size1 () == 0)
    {
      for (size_t column = 0; column < matroid.size2 (); ++column)
      {
        boost::add_edge (boost::vertex (0, *graph), boost::vertex (0, *graph), matroid.name2 (0), *graph);
      }
    }
    else if (matroid.size2 () == 0)
    {
      for (size_t row = 0; row < matroid.size1 (); ++row)
      {
        boost::add_edge (boost::vertex (row, *graph), boost::vertex (row + 1, *graph), matroid.name1 (row), *graph);
      }
    }
    else if (matroid.size1 () >= 3 && matroid.size2 () == 1)
    {
      size_t current_edge_vertex = 0;
      size_t current_free_vertex = matroid.size1 ();

      for (size_t row = 0; row < matroid.size1 (); ++row)
      {
        if (matrix (row, 0) == 0)
        {
          boost::add_edge (boost::vertex (current_free_vertex - 1, *graph), boost::vertex (current_free_vertex, *graph), matroid.name1 (row), *graph);
          --current_free_vertex;
        }
        else
        {
          boost::add_edge (boost::vertex (current_edge_vertex, *graph), boost::vertex (current_edge_vertex + 1, *graph), matroid.name1 (row), *graph);
          ++current_edge_vertex;
        }
      }

      assert (current_edge_vertex == current_free_vertex);

      boost::add_edge (boost::vertex (0, *graph), boost::vertex (current_edge_vertex, *graph), matroid.name2 (0), *graph);
    }
    else if (matroid.size1 () >= 3 && matroid.size2 () == 2)
    {
      size_t count[4] = { 0, 0, 0, 0 };
      for (size_t row = 0; row < matroid.size1 (); ++row)
      {
        count[((matrix (row, 0) != 0) ? 2 : 0) + ((matrix (row, 1) != 0) ? 1 : 0)]++;
      }

      size_t vertex[4];
      vertex[0] = 0;
      vertex[1] = count[1];
      vertex[2] = vertex[1] + count[3];
      vertex[3] = vertex[2] + count[2];

      boost::add_edge (boost::vertex (vertex[0], *graph), boost::vertex (vertex[2], *graph), matroid.name2 (1), *graph);
      boost::add_edge (boost::vertex (vertex[1], *graph), boost::vertex (vertex[3], *graph), matroid.name2 (0), *graph);

      for (size_t row = 0; row < matroid.size1 (); ++row)
      {
        int element = matroid.name1 (row);
        switch (((matrix (row, 0) != 0) ? 2 : 0) + ((matrix (row, 1) != 0) ? 1 : 0))
        {
        case 1: /// 0 1
          boost::add_edge (boost::vertex (vertex[0], *graph), boost::vertex (vertex[0] + 1, *graph), element, *graph);
          ++vertex[0];
        break;
        case 2: /// 1 0
          boost::add_edge (boost::vertex (vertex[2], *graph), boost::vertex (vertex[2] + 1, *graph), element, *graph);
          ++vertex[2];
        break;
        case 3: /// 1 1
          boost::add_edge (boost::vertex (vertex[1], *graph), boost::vertex (vertex[1] + 1, *graph), element, *graph);
          ++vertex[1];
        break;
        default:
          assert (matrix (row, 0) ==0 && matrix (row, 1) == 0);
          boost::add_edge (boost::vertex (vertex[3], *graph), boost::vertex (vertex[3] + 1, *graph), element, *graph);
          ++vertex[3];
        }
      }
    }
    else
    {
      assert (matroid.size1() <= 2);

      for (size_t row = 0; row < matroid.size1 (); ++row)
      {
        boost::add_edge (boost::vertex (row, *graph), boost::vertex (row + 1, *graph), matroid.name1 (row), *graph);
      }

      for (size_t column = 0; column < matroid.size2 (); ++column)
      {
        size_t end_vertex, start_vertex = 0;
        while (start_vertex < matroid.size1 () && matrix (start_vertex, column) == 0)
        {
          ++start_vertex;
        }
        end_vertex = start_vertex;
        while (end_vertex < matroid.size1 () && matrix (end_vertex, column) == 1)
        {
          ++end_vertex;
        }

        boost::add_edge (boost::vertex (start_vertex, *graph), boost::vertex (end_vertex, *graph), matroid.name2 (column), *graph);
      }
    }

    return graph;
  }

  /**
   * 
   * 
   * @param matroid
   * @param matrix
   * @param construct_decomposition
   * @return
   */

  template <typename MatroidType, typename MatrixType>
  std::pair <bool, decomposed_matroid*> decompose_binary_matroid (MatroidType& matroid, MatrixType& matrix, matroid_element_set extra_elements,
      bool construct_decomposition)
  {
    std::cout << "Starting decomposition of binary matroid of size " << matrix.size1 () << " x " << matrix.size2 () << std::endl;
    //    matroid_print (matroid, matrix);

    //    {
    //      integer_matrix copy = matrix;
    //      sign_matrix (copy);
    //      bool copy_result = ghouila_houri_is_totally_unimodular (copy);
    //      std::cout << "The matroid is " << (copy_result ? "" : "NOT ") << "regular." << std::endl;
    //    }

    if (matroid.size1 () <= 2 || matroid.size2 () <= 2)
    {
      //      std::cout << "Very small matroid." << std::endl;
      if (construct_decomposition)
      {
        matroid_transposed <MatroidType> transposed_matroid (matroid);
        matrix_transposed <MatrixType> transposed_matrix (matrix);
        matroid_graph* g = construct_small_matroid_graph (matroid, matrix);
        matroid_graph* c = construct_small_matroid_graph (transposed_matroid, transposed_matrix);

        return std::make_pair (true, new decomposed_matroid_leaf (g, c, false, matroid_elements (matroid), extra_elements));
      }
      else
      {
        return std::make_pair (true, (decomposed_matroid*) NULL);
      }
    }

    typedef matroid_permuted <MatroidType> permuted_matroid_type;
    typedef matrix_permuted <MatrixType> permuted_marix_type;

    permuted_matroid_type permuted_matroid (matroid);
    permuted_marix_type permuted_matrix (matrix);

    // Identifies a W_3 minor in the upper left corner or finds a separation
    std::cout << "Searching for W3 minor in " << permuted_matroid.size1 () << " x " << permuted_matroid.size2 () << " matroid... " << std::flush;

    //    {
    //      integer_matrix copy = permuted_matrix;
    //      sign_matrix (copy);
    //      bool copy_result = ghouila_houri_is_totally_unimodular (copy);
    //      std::cout << "  [The matroid is " << (copy_result ? "" : "NOT ") << "regular.]  ";
    //    }

    separation sep = find_wheel_minor (permuted_matroid, permuted_matrix, extra_elements);

    //    {
    //      integer_matrix copy = permuted_matrix;
    //      sign_matrix (copy);
    //      bool copy_result = ghouila_houri_is_totally_unimodular (copy);
    //      std::cout << "  [The matroid is " << (copy_result ? "" : "NOT ") << "regular.]  ";
    //    }

    if (sep.is_valid ()) // separation

    {
      std::cout << "found a " << (sep.rank () + 1) << "-separation instead." << std::endl;

      //        std::cout << "find_wheel_minor found a separation" << std::endl;
      //        std::cout << "upper-left is " << sep.split ().first << " x " << sep.split ().second << std::endl;
      //        if (sep.rank () == 1)
      //        {
      //            separation::witness_type w = sep.witness ();
      //            std::cout << "1-separation witness is at " << w.first << "," << w.second << ", which is a "
      //                    << permuted_matrix (w.first, w.second) << std::endl;
      //        }
      //        matroid_print (permuted_matroid, permuted_matrix);

      integer_matroid upper_left_matroid;
      integer_matrix upper_left_matrix;
      integer_matroid lower_right_matroid;
      integer_matrix lower_right_matrix;

      sep.create_components (permuted_matroid, permuted_matrix, upper_left_matroid, upper_left_matrix, lower_right_matroid, lower_right_matrix);

      //      std::cout << "separation successful. Looking at part 1\n" << std::endl;
      //      matroid_print (upper_left_matroid, upper_left_matrix);

      std::pair <bool, decomposed_matroid*> upper_left_result = decompose_binary_matroid (upper_left_matroid, upper_left_matrix, extra_elements,
          construct_decomposition);

      if (!construct_decomposition && !upper_left_result.first)
        return std::pair <bool, decomposed_matroid*> (false, NULL);

      //      std::cout << "\nLooking at part2\n" << std::endl;
      //      matroid_print (lower_right_matroid, lower_right_matrix);

      std::pair <bool, decomposed_matroid*> lower_right_result = decompose_binary_matroid (lower_right_matroid, lower_right_matrix, extra_elements,
          construct_decomposition);

      if (construct_decomposition)
      {
        int type = sep.rank () == 0 ? decomposed_matroid_separator::ONE_SEPARATION : decomposed_matroid_separator::TWO_SEPARATION;

        //        std::cout << "\nSTATUS OF 1-2-separation = " << lower_right_result.first << " && " << upper_left_result.first
        //            << "\n\n";

        return std::pair <bool, decomposed_matroid*> (lower_right_result.first && upper_left_result.first, new decomposed_matroid_separator (
            upper_left_result.second, lower_right_result.second, type, matroid_elements (permuted_matroid), extra_elements));
      }
      else
        return std::pair <bool, decomposed_matroid*> (lower_right_result.first && upper_left_result.first, NULL);
    }

    std::cout << "done." << std::endl;

    //    std::cout << "W3-minor should now be visible.\n";
    //    matroid_print (permuted_matroid, permuted_matrix);

    nested_minor_sequence nested_minors;

    return decompose_minor_sequence (permuted_matroid, permuted_matrix, nested_minors, extra_elements, construct_decomposition);
  }

}

#endif /* ALGORITHM_HPP_ */
