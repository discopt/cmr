/**
 *          Copyright Matthias Walter 2010.
 * Distributed under the Boost Software License, Version 1.0.
 *    (See accompanying file LICENSE_1_0.txt or copy at
 *          http://www.boost.org/LICENSE_1_0.txt)
 **/

#ifndef ENUMERATION_HPP_
#define ENUMERATION_HPP_

#include <vector>

#include "separation.hpp"
#include "permutations.hpp"
#include "comparators.hpp"
#include "partition.hpp"
#include "binary_linear_space.hpp"
#include "bipartite_graph_bfs.hpp"
#include "matrix_modified.hpp"
#include "logger.hpp"

namespace tu {
  namespace detail {

    /**
     * This functor is a modifier for a modified matrix. It takes a split point.
     * The part left and below from the point is made zero.
     */

    struct empty_bottom_left_modifier
    {
      typedef int value_type;

      /**
       * Constructs the modifier.
       *
       * @param split A given split point for this modifier
       */

      empty_bottom_left_modifier (size_pair_t split) :
        _split(split)
      {

      }

      /**
       * Semantic of this functor: Make the bottom left part zero.
       *
       * @param i A row
       * @param j A column
       * @param value The value of the corresponding entry.
       * @return The filtered value
       */

      value_type operator() (size_t i, size_t j, value_type value)
      {
        if (i < _split.first || j >= _split.second)
          return value;
        else
          return 0;
      }

    private:
      size_pair_t _split;

    };

    /**
     * Finds paths in the top-left BFS graph between two bottom-left independet columns,
     * shortens them and swaps rows and columns such that the columns and the connecting
     * row are near the split point.
     *
     * @param matroid Given matroid
     * @param matrix A representation matrix for the matroid
     * @param split Split of a given (3|4)-separation.
     */

    template <typename MatroidType, typename MatrixType>
    inline void find_column_witnesses (MatroidType& matroid, MatrixType& matrix, size_pair_t split, matroid_element_set& extra_elements)
    {
      binary_linear_space vector_space(matroid.size1() - split.first);
      binary_linear_space::vector_type vector(matroid.size1() - split.first);
      std::vector <size_t> type(split.second);

      typedef std::set <size_t> column_set_t;
      column_set_t types[4];

      /// Set up matrix for BFS
      empty_bottom_left_modifier modifier(split);
      matrix_modified <MatrixType, empty_bottom_left_modifier> modified_matrix(matrix, modifier);
      bipartite_graph_dimensions dimensions(modified_matrix.size1(), modified_matrix.size2());

      for (size_t column = 0; column < split.second; ++column)
      {
        copy_partial_column(matrix, vector, column, split.first, matrix.size1());
        vector_space.insert_checked(vector);
        size_t type = vector_space.find(vector);

        assert (type < 4);
        types[type].insert(dimensions.column_to_index(column));
      }

      size_t reached_column = 0;
      size_t starting_column = 0;
      size_t connecting_row = 0;

      /// Try to find path from type 1 to the 2 and 3.

      size_t reached_node;
      std::vector <bipartite_graph_bfs_node> nodes;
      column_set_t start_nodes = types[1];
      column_set_t end_nodes = types[2];
      std::copy(types[3].begin(), types[3].end(), std::inserter(end_nodes, end_nodes.end()));

      if (bipartite_graph_bfs(modified_matrix, dimensions, start_nodes, end_nodes, false, nodes))
      {
        /// Find reached node
        reached_node = 0;
        for (column_set_t::const_iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter)
        {
          if (nodes[*iter].is_reachable())
          {
            reached_node = *iter;
            break;
          }
        }
      }
      else if (bipartite_graph_bfs(modified_matrix, dimensions, types[2], types[3], false, nodes))
      {
        /// Find reached node
        reached_node = 0;
        for (column_set_t::const_iterator iter = types[3].begin(); iter != types[3].end(); ++iter)
        {
          if (nodes[*iter].is_reachable())
          {
            reached_node = *iter;
            break;
          }
        }
      }
      else
      {
        throw std::runtime_error("tu::find_column_witnesses failed to find paths between columns of different type.");
      }

      /// We now have a path between two columns of different type and apply the path shortening technique to make it short.
      reached_column = dimensions.index_to_column(reached_node);
      for (size_t current_node = nodes[reached_node].predecessor; nodes[current_node].distance > 0; current_node = nodes[current_node].predecessor)
      {
        size_t next_node = nodes[current_node].predecessor;
        if (nodes[current_node].distance > 1 && dimensions.is_row(current_node))
        {
          assert (dimensions.is_column(next_node));

          matroid_binary_pivot(matroid, matrix, dimensions.index_to_row(current_node), dimensions.index_to_column(next_node));
          extra_elements.insert(matroid.name1(dimensions.index_to_row(current_node)));
          extra_elements.insert(matroid.name2(dimensions.index_to_column(next_node)));
        }
        else if (nodes[current_node].distance == 1)
        {
          connecting_row = dimensions.index_to_row(current_node);
          starting_column = dimensions.index_to_column(next_node);
        }
      }

      /// Swap everything to be near split point.
      matroid_permute1(matroid, matrix, connecting_row, split.first - 1);
      if (starting_column > reached_column)
        std::swap(starting_column, reached_column);
      matroid_permute2(matroid, matrix, reached_column, split.second - 1);
      matroid_permute2(matroid, matrix, starting_column, split.second - 2);
    }

    /**
     * Finds paths in the bottom-right BFS graph between two bottom-left independet rows,
     * shortens them and swaps rows and columns such that the rows and the connecting
     * column are near the split point.
     *
     * @param matroid Given matroid
     * @param matrix A representation matrix for the matroid
     * @param split Split of a given (3|4)-separation.
     */

    template <typename MatroidType, typename MatrixType>
    inline void find_row_witnesses (MatroidType& matroid, MatrixType& matrix, size_pair_t split, matroid_element_set& extra_elements)
    {
      binary_linear_space vector_space(split.second);
      binary_linear_space::vector_type vector(split.second);
      std::vector <size_t> type(matroid.size1() - split.first);

      typedef std::set <size_t> row_set_t;
      row_set_t types[4];

      /// Set up matrix for BFS
      empty_bottom_left_modifier modifier(split);
      matrix_modified <MatrixType, empty_bottom_left_modifier> modified_matrix(matrix, modifier);
      bipartite_graph_dimensions dimensions(modified_matrix.size1(), modified_matrix.size2());

      for (size_t row = split.first; row < matrix.size1(); ++row)
      {
        copy_partial_row(matrix, vector, row, 0, split.second);
        vector_space.insert_checked(vector);
        size_t type = vector_space.find(vector);
        assert (type < 4);
        types[type].insert(dimensions.row_to_index(row));
      }

      size_t reached_node;
      size_t reached_row = 0;
      size_t starting_row = 0;
      size_t connecting_column = 0;

      /// Try to find path from type 1 to the 2 and 3.

      std::vector <bipartite_graph_bfs_node> nodes;
      row_set_t start_nodes = types[1];
      row_set_t end_nodes = types[2];
      std::copy(types[3].begin(), types[3].end(), std::inserter(end_nodes, end_nodes.end()));

      if (bipartite_graph_bfs(modified_matrix, dimensions, start_nodes, end_nodes, false, nodes))
      {
        /// Find reached node
        reached_node = 0;
        for (row_set_t::const_iterator iter = end_nodes.begin(); iter != end_nodes.end(); ++iter)
        {
          if (nodes[*iter].is_reachable())
          {
            reached_node = *iter;
            break;
          }
        }
      }
      else if (bipartite_graph_bfs(modified_matrix, dimensions, types[2], types[3], false, nodes))
      {
        /// Find reached node
        reached_node = 0;
        for (row_set_t::const_iterator iter = types[3].begin(); iter != types[3].end(); ++iter)
        {
          if (nodes[*iter].is_reachable())
          {
            reached_node = *iter;
            break;
          }
        }
      }
      else
      {
        throw std::runtime_error("tu::find_row_witnesses failed to find paths between rows of different type.");
      }

      /// We now have a path between two columns of different type and apply the path shortening technique to make it short.
      reached_row = dimensions.index_to_row(reached_node);
      for (size_t current_node = nodes[reached_node].predecessor; nodes[current_node].distance > 0; current_node = nodes[current_node].predecessor)
      {
        size_t next_node = nodes[current_node].predecessor;
        if (nodes[current_node].distance > 1 && dimensions.is_column(current_node))
        {
          assert (dimensions.is_row(next_node));

          matroid_binary_pivot(matroid, matrix, dimensions.index_to_row(next_node), dimensions.index_to_column(current_node));
          extra_elements.insert(matroid.name1(dimensions.index_to_row(next_node)));
          extra_elements.insert(matroid.name2(dimensions.index_to_column(current_node)));
        }
        else if (nodes[current_node].distance == 1)
        {
          connecting_column = dimensions.index_to_column(current_node);
          starting_row = dimensions.index_to_row(next_node);
        }
      }

      matroid_permute2(matroid, matrix, connecting_column, split.second);
      if (starting_row < reached_row)
        std::swap(starting_row, reached_row);
      matroid_permute1(matroid, matrix, reached_row, split.first);
      matroid_permute1(matroid, matrix, starting_row, split.first + 1);
    }

    /**
     * Finds witnesses for a 3-separation.
     *
     * @param matroid A given matroid
     * @param matrix A representation matrix for the matroid
     * @param split Split of the 3-separation
     * @param extra_elements Set of matroid elements to be filled with pivot-elements
     * @return The final separation
     */

    template <typename MatroidType, typename MatrixType>
    inline separation find_witnesses (MatroidType& matroid, MatrixType& matrix, size_pair_t split, matroid_element_set& extra_elements)
    {
      find_column_witnesses(matroid, matrix, split, extra_elements);
      find_row_witnesses(matroid, matrix, split, extra_elements);

      if (matrix(split.first, split.second - 2) == 0 || matrix(split.first + 1, split.second - 1) == 0)
        matroid_permute1(matroid, matrix, split.first, split.first + 1);

      return separation(split, separation::witness_type(split.first, split.second - 2), separation::witness_type(split.first + 1, split.second - 1));
    }

    /**
     * Does a pivot in top-right part of a given (3|4)-separation to shift the rank.
     *
     * @param matroid Given matroid
     * @param matrix Represenation matrix of the given matroid
     * @param split Split of a given (3|4)-separation
     */

    template <typename MatroidType, typename MatrixType>
    inline void pivot_top_right (MatroidType& matroid, MatrixType& matrix, size_pair_t& split, matroid_element_set& extra_elements)
    {
      for (size_t row = 0; row < split.first; ++row)
      {
        for (size_t column = split.second; column < matrix.size2(); ++column)
        {
          if (matrix(row, column) != 0)
          {
            matroid_binary_pivot(matroid, matrix, row, column);
            extra_elements.insert(matroid.name1(row));
            extra_elements.insert(matroid.name2(column));
            matroid_permute1(matroid, matrix, row, split.first - 1);
            matroid_permute2(matroid, matrix, column, split.second);
            split.first--;
            split.second++;
            return;
          }
        }
      }
    }

    /**
     * Does a pivot in bottom-left part of a given (3|4)-separation to shift the rank.
     *
     * @param matroid Given matroid
     * @param matrix A representation matrix for the matroid
     * @param split Split of a given (3|4)-separation
     */

    template <typename MatroidType, typename MatrixType>
    inline void pivot_bottom_left (MatroidType& matroid, MatrixType& matrix, size_pair_t& split, matroid_element_set& extra_elements)
    {
      for (size_t row = split.first; row < matrix.size1(); ++row)
      {
        for (size_t column = 0; column < split.second; ++column)
        {
          if (matrix(row, column) != 0)
          {
            matroid_binary_pivot(matroid, matrix, row, column);
            extra_elements.insert(matroid.name1(row));
            extra_elements.insert(matroid.name2(column));
            matroid_permute1(matroid, matrix, row, split.first);
            matroid_permute2(matroid, matrix, column, split.second - 1);
            split.first++;
            split.second--;
            return;
          }
        }
      }
    }

    /**
     * Modifies a (3|4)-separation to have rank 2 in bottom-left area and top-right
     * non-degenerate.
     *
     * @param matroid Given matroid
     * @param matrix A representation matrix for the matroid
     * @param split Split of a given (3|4)-separation
     * @param ranks Rank distribution of the given separation
     */

    template <typename MatroidType, typename MatrixType>
    inline void normalize_3_4_separation (MatroidType& matroid, MatrixType& matrix, size_pair_t& split, rank_distribution ranks,
        matroid_element_set& extra_elements)
    {
      if (ranks == RANK_BL_2_TR_0)
        pivot_bottom_left(matroid, matrix, split, extra_elements);
      else if (ranks == RANK_BL_0_TR_2)
        pivot_top_right(matroid, matrix, split, extra_elements);

      /// If top-right has few rows, mirror the matrix, swapping bottom-left with top-right
      if (2 * split.first < matroid.size1())
      {
        for (size_t row = 0; row < matroid.size1() / 2; ++row)
          matroid_permute1(matroid, matrix, row, matroid.size1() - 1 - row);
        for (size_t column = 0; column < matroid.size2() / 2; ++column)
          matroid_permute2(matroid, matrix, column, matroid.size2() - 1 - column);
        split = size_pair_t(matroid.size1() - split.first, matroid.size2() - split.second);
      }

      /// Now as top-right has enough rows, we can pivot without degenerating it.
      pivot_top_right(matroid, matrix, split, extra_elements);
    }

    /**
     * Calls the partition algorithm to extend a given separation of a minor to the whole matroid.
     *
     * @param matroid The original matrix
     * @param matrix A representation matrix for the original matroid
     * @param worker_matrix Matrix for the partition algorithm to work on
     * @param top_left_size Size of the first part of the minor-separation
     * @param bottom_right_size Size of the second part of the minor-separation
     * @param separation Separation to return
     * @param extra_elements Set of matroid elements to be filled with pivot-elements
     * @return true if and only if The partitioning was successful
     */

    template <typename MatroidType, typename MatrixType>
    inline bool extend_to_3_4_separation (MatroidType& matroid, MatrixType& matrix, matrix_permuted <const integer_matrix>& worker_matrix,
        size_pair_t top_left_size, size_pair_t bottom_right_size, separation& separation, matroid_element_set& extra_elements)
    {
      rank_distribution ranks =
          partition(worker_matrix, top_left_size.first, top_left_size.second, bottom_right_size.first, bottom_right_size.second);
      if (ranks == RANK_TOO_HIGH)
        return false;

      if (top_left_size.first + top_left_size.second < 4 || bottom_right_size.first + bottom_right_size.second < 4)
      {
        return false;
      }

      /// Now we really found a separation. We apply the worker-matrix-permutation to the orginal matrix.

      permutation row_permutation = matrix.perm1() * worker_matrix.perm1();
      permutation column_permutation = matrix.perm2() * worker_matrix.perm2();

      matroid.perm1() = row_permutation;
      matrix.perm1() = row_permutation;
      matroid.perm2() = column_permutation;
      matrix.perm2() = column_permutation;

      assert (matrix_equals(matrix, worker_matrix));

      /// Normalize it

      normalize_3_4_separation(matroid, matrix, top_left_size, ranks, extra_elements);

      /// Make the witnesses visible.

      separation = find_witnesses(matroid, matrix, top_left_size, extra_elements);

      return true;
    }

    /**
     * Creates a permutation which sorts a given vector. It further counts the
     * number of positive and negativ values.
     *
     * @param permutation Permutation to return
     * @param mapping The given vector
     * @return (positive entries, negativ entries)
     */

    template <typename MappingValue>
    inline size_pair_t apply_mapping (permutation& permutation, std::vector <MappingValue>& mapping)
    {
      vector_less <MappingValue> less(mapping);
      sort(permutation, less);

      size_t negative = 0;
      while (negative < permutation.size() && mapping[permutation(negative)] < 0)
        ++negative;

      size_t positive = 0;
      while (positive < permutation.size() && mapping[permutation(permutation.size() - 1 - positive)] > 0)
        ++positive;

      return std::make_pair(negative, positive);
    }

    /**
     * Finds the index of the first minor in a given sequence of nested minors which has at least 8 elements.
     *
     * @param nested_minors Given sequence of nested minors
     * @param minor_size Size of the minor found
     * @return Index of the minor found
     */

    template <typename NestedMinorSequence>
    inline size_t find_first_8_element_minor (const NestedMinorSequence& nested_minors, size_pair_t& minor_size)
    {
      size_t current_height = 3;
      size_t current_width = 3;
      for (size_t i = 0; i < nested_minors.size(); ++i)
      {
        if (current_height + current_width >= 8)
        {
          minor_size = std::make_pair(current_height, current_width);
          return i;
        }
        current_height += nested_minors.get_extension_height(i);
        current_width += nested_minors.get_extension_width(i);
      }

      assert (current_height + current_width >= 8);

      minor_size = std::make_pair(current_height, current_width);
      return nested_minors.size();
    }

    /**
     * Enumerates all necessary partitions for a fixed extension inside a sequence of nested minors.
     * It works on a separate matrix with mapping vectors. Those contain either entries +1 and -1
     * to show to which part of the minor-separation a row or column belongs and a 0 to designate
     * a row or column to be unspecified.
     *
     * @param matroid The given matroid
     * @param matrix Representation matrix for the given matroid
     * @param worker_matrix Matrix to work on with the right permutation
     * @param row_mapping Mapping vector for rows
     * @param column_mapping Mapping vector for columns
     * @param minor_size Size of the current minor in the sequence of nested minors
     * @param ext_height Height of the extension
     * @param ext_width Width of the extension
     * @param separation Separation to be returned
     * @param extra_elements Set of matroid-elements to be filled with pivot-elements
     * @param enumeration Iteration step in the whole enumeration.
     * @param next_enumeration Next enumeration
     * @param max_enumerations Total number of enumerations
     * @param next_percent Next percent value for a percent-jump
     * @param log Logger
     * @param cut Length of the line to cut off the progress-part
     * @return true if and only if the partitioning algorithm was successful
     */

    template <typename MatroidType, typename MatrixType, typename MappingValue>
    inline bool enumerate_extension (MatroidType& matroid, MatrixType& matrix, matrix_permuted <const integer_matrix>& worker_matrix, std::vector <
        MappingValue>& row_mapping, std::vector <MappingValue>& column_mapping, size_pair_t minor_size, size_t ext_height, size_t ext_width,
        separation& separation, matroid_element_set& extra_elements, unsigned long long& enumeration, unsigned long long& next_enumeration,
        unsigned long long max_enumerations, unsigned int& next_percent, logger& log, size_t cut)
    {
      const size_t minor_length = minor_size.first + minor_size.second;
      const size_t ext_length = ext_height + ext_width;

      for (size_t minor_enum_iter = 0; minor_enum_iter <= minor_length; ++minor_enum_iter)
      {
        /// Type and index of enumerated element in smaller minor
        bool minor_enum_is_row = minor_enum_iter < minor_size.first;
        bool minor_enum_is_column = !minor_enum_is_row && minor_enum_iter != minor_length;
        size_t minor_enum_index = minor_enum_is_column ? minor_enum_iter - minor_size.first : minor_enum_iter;

        for (size_t bits = 1; bits < (size_t) (1 << ext_length); ++bits)
        {
          /// Set all mappings to 1 at first
          std::fill(row_mapping.begin(), row_mapping.begin() + minor_size.first + ext_height, 1);
          std::fill(row_mapping.begin() + minor_size.first + ext_height, row_mapping.end(), 0);
          std::fill(column_mapping.begin(), column_mapping.begin() + minor_size.second + ext_width, 1);
          std::fill(column_mapping.begin() + minor_size.second + ext_width, column_mapping.end(), 0);

          /// Set minor enumerated to -1
          size_t num_enumerated = 1;
          if (minor_enum_is_row)
            row_mapping[minor_enum_index] = -1;
          else if (minor_enum_is_column)
            column_mapping[minor_enum_index] = -1;
          else
            num_enumerated = 0;

          /// Iterator over extension bits and set row/column mapping to -1 if the bit is set
          for (size_t ext_enum_iter = 0; ext_enum_iter < ext_length; ++ext_enum_iter)
          {
            if (((bits >> ext_enum_iter) & 1) == 1)
            {
              ++num_enumerated;
              if (ext_enum_iter < ext_height)
                row_mapping[minor_size.first + ext_enum_iter] = -1;
              else
                column_mapping[minor_size.second + ext_enum_iter - ext_height] = -1;
            }
          }

          if (num_enumerated < 2)
            continue;

          /// Transform the mappings into row/column permutations.
          size_pair_t heights = detail::apply_mapping(worker_matrix.perm1(), row_mapping);
          size_pair_t widths = detail::apply_mapping(worker_matrix.perm2(), column_mapping);

          ++enumeration;

          if (log.is_updating() && enumeration == next_enumeration)
          {
            log.erase(cut);
            log.line() << next_percent << "%";
            std::cout << log;
            next_percent++;
            next_enumeration = (max_enumerations * next_percent) / 100;
          }

          if (detail::extend_to_3_4_separation(matroid, matrix, worker_matrix, size_pair_t(heights.first, widths.first), size_pair_t(heights.second,
              widths.second), separation, extra_elements))
          {
            return true;
          }
        }
      }

      return false;
    }

  }

  /**
   * Enumerates the partitions along the sequence of nested minors.
   *
   * @param matroid Given matroid
   * @param matrix Representation matrix of the given matroid
   * @param nested_minors Sequence of nested minors
   * @param extra_elements Set of matroid-elements to be filled with pivot-elements
   * @param log Logger
   * @return true if and only if the partition algorithm was successful
   */

  template <typename MatroidType, typename MatrixType, typename NestedMinorSequence>
  inline separation enumerate_separations (MatroidType& matroid, MatrixType& matrix, const NestedMinorSequence& nested_minors,
      matroid_element_set& extra_elements, logger& log)
  {
    typedef signed char mapping_value_t;

    /// Every regular 3-connected matroid which is non-graphic, non-cographic and not isomorphic to R10 must contain R12

    if (matrix.size1() + matrix.size2() < 12)
    {
      if (log.is_updating())
      {
        log.line() << ", TOO SMALL --> IRREGULAR";
        std::cout << log << std::endl;
        log.clear();
      }
      else if (log.is_verbose())
      {
        std::cout << "Matroid is too small to contain a (3|4)-separation and must be irregular." << std::endl;
      }

      return separation();
    }

    /// The first minors must be enumerated fully until we have one with at least 8 elements.
    size_pair_t minor_size;
    size_t minor_index = detail::find_first_8_element_minor(nested_minors, minor_size);
    assert (minor_size.first + minor_size.second >= 8);
    assert (minor_size.first + minor_size.second <= 10);

    /// preparations
    const integer_matrix worker_matrix_base(matrix);

    matrix_permuted <const integer_matrix> worker_matrix(worker_matrix_base);
    std::vector <mapping_value_t> row_mapping(worker_matrix_base.size1(), 0);
    std::vector <mapping_value_t> column_mapping(worker_matrix_base.size2(), 0);
    separation result;

    unsigned long long enumeration = 0;

    /// Calculate number of enumerations
    size_t cut = 0, full_cut = 0;
    unsigned long long max_enumerations = 0;

    if (log.is_updating() || log.is_verbose())
    {
      max_enumerations = 1L << (minor_size.first + minor_size.second);
      size_t h = minor_size.first;
      size_t w = minor_size.second;
      for (size_t i = minor_index; i < nested_minors.size(); ++i)
      {
        switch (nested_minors.get_extension(i))
        {
        case nested_minor_sequence::ONE_COLUMN:
        case nested_minor_sequence::ONE_ROW:
          max_enumerations += h + w;
        break;
        case nested_minor_sequence::ONE_ROW_ONE_COLUMN:
          max_enumerations += 3 * (h + w) + 1;
        break;
        case nested_minor_sequence::ONE_ROW_TWO_COLUMNS:
        case nested_minor_sequence::TWO_ROWS_ONE_COLUMN:
          max_enumerations += 7 * (h + w) + 4;
        break;
        default:
          assert (false);
        }
        h += nested_minors.get_extension_height(i);
        w += nested_minors.get_extension_width(i);
      }

      if (log.is_updating())
      {
        full_cut = log.size();
        log.line() << ", ENUMERATING " << max_enumerations << " PARTITIONS: ";
        cut = log.size();
        log.line() << "0%";
        std::cout << log;
      }
      else if (log.is_verbose())
      {
        std::cout << "Enumerating " << max_enumerations << " partitions along the sequence." << std::endl;
      }
    }

    unsigned long long next_enumeration = max_enumerations / 100;
    unsigned int next_percent = 1;

    /// Full enumeration
    size_t limit = (1 << (minor_size.first + minor_size.second));
    for (size_t bits = 0; bits < limit; ++bits)
    {
      ++enumeration;

      for (size_t row = 0; row < minor_size.first; ++row)
      {
        row_mapping[row] = ((bits >> row) & 0x1) ? 1 : -1;
      }
      for (size_t column = 0; column < minor_size.second; ++column)
      {
        column_mapping[column] = ((bits >> (minor_size.first + column)) & 0x1) ? 1 : -1;
      }

      /// Transform the mappings into row/column permutations.

      size_pair_t widths = detail::apply_mapping(worker_matrix.perm2(), column_mapping);
      size_pair_t heights = detail::apply_mapping(worker_matrix.perm1(), row_mapping);

      if (detail::extend_to_3_4_separation(matroid, matrix, worker_matrix, size_pair_t(heights.first, widths.first), size_pair_t(heights.second,
          widths.second), result, extra_elements))
      {
        return result;
      }

      if (log.is_updating() && enumeration == next_enumeration)
      {
        log.erase(cut);
        log.line() << next_percent << '%';
        std::cout << log;
        next_percent++;
        next_enumeration = (max_enumerations * next_percent) / 100;
      }
    }

    if (log.is_verbose())
    {
      std::cout << "Complete enumeration of " << (minor_index + 1) << " minors done." << std::endl;
    }

    /// Clever enumeration along the sequence
    for (size_t i = minor_index; i < nested_minors.size(); ++i)
    {
      size_t extension_height = nested_minors.get_extension_height(i);
      size_t extension_width = nested_minors.get_extension_width(i);
      if (detail::enumerate_extension(matroid, matrix, worker_matrix, row_mapping, column_mapping, minor_size, extension_height, extension_width,
          result, extra_elements, enumeration, next_enumeration, max_enumerations, next_percent, log, cut))
      {
        return result;
      }
      minor_size.first += extension_height;
      minor_size.second += extension_width;

      if (log.is_verbose())
      {
        std::cout << "Clever enumeration done for nested minor " << (i + 2) << "." << std::endl;
      }
    }

    if (log.is_updating())
    {
      log.erase(full_cut);
      log.line() << ", ENUMERATED " << enumeration << " PARTITIONS --> IRREGULAR";
      std::cout << log << std::endl;
      log.clear();
    }

    if (log.is_verbose())
    {
      std::cout << "Matroid does not contain a (3|4)-separation and must be irregular." << std::endl;
    }

    return separation();
  }

}

#endif /* ENUMERATION_HPP_ */
