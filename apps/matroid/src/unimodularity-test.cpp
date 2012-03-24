#include <iostream>

#include <polymake/client.h>

#include <polymake/Integer.h>
#include <polymake/Matrix.h>
#include <polymake/Set.h>
#include <polymake/Map.h>
#include <polymake/Graph.h>

#include "matroid_decomposition.hpp"
#include "total_unimodularity.hpp"
#include "unimodularity.hpp"

namespace polymake
{
  namespace matroid
  {
    perl::Object construct_graph(const unimod::matroid_graph* graph)
    {
      Graph<> g(boost::num_vertices(*graph));

      typedef boost::graph_traits<unimod::matroid_graph> traits;
      traits::edge_iterator edge_iter, edge_end;
      for (boost::tie(edge_iter, edge_end) = boost::edges(*graph); edge_iter != edge_end; ++edge_iter)
      {
        traits::edge_descriptor e = *edge_iter;
        traits::vertex_descriptor s,t;
        s = boost::source(e, *graph);
        t = boost::target(e, *graph);
        g.edge((int)s, (int)t);
      }

      //EdgeMap<Undirected,int> edge_elements(g);
      //for (typename Entire< Edges< Graph<Undirected> > >::const_iterator e_it=entire(edges(g)); !e_it.at_end(); ++e_it)
      //{
      //  edge_elements[*e_it] = 0;
      //  boost::edge(boost::vertex(e_it.from_node(), *graph), boost::vertex(e_it.to_node(), *graph));
      //}

      perl::Object result("graph::Graph<Undirected>");
      result.take("ADJACENCY") << g;

      return result;
    }

    perl::Object construct_decomposition_tree(unimod::decomposed_matroid* tree)
    {
      unimod::decomposed_matroid_leaf* leaf = dynamic_cast<unimod::decomposed_matroid_leaf*>(tree);
      unimod::decomposed_matroid_separator* separator = dynamic_cast<unimod::decomposed_matroid_separator*>(tree);

      perl::Object result(leaf != NULL ? "MatroidDecompositionLeaf" : "MatroidDecompositionSeparation");
      result.set_description() << "matroid decomposition of a " << (tree->is_regular() ? "non-regular" : "regular" ) << " matroid.";

      Set<int> elements;
      std::copy(tree->elements().begin(), tree->elements().end(), std::inserter(elements, elements.begin()));
      result.take("REGULAR") << tree->is_regular();    
      result.take("ELEMENTS") << elements;

      if (leaf != NULL)
      {
        result.take("R10") << leaf->is_R10();
        result.take("GRAPHIC") << (leaf->is_graphic() ? 1 : 0);
        if (leaf->is_graphic())
          result.take("GRAPH") << construct_graph(leaf->graph());
        result.take("COGRAPHIC") << (leaf->is_cographic() ? 1 : 0);
        if (leaf->is_cographic())
          result.take("COGRAPH") << construct_graph(leaf->cograph());
      }
      else
      {
        result.take("SEPARATION_TYPE") << separator->separation_type();
        result.take("FIRST_CHILD") << construct_decomposition_tree(separator->first());
        result.take("SECOND_CHILD") << construct_decomposition_tree(separator->second());
      }

      return result;
    }
  
    perl::ListReturn is_totally_unimodular(const Matrix<Integer>& matrix, perl::OptionSet options)
    {
      unimod::integer_matrix input_matrix(matrix.rows(), matrix.cols());
      for (size_t r = 0; r < input_matrix.size1(); ++r)
      {
        for (size_t c = 0; c < input_matrix.size2(); ++c)
        {
          input_matrix(r, c) = matrix(r, c).to_int();
        }
      }

      unimod::decomposed_matroid* tree = NULL;
      perl::ListReturn result;
      bool answer;
      if (options["violator"])
      {
        unimod::submatrix_indices violator;
        answer = unimod::is_totally_unimodular(input_matrix, tree, violator);
        result << answer;
        if (tree)
        {
          result << construct_decomposition_tree(tree);
          delete tree;
        }
        else
          result << 0;
        if (answer)
          result << 0;
        else
        {
          std::pair<Set<int>, Set<int> > indices;
          std::copy(violator.rows.begin(), violator.rows.end(), std::inserter(indices.first, indices.first.begin()));
          std::copy(violator.columns.begin(), violator.columns.end(), std::inserter(indices.second, indices.second.begin()));
          result << indices;
        }
      }
      else
      {
        answer = unimod::is_totally_unimodular(input_matrix, tree);
        result << answer;
        if (tree)
        {
          result << construct_decomposition_tree(tree);
          delete tree;
        }
        else
          result << 0;
        delete tree;
      }

      return result;
    }
    
    /// total unimodularity

    UserFunction4perl("# Tests a given //matrix// for total unimodularity."
      "# @param Matrix matrix"
      "# @return Bool",
      &is_totally_unimodular, "is_totally_unimodular(Matrix { violator => 1 })");

    bool is_unimodular_plain(const Matrix <Integer>& matrix)
    {
      size_t rank;
      unimod::integer_matrix input_matrix(matrix.rows(), matrix.cols());
      for (size_t r = 0; r < input_matrix.size1(); ++r)
      {
        for (size_t c = 0; c < input_matrix.size2(); ++c)
        {
          input_matrix(r, c) = matrix(r, c).to_int();
        }
      }

      return unimod::is_unimodular(input_matrix, rank);
    }

    bool is_strongly_unimodular_plain(const Matrix <Integer>& matrix)
    {
      size_t rank;
      unimod::integer_matrix input_matrix(matrix.rows(), matrix.cols());
      for (size_t r = 0; r < input_matrix.size1(); ++r)
      {
        for (size_t c = 0; c < input_matrix.size2(); ++c)
        {
          input_matrix(r, c) = matrix(r, c).to_int();
        }
      }

      return unimod::is_strongly_unimodular(input_matrix, rank);
    }

    perl::ListReturn is_k_modular_compute(const Matrix <Integer>& matrix)
    {
      size_t rank;
      unsigned int k;
      unimod::integer_matrix input_matrix(matrix.rows(), matrix.cols());
      for (size_t r = 0; r < input_matrix.size1(); ++r)
      {
        for (size_t c = 0; c < input_matrix.size2(); ++c)
        {
          input_matrix(r, c) = matrix(r, c).to_int();
        }
      }

      perl::ListReturn result;
      if (unimod::is_k_modular(input_matrix, rank, k))
        result << (int) k;
      return result;
    }

    perl::ListReturn is_strongly_k_modular_compute(const Matrix <Integer>& matrix)
    {
      size_t rank;
      unsigned int k;
      unimod::integer_matrix input_matrix(matrix.rows(), matrix.cols());
      for (size_t r = 0; r < input_matrix.size1(); ++r)
      {
        for (size_t c = 0; c < input_matrix.size2(); ++c)
        {
          input_matrix(r, c) = matrix(r, c).to_int();
        }
      }

      if (unimod::is_strongly_k_modular(input_matrix, rank, k))
      {

      }
      perl::ListReturn result;
      if (unimod::is_strongly_k_modular(input_matrix, rank, k))
        result << (int) k;
      return result;
    }

  /// unimodularity

  UserFunction4perl("# Tests a given //matrix// for unimodularity, without certificates."
      "# @param Matrix<Integer> matrix"
      "# @return Bool",
      &is_unimodular_plain,"is_unimodular(Matrix<Integer>)");

  /// strong unimodularity

  UserFunction4perl("# Tests a given //matrix// for strong unimodularity, without certificates."
      "# @param Matrix<Integer> matrix"
      "# @return Bool",
      &is_strongly_unimodular_plain,"is_strongly_unimodular(Matrix<Integer>)");

  /// k-modularity

  UserFunction4perl("# Tests a given //matrix// for k-modularity, without certificates. It also computes k."
      "# @param Matrix<Integer> matrix"
      "# @return ListReturn k or undefined",
      &is_k_modular_compute,"is_k_modular(Matrix<Integer>)");

  /// strong k-modularity

  UserFunction4perl("# Tests a given //matrix// for strong k-modularity, without certificates."
      "# @param Matrix<Integer> matrix"
      "# @return ListReturn k or undefined",
      &is_strongly_k_modular_compute,"is_strongly_k_modular(Matrix<Integer>)");
  
  /// Method to convert matroid elemtents to row or column indices.
  
  Function4perl(&elements_to_indices, "elements_to_indices");

  Set<int> elements_to_indices(const Set<int> set, int sign)
  {
    Set<int> result;
    for (Set<int>::const_iterator iter = set.begin(); iter != set.end(); ++iter)
    {
      int index = *iter * sign - 1;
      if (index >= 0)
        result.insert(index);
    }

    return result;
  }

  
}}

