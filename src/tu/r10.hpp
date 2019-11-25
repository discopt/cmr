#pragma once

#include <boost/graph/isomorphism.hpp>

namespace tu
{

  /**
   * Singleton class to test a bipartite graph to be ismorphic to the graphs representing R10.
   */

  class bipartite_r10_graphs
  {
  public:
    typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> graph_t;

  private:

    /**
     * Constructs the R10 graph class.
     */

    bipartite_r10_graphs()
    {
      typedef std::pair <int, int> E;

      E g1_edges[] =
      { E(0, 5), E(0, 8), E(0, 9), E(1, 5), E(1, 6), E(1, 9), E(2, 6), E(2, 7), E(2, 9), E(3, 7), E(3, 8), E(3, 9), E(4, 5), E(4, 6), E(4, 7),
          E(4, 8), E(4, 9) };
      g1 = graph_t(&g1_edges[0], &g1_edges[0] + sizeof(g1_edges) / sizeof(E), 10);

      E g2_edges[] =
      { E(0, 5), E(0, 8), E(0, 9), E(1, 5), E(1, 6), E(1, 8), E(2, 6), E(2, 7), E(2, 9), E(3, 7), E(3, 8), E(3, 9), E(4, 5), E(4, 6), E(4, 7) };
      g2 = graph_t(&g2_edges[0], &g2_edges[0] + sizeof(g2_edges) / sizeof(E), 10);
    }

    /**
     * Singleton instance function.
     *
     * @return The unique instance
     */

    static bipartite_r10_graphs& instance()
    {
      static bipartite_r10_graphs* instance = NULL;
      if (instance == NULL)
        instance = new bipartite_r10_graphs();
      return *instance;
    }

  public:

    /**
     * Checks a given graph to be isomorphic to the ones stored.
     *
     * @param graph Given graph
     * @return true if and only if it is isomorphic to any of the two stored ones.
     */

    template <typename Graph>
    static bool is_r10_graph(const Graph& graph)
    {
      return boost::isomorphism(graph, instance().g1) || boost::isomorphism(graph, instance().g2);
    }

  private:
    graph_t g1, g2;
  };

  /**
   * Tests a given matrix to be matrix-isomorphic to one of the R10-representing matrices
   * by examining the corresponding bipartite graphs.
   *
   * @param matrix A given matrix
   * @return true if and only if this matrix is a represenation matrix for R10
   */

  template <typename MatrixType>
  inline bool is_r10(MatrixType matrix)
  {
    if (matrix.size1() != 5 || matrix.size2() != 5)
      return false;

    bipartite_r10_graphs::graph_t graph(10);

    for (size_t row = 0; row < 5; ++row)
    {
      for (size_t column = 0; column < 5; ++column)
      {
        if (matrix(row, column))
        {
          boost::add_edge(boost::vertex(row, graph), boost::vertex(5 + column, graph), graph);
        }
      }
    }

    return bipartite_r10_graphs::is_r10_graph(graph);
  }
} /* namespace tu */
