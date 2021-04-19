#include <gtest/gtest.h>

#include <stdlib.h>

#include "common.h"
#include <tu/graphic.h>
#include <tu/tdec.h>

void testGraphicMatrix(
  TU* tu,             /**< \ref TU environment. */
  TU_CHRMAT* matrix,  /**< Matrix to be used for testing. */
  int mergeLeafBonds  /**< Leaf bonds of the t-decomposition are merged (1: at the end; 2: after each column). */
)
{
  TU_GRAPH* graph = NULL;
  ASSERT_TU_CALL( TUgraphCreateEmpty(tu, &graph, 0, 0) );
  TU_GRAPH_EDGE* basis = NULL;
  ASSERT_TU_CALL( TUallocBlockArray(tu, &basis, matrix->numRows) );
  TU_GRAPH_EDGE* cobasis = NULL;
  ASSERT_TU_CALL( TUallocBlockArray(tu, &cobasis, matrix->numColumns) );
  bool isGraphic;
  TU_CHRMAT* transpose = NULL;
  ASSERT_TU_CALL( TUchrmatTranspose(tu, matrix, &transpose) );

  ASSERT_TU_CALL( testGraphicnessTDecomposition(tu, matrix, transpose, &isGraphic, graph, basis,
    cobasis, NULL, mergeLeafBonds) );

  ASSERT_TRUE( isGraphic );
  ASSERT_TRUE( basis );
  ASSERT_TRUE( cobasis );

  ASSERT_TU_CALL( TUchrmatFree(tu, &transpose) );

  TU_CHRMAT* result = NULL;
  ASSERT_TU_CALL( TUconvertGraphToBinaryMatrix(tu, graph, &result, matrix->numRows, basis,
    matrix->numColumns, cobasis) );

  ASSERT_TRUE( result );

  if (TUchrmatCheckEqual(matrix, result))
  {
    printf("The representation matrix of represented graph is equal to input matrix.\n");
  }
  else
  {
    printf("Input matrix:\n");
    ASSERT_TU_CALL( TUchrmatPrintDense(stdout, matrix, ' ', true) );
  
    printf("Graph:\n");
    ASSERT_TU_CALL( TUgraphPrint(stdout, graph) );

    printf("Representation matrix:\n");
    ASSERT_TU_CALL( TUchrmatPrintDense(stdout, result, ' ', true) );
  }

  ASSERT_TU_CALL( TUgraphFree(tu, &graph) );
  ASSERT_TU_CALL( TUfreeBlockArray(tu, &basis) );
  ASSERT_TU_CALL( TUfreeBlockArray(tu, &cobasis) );
  ASSERT_TU_CALL( TUchrmatFree(tu, &result) );
}

void testMatrix(
  TU* tu,             /**< \ref TU environment. */
  TU_CHRMAT* matrix,  /**< Matrix to be used for testing. */
  int mergeLeafBonds  /**< Leaf bonds of the t-decomposition are merged (1: at the end; 2: after each column). */
)
{
  TU_GRAPH* graph = NULL;
  ASSERT_TU_CALL( TUgraphCreateEmpty(tu, &graph, 0, 0) );
  bool isGraphic;
  TU_CHRMAT* transpose = NULL;
  ASSERT_TU_CALL( TUchrmatTranspose(tu, matrix, &transpose) );

  ASSERT_TU_CALL( testGraphicnessTDecomposition(tu, matrix, transpose, &isGraphic, graph, NULL,
    NULL, NULL, mergeLeafBonds) );

  ASSERT_TU_CALL( TUchrmatFree(tu, &transpose) );

  ASSERT_TU_CALL( TUgraphFree(tu, &graph) );
}

TEST(Graphic, RootBond)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  /* Just a bond. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "1 1 "
      "1 "
    ) );
    testGraphicMatrix(tu, A, 0);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* A larger bond */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "1 3 "
      "1 1 1 "
    ) );
    testGraphicMatrix(tu, A, 0);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }
  
  /* A root bond (attached to a polygon) with two child markers, each containing a (single-edge-) path. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "3 3 "
      "1 1  1 "
      "1 0  1 "
      "0 1  1 "
    ) );
    testGraphicMatrix(tu, A, 0);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* A bond with an additional edge and two child markers, each containing a (single-edge-) path. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "3 4 "
      "1 1 1  1 "
      "0 1 0  1 "
      "0 0 1  1 "
    ) );
    testGraphicMatrix(tu, A, 0);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* A root bond with a K_4 prime as a child, whose tree nodes are parallel edges. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "3 5 "
      "1 1 1 1  1 "
      "0 0 1 1  1 "
      "0 1 0 1  1 "
    ) );
    testGraphicMatrix(tu, A, 0);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, RootPolygon)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  /* Just a polygon */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "4 1 "
      "1 "
      "1 "
      "1 "
      "1 "
    ) );
    testGraphicMatrix(tu, A, 0);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* Polygon with a path edge and two attached bond/polygons containing the path ends. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "6 4 "
      "1 0 0  0 " /* edge of polygon */
      "1 0 0  1 " /* second edge of polygon */
      "1 1 0  0 " /* edge of first connecting bond. */
      "0 1 0  1 " /* edge of first triangle */
      "1 0 1  0 " /* edge of second connecting bond */
      "0 0 1  1 " /* edge of second triangle */
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* Polygon with no path edge and two attached bond/polygons containing the path ends. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "6 4 "
      "1 0 0  0 " /* edge of polygon */
      "1 0 0  0 " /* second edge of polygon */
      "1 1 0  0 " /* edge of first connecting bond. */
      "0 1 0  1 " /* edge of first triangle */
      "1 0 1  0 " /* edge of second connecting bond */
      "0 0 1  1 " /* edge of second triangle */
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }
  
  /* Polygon with a K_4. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "4 4 "
      "1 1 1 1 " /* first edge of polygon */
      "1 1 1 1 " /* second edge of polygon, where K_4 is attached. */
      "0 1 0 1 " /* */
      "0 0 1 1 " /* */
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, RootPrime)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  /* A K_4 with an attached polygon (via bond) containing an end node. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "5 5 "
      "1 1 0 0  1 "
      "1 0 1 1  0 "
      "0 1 1 1  1 "
      "0 0 0 1  1 "
      "0 0 0 1  1 "
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* A K_4 with two attached polygons (via bonds), each containing an end node. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "7 6 "
      "1 1 0 0 1  1 "
      "1 0 1 1 0  0 "
      "0 1 1 1 1  1 "
      "0 0 0 1 0  1 "
      "0 0 0 1 0  1 "
      "0 0 0 0 1  1 "
      "0 0 0 0 1  1 "
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* A type-4 prime with two attached polygons (via bonds), each containing an end node. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "8 9 "
      "1 1 1 1 1 1  1  1  1 "
      "1 1 1 1 1 1  1  1  1 "
      "1 1 1 1 1 1  1  1  1 "
      "0 1 0 1 1 0  1  1  1 "
      "0 0 1 1 1 1  1  1  1 "
      "0 0 0 0 1 1  1  0  0 "
      "0 0 0 0 0 0  1  0  1 "
      "0 0 0 0 0 0  0  1  1 "
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, InternalBond)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  /* A triangle linked to a small bond that is linked to a triangle.
   * Path contains one edge of each member. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "3 3 "
      "1 0  1 " /* edge of first triangle */
      "1 1  1 " /* bond edge */
      "0 1  1 " /* edge of second triangle */
    ) );
    testGraphicMatrix(tu, A, 0);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* A triangle linked to a large bond that is linked to a triangle.
   * Path contains one edge of each member. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "3 4 "
      "1 0 0  1 " /* edge of first triangle */
      "1 1 1  1 " /* bond edge */
      "0 0 1  1 " /* edge of second triangle */
    ) );
    testGraphicMatrix(tu, A, 0);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, InternalPrime)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  /* A triangle linked to a prime linked to a triangle with prime path of length 1. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "5 5 "
      "1 0 0  0  1 " /* edge of triangle root */
      "1 1 0  0  1 " /* edge of prime */
      "1 1 1  1  0 " /* edge of prime */
      "1 0 1  1  0 " /* edge of prime */
      "0 0 0  1  1 " /* edge of triangle leaf */
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* A triangle linked to a prime linked to a triangle with prime path of length 0. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "5 5 "
      "1 0 0  0  1 " /* edge of triangle root */
      "1 1 0  1  0 " /* edge of prime */
      "1 1 1  1  0 " /* edge of prime */
      "1 0 1  0  0 " /* edge of prime */
      "0 0 0  1  1 " /* edge of triangle leaf */
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* A triangle linked to a prime linked to a triangle with prime path of length 2 that ends at the second parent
   * marker node = child marker node. The path edges of the prime member are 3-star. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "5 6 "
      "1  0 0 0  0  1 " /* edge of triangle root */
      "1  1 1 0  0  1 " /* edge of prime */
      "1  1 0 1  1  1 " /* edge of prime */
      "0  0 1 1  1  0 " /* edge of prime */
      "0  0 0 0  1  1 " /* edge of triangle leaf */
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* A triangle linked to a prime linked to a triangle with prime path of length 3 that passes by at the second parent
   * marker node. The path edges of the prime member are 3-star. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "5 6 "
      "1  0 0 0  0  1 " /* edge of triangle root */
      "1  1 1 0  0  1 " /* edge of prime */
      "1  1 1 1  1  1 " /* edge of prime */
      "0  0 1 1  1  1 " /* edge of prime */
      "0  0 0 0  1  1 " /* edge of triangle leaf */
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* A triangle linked to a prime that is linked to two other triangles, next to the parent marker edge. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "6 7 "
      "1 1 1 1  0  0  1 " /* first edge of triangle root */
      "1 1 1 1  0  0  1 " /* second edge of triangle root */
      "0 1 0 1  1  0  0 " /* edge of prime */
      "0 0 1 1  0  1  0 " /* edge of prime */
      "0 0 0 0  1  0  1 " /* edge of triangle */
      "0 0 0 0  0  1  1 " /* edge of triangle */
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* A triangle linked to a prime that is linked to two other triangles, one next to the parent marker edge and one
   * via a path of length 1. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "6 7 "
      "1 1 1 1  1  0  1 " /* first edge of triangle root */
      "1 1 1 1  1  0  1 " /* second edge of triangle root */
      "0 1 0 1  1  0  1 " /* edge of prime */
      "0 0 1 1  0  1  0 " /* edge of prime */
      "0 0 0 0  1  0  1 " /* edge of triangle */
      "0 0 0 0  0  1  1 " /* edge of triangle */
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, LeafPrime)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  /* A triangle linked to a prime */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "4 4 "
      "1 0 0  1 " /* edge of triangle root */
      "1 1 0  1 " /* edge of prime */
      "1 1 1  1 " /* edge of prime */
      "1 0 1  0 " /* edge of prime */
    ) );
    testGraphicMatrix(tu, A, 2);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, BixbyWagnerAppendix)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "9 7 "
      "1 0 0 0 0 0 1 " // 0
      "0 1 1 1 0 0 1 " // 1
      "0 1 1 1 0 0 1 " // 2
      "0 1 1 1 0 1 1 " // 3
      "0 0 0 0 1 1 1 " // 4
      "0 1 1 1 1 0 0 " // 5
      "1 1 1 1 0 0 0 " // 6
      "1 1 0 1 0 0 0 " // 7
      "1 0 0 0 0 0 0 " // 8
    ) );
    testGraphicMatrix(tu, A, 0);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }
  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, Specials)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  /* Bond matrix plus empty row. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "2 1 "
      "1 "
      "0 "
    ) );
    testGraphicMatrix(tu, A, 0);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  /* Bond matrix plus empty column. */
  {
    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( stringToCharMatrix(tu, &A, "1 2 "
      "1 0 "
    ) );
    testGraphicMatrix(tu, A, 0);
    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, RandomMatrix)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  
  srand(0);
  const int numMatrices = 1000;
  const int numRows = 30;
  const int numColumns = 100;
  const double probability = 0.3;

  for (int i = 0; i < numMatrices; ++i)
  {
    TU_CHRMAT* A = NULL;
    TUchrmatCreate(tu, &A, numRows, numColumns, numRows * numColumns);

    A->numNonzeros = 0;
    for (int row = 0; row < numRows; ++row)
    {
      A->rowStarts[row] = A->numNonzeros;
      for (int column = 0; column < numColumns; ++column)
      {
        if ((rand() * 1.0 / RAND_MAX) < probability)
        {
          A->entryColumns[A->numNonzeros] = column;
          A->entryValues[A->numNonzeros] = 1;
          A->numNonzeros++;
        }
      }
    }
    A->rowStarts[numRows] = A->numNonzeros;

    testMatrix(tu, A, 0);

    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}

TEST(Graphic, RandomGraph)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );
  
  srand(0);
  const int numGraphs = 1000;
  const int numNodes = 5;
  const int numEdges = 30;

  TU_GRAPH_NODE* nodes = NULL;
  ASSERT_TU_CALL( TUallocBlockArray(tu, &nodes, numNodes) );
  for (int i = 0; i < numGraphs; ++i)
  {
    TU_GRAPH* graph = NULL;
    ASSERT_TU_CALL( TUgraphCreateEmpty(tu, &graph, numNodes, numEdges) );
    for (int v = 0; v < numNodes; ++v)
      ASSERT_TU_CALL( TUgraphAddNode(tu, graph, &nodes[v]) );

    for (int e = 0; e < numEdges; ++e)
    {
      int u = (rand() * 1.0 / RAND_MAX) * numNodes;
      int v = (rand() * 1.0 / RAND_MAX) * numNodes;
      ASSERT_TU_CALL( TUgraphAddEdge(tu, graph, nodes[u], nodes[v], NULL) );
    }

    TU_CHRMAT* A = NULL;
    ASSERT_TU_CALL( TUconvertGraphToBinaryMatrix(tu, graph, &A, 0, NULL, 0, NULL) );

    ASSERT_TU_CALL( TUgraphFree(tu, &graph) );

    testGraphicMatrix(tu, A, 0);

    ASSERT_TU_CALL( TUchrmatFree(tu, &A) );
  }
  ASSERT_TU_CALL( TUfreeBlockArray(tu, &nodes) );

  ASSERT_TU_CALL( TUfreeEnvironment(&tu) );
}
