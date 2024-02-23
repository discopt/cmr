#include <gtest/gtest.h>

#include <stdlib.h>

#include "common.h"

#include <cmr/network.h>

void testNetworkMatrix(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix, /**< Matrix to be used for testing. */
  bool knownNetwork   /**< Whether \p matrix is known to be network. */
)
{
  bool isNetwork;
  CMR_GRAPH* graph = NULL;
  CMR_GRAPH_EDGE* basis = NULL;
  CMR_GRAPH_EDGE* cobasis = NULL;
  bool* edgesReversed = NULL;
  CMR_SUBMAT* violatorSubmatrix = NULL;

  ASSERT_CMR_CALL( CMRnetworkTestMatrix(cmr, matrix, &isNetwork, NULL, &graph, &basis, &cobasis, &edgesReversed,
    &violatorSubmatrix, NULL, DBL_MAX) );

  ASSERT_EQ( isNetwork, knownNetwork );
  if (isNetwork)
  {
    ASSERT_TRUE( basis );
    ASSERT_TRUE( cobasis );

    CMR_CHRMAT* result = NULL;
    bool isCorrectBasis;
    ASSERT_CMR_CALL( CMRnetworkComputeMatrix(cmr, graph, &result, NULL, edgesReversed, matrix->numRows,
      basis, matrix->numColumns, cobasis, &isCorrectBasis) );
    ASSERT_TRUE( isCorrectBasis );
    ASSERT_TRUE( result );

    if (!CMRchrmatCheckEqual(matrix, result))
    {
      printf("Input matrix:\n");
      ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );

      printf("Graph:\n");
      ASSERT_CMR_CALL( CMRgraphPrint(graph, stdout) );

      printf("Representation matrix:\n");
      ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, result, stdout, '0', true) );

      printf("Basis:");
      for (size_t r = 0; r < matrix->numRows; ++r)
        printf(" %d", basis[r]);
      printf("\n");

      printf("Cobasis:");
      for (size_t c = 0; c < matrix->numColumns; ++c)
        printf(" %d", cobasis[c]);
      printf("\n");
    }

    ASSERT_TRUE( CMRchrmatCheckEqual(matrix, result) );

    ASSERT_CMR_CALL( CMRgraphFree(cmr, &graph) );
    ASSERT_CMR_CALL( CMRfreeBlockArray(cmr, &basis) );
    ASSERT_CMR_CALL( CMRfreeBlockArray(cmr, &cobasis) );
    ASSERT_CMR_CALL( CMRfreeBlockArray(cmr, &edgesReversed) );
    ASSERT_CMR_CALL( CMRchrmatFree(cmr, &result) );
    ASSERT_EQ( violatorSubmatrix, (CMR_SUBMAT*) NULL );
  }
  else
  {
    ASSERT_EQ( basis, (CMR_GRAPH_EDGE*) NULL );
    ASSERT_EQ( cobasis, (CMR_GRAPH_EDGE*) NULL );
    ASSERT_EQ( edgesReversed, (bool*) NULL );

    ASSERT_CMR_CALL( CMRsubmatPrint(cmr, violatorSubmatrix, matrix->numRows, matrix->numColumns, stdout) );

    ASSERT_CMR_CALL( CMRsubmatFree(cmr, &violatorSubmatrix) );
  }
}

TEST(Network, Basic)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 7 "
    "-1  0  0  0  1 -1  0 "
    " 1  0  0  1 -1  1  0 "
    " 0 -1  0 -1  1 -1  0 "
    " 0  1  0  0  0  0  1 "
    " 0  0  1 -1  1  0  1 "
    " 0  0 -1  1 -1  0  0 "
  ) );

  testNetworkMatrix(cmr, matrix, true);
  
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}

TEST(Network, NonCamion)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );
  CMR_CHRMAT* matrix = NULL;
  ASSERT_CMR_CALL( stringToCharMatrix(cmr, &matrix, "6 7 "
    "-1  0  0  0  1 -1  0 "
    " 1  0  0  1 -1  1  0 "
    " 0 -1  0 -1  1 -1  0 "
    " 0  1  0  0  0  0 -1 "
    " 0  0  1 -1  1  0  1 "
    " 0  0 -1  1 -1  0  0 "
  ) );

  testNetworkMatrix(cmr, matrix, false);

  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
