#include <gtest/gtest.h>

#include <stdlib.h>

#include "common.h"

#include <cmr/network.h>

void testNetworkMatrix(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix  /**< Matrix to be used for testing. */
)
{
  bool isGraphic;
  CMR_GRAPH* graph = NULL;
  CMR_GRAPH_EDGE* basis = NULL;
  CMR_GRAPH_EDGE* cobasis = NULL;
  bool* edgesReversed = NULL;

  ASSERT_CMR_CALL( CMRtestNetworkMatrix(cmr, matrix, &isGraphic, &graph, &basis, &cobasis, &edgesReversed, NULL,
    NULL, DBL_MAX) );

  ASSERT_TRUE( isGraphic );
  ASSERT_TRUE( basis );
  ASSERT_TRUE( cobasis );

  CMR_CHRMAT* result = NULL;
  bool isCorrectBasis;
  ASSERT_CMR_CALL( CMRcomputeNetworkMatrix(cmr, graph, &result, NULL, edgesReversed, matrix->numRows,
    basis, matrix->numColumns, cobasis, &isCorrectBasis) );
  ASSERT_TRUE( isCorrectBasis );
  ASSERT_TRUE( result );

  if (!CMRchrmatCheckEqual(matrix, result))
  {
    printf("Input matrix:\n");
    ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', true) );
  
    printf("Graph:\n");
    ASSERT_CMR_CALL( CMRgraphPrint(stdout, graph) );

    printf("Representation matrix:\n");
    ASSERT_CMR_CALL( CMRchrmatPrintDense(cmr, result, stdout, '0', true) );

    printf("Basis:");
    for (int r = 0; r < matrix->numRows; ++r)
      printf(" %d", basis[r]);
    printf("\n");

    printf("Cobasis:");
    for (int c = 0; c < matrix->numColumns; ++c)
      printf(" %d", cobasis[c]);
    printf("\n");
  }

  ASSERT_TRUE( CMRchrmatCheckEqual(matrix, result) );
  ASSERT_TRUE( isGraphic );

  ASSERT_CMR_CALL( CMRgraphFree(cmr, &graph) );
  ASSERT_CMR_CALL( CMRfreeBlockArray(cmr, &basis) );
  ASSERT_CMR_CALL( CMRfreeBlockArray(cmr, &cobasis) );
  ASSERT_CMR_CALL( CMRfreeBlockArray(cmr, &edgesReversed) );
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &result) );
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

  testNetworkMatrix(cmr, matrix);
  
  ASSERT_CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  ASSERT_CMR_CALL( CMRfreeEnvironment(&cmr) );
}
