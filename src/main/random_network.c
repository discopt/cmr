#include <tu/graph.h>
#include <tu/matrix.h>
#include <tu/graphic.h>

void printUsage(const char* program)
{
  printf("Usage: %s #ROWS #COLUMNS\n", program);
}

int main(int argc, const char** argv)
{
  int numNodes, numEdges;

  if (argc != 3)
  {
    puts("Invalid number of arguments.");
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }
  if (sscanf(argv[1], "%d", &numNodes) != 1)
  {
    printf("Error: Invalid first argument <%s>.\n", argv[1]);
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }
  if (sscanf(argv[2], "%d", &numEdges) != 1)
  {
    printf("Error: Invalid second argument <%s>.\n", argv[2]);
    printUsage(argv[0]);
    return EXIT_FAILURE;
  }
  ++numNodes;

  CMR* cmr = NULL;
  CMRcreateEnvironment(&cmr);

  /* Init transpose of matrix. */
  CMR_CHRMAT* transposed = NULL;
  CMRchrmatCreate(cmr, &transposed, numEdges, numNodes-1, numEdges * (numNodes-1));
  transposed->numNonzeros = 0;

  /* Create random arborescence. */
  int* nextTreeNode = NULL;
  int* treeDistance = NULL;
  CMRallocBlockArray(cmr, &nextTreeNode, numNodes);
  CMRallocBlockArray(cmr, &treeDistance, numNodes);
  nextTreeNode[0] = 0;
  treeDistance[0] = 0;
  for (int v = 1; v < numNodes; ++v)
  {
    int w = (int)(rand() * 1.0 * v / RAND_MAX);
    nextTreeNode[v] = w;
    treeDistance[v] = treeDistance[w] + 1;
  }

  int* column = NULL;
  CMRallocBlockArray(cmr, &column, numNodes - 1);
  for (int e = 0; e < numEdges; ++e)
  {
    for (int v = 1; v < numNodes; ++v)
      column[v-1] = 0;
    int first = (int)(rand() * 1.0 * numNodes / RAND_MAX);
    int second = (int)(rand() * 1.0 * numNodes / RAND_MAX);
    while (treeDistance[first] > treeDistance[second])
    {
      column[first-1] = 1;
      first = nextTreeNode[first];
    }
    while (treeDistance[second] > treeDistance[first])
    {
      column[second-1] = 1;
      second = nextTreeNode[second];
    }
    while (first != second && first)
    {
      column[first-1] = 1;
      first = nextTreeNode[first];
      column[second-1] = 1;
      second = nextTreeNode[second];
    }
    transposed->rowSlice[e] = transposed->numNonzeros;
    for (int v = 1; v < numNodes; ++v)
    {
      if (column[v-1])
      {
        transposed->entryColumns[transposed->numNonzeros] = v-1;
        transposed->entryValues[transposed->numNonzeros] = 1;
        transposed->numNonzeros++;
      }
    }
  }

  CMRfreeBlockArray(cmr, &column);
  CMRfreeBlockArray(cmr, &nextTreeNode);
  CMRfreeBlockArray(cmr, &treeDistance);

  CMR_CHRMAT* matrix = NULL;
  CMRchrmatTranspose(cmr, transposed, &matrix);

  /* Print matrix. */

  CMRchrmatPrintDense(stdout, matrix, ' ', true);

  /* Check for graphicness. */

  CMR_GRAPH* graph = NULL;
  CMR_GRAPH_EDGE* basis = NULL;
  CMR_GRAPH_EDGE* cobasis = NULL;
  CMR_SUBMAT* submatrix = NULL;
  bool isGraphic;

  CMR_CALL( CMRtestBinaryGraphic(cmr, transposed, &isGraphic, &graph, &basis, &cobasis, &submatrix) );

  if (graph)
  {
    printf("Represented graph:\n");
    CMRgraphPrint(stdout, graph);
    if (basis)
    {
      for (int r = 0; r < matrix->numRows; ++r)
        printf("Row %d corresponds to edge %d.\n", r, basis[r]);
    }
    if (cobasis)
    {
      for (int c = 0; c < matrix->numColumns; ++c)
        printf("Col %d corresponds to edge %d.\n", c, cobasis[c]);
    }
  }

  /* Cleanup */

  CMRchrmatFree(cmr, &transposed);
  CMRchrmatFree(cmr, &matrix);
  CMRfreeEnvironment(&cmr);

  return EXIT_SUCCESS;
}
