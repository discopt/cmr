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

  TU* tu = NULL;
  TUcreateEnvironment(&tu);

  /* Init transpose of matrix. */
  TU_CHRMAT* transposed = NULL;
  TUchrmatCreate(tu, &transposed, numEdges, numNodes-1, numEdges * (numNodes-1));
  transposed->numNonzeros = 0;

  /* Create random arborescence. */
  int* nextTreeNode = NULL;
  int* treeDistance = NULL;
  TUallocBlockArray(tu, &nextTreeNode, numNodes);
  TUallocBlockArray(tu, &treeDistance, numNodes);
  nextTreeNode[0] = 0;
  treeDistance[0] = 0;
  for (int v = 1; v < numNodes; ++v)
  {
    int w = (int)(rand() * 1.0 * v / RAND_MAX);
    nextTreeNode[v] = w;
    treeDistance[v] = treeDistance[w] + 1;
  }

  int* column = NULL;
  TUallocBlockArray(tu, &column, numNodes - 1);
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
    transposed->rowStarts[e] = transposed->numNonzeros;
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
  
  TUfreeBlockArray(tu, &column);
  TUfreeBlockArray(tu, &nextTreeNode);
  TUfreeBlockArray(tu, &treeDistance);

  TU_CHRMAT* matrix = NULL;
  TUchrmatTranspose(tu, transposed, &matrix);
  TUchrmatFree(tu, &transposed);

  /* Print matrix. */

  TUchrmatPrintDense(stdout, matrix, ' ', true);

  /* Check for graphicness. */

  TU_LISTGRAPH* graph = NULL;
  TU_LISTGRAPH_EDGE* basis = NULL;
  TU_LISTGRAPH_EDGE* cobasis = NULL;
  TU_SUBMAT* submatrix = NULL;

  TUtestGraphicnessChr(tu, matrix, &graph, &basis, &cobasis, &submatrix);

  if (graph)
  {
    printf("Represented graph:\n");
    TUlistgraphPrint(stdout, graph);
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

  TUchrmatFree(tu, &matrix);
  TUfreeEnvironment(&tu);

  return EXIT_SUCCESS;
}
