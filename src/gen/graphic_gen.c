#include <cmr/graph.h>
#include <cmr/matrix.h>
#include <cmr/graphic.h>
#include <time.h>
#include <sys/time.h>

int printUsage(const char* program)
{
  fprintf(stderr, "Usage: %s [OPTIONS] #ROWS #COLUMNS\n\n", program);
  fputs("Creates a random #ROWS-by-#COLUMNS network matrix.\n", stderr);
  fputs("Options:\n", stderr);
  return EXIT_FAILURE;
}

int main(int argc, const char** argv)
{
  int numNodes, numEdges;
  struct timeval time; 
  gettimeofday(&time, NULL);
  srand((time.tv_sec * 1000) + (time.tv_usec / 1000));

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
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Init transpose of matrix. */
  CMR_CHRMAT* transposed = NULL;
  CMR_CALL( CMRchrmatCreate(cmr, &transposed, numEdges, numNodes-1, numEdges * (numNodes-1)) );
  transposed->numNonzeros = 0;

  /* Create random arborescence. */
  int* nextTreeNode = NULL;
  int* treeDistance = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &nextTreeNode, numNodes) );
  CMR_CALL( CMRallocBlockArray(cmr, &treeDistance, numNodes) );
  nextTreeNode[0] = 0;
  treeDistance[0] = 0;
  for (int v = 1; v < numNodes; ++v)
  {
    int w = (int)(rand() * 1.0 * v / RAND_MAX);
    nextTreeNode[v] = w;
    treeDistance[v] = treeDistance[w] + 1;
  }

  int* column = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &column, numNodes - 1) );
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
  transposed->rowSlice[transposed->numRows] = transposed->numNonzeros;

  CMR_CALL( CMRfreeBlockArray(cmr, &column) );
  CMR_CALL( CMRfreeBlockArray(cmr, &nextTreeNode) );
  CMR_CALL( CMRfreeBlockArray(cmr, &treeDistance) );

  CMR_CHRMAT* matrix = NULL;
  CMR_CALL( CMRchrmatTranspose(cmr, transposed, &matrix) );

  /* Print matrix. */

  CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );

  /* Check for graphicness. */

//   CMR_GRAPH* graph = NULL;
//   CMR_GRAPH_EDGE* basis = NULL;
//   CMR_GRAPH_EDGE* cobasis = NULL;
//   CMR_SUBMAT* submatrix = NULL;
//   bool isGraphic;
// 
//   CMR_CALL( CMRtestGraphicMatrix(cmr, matrix, &isGraphic, &graph, &basis, &cobasis, &submatrix) );
// 
//   if (graph)
//   {
//     fprintf(stderr, "Represented graph:\n");
//     CMRgraphPrint(stderr, graph);
//     if (basis)
//     {
//       for (int r = 0; r < matrix->numRows; ++r)
//         fprintf(stderr, "Row %d corresponds to edge %d.\n", r, basis[r]);
//     }
//     if (cobasis)
//     {
//       for (int c = 0; c < matrix->numColumns; ++c)
//         fprintf(stderr, "Col %d corresponds to edge %d.\n", c, cobasis[c]);
//     }
//   }

  /* Cleanup */

  CMRchrmatFree(cmr, &transposed);
  CMRchrmatFree(cmr, &matrix);
  CMRfreeEnvironment(&cmr);

  return EXIT_SUCCESS;
}
