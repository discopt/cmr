#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <tu/matrix.h>
#include <tu/graphic.h>
#include <tu/graph.h>

typedef enum
{
  GRAPH_FORMAT_UNDEFINED = 0,
  GRAPH_FORMAT_EDGELIST = 1
} GraphFormat;

typedef enum
{
  MATRIX_FORMAT_DENSE = 1,
  MATRIX_FORMAT_SPARSE = 2
} MatrixFormat;

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... GRAPH\n\n", program);
  puts("Converts GRAPH to a representation matrix.");
  puts("\nOptions:");
  puts("  -i FORMAT  Format of input GRAPH file, among {edgelist}; default: edgelist.");
  puts("  -o FORMAT  Format of output matrix file, among {sparse,dense}; default: dense.");
  puts("  -b         Convert to binary representation matrix (default: ternary).");

  return EXIT_FAILURE;
}

TU_ERROR run(const char* instanceFileName, GraphFormat inputFormat, MatrixFormat outputFormat, bool binary)
{
  FILE* instanceFile = fopen(instanceFileName, "r");
  if (!instanceFile)
    return TU_ERROR_INPUT;

  TU* tu = NULL;
  TU_CALL( TUcreateEnvironment(&tu) );

  /* Read edge list. */

  TU_GRAPH* graph = NULL;
  Element* edgeElements = NULL;
  if (inputFormat == GRAPH_FORMAT_EDGELIST)
  {
    TU_CALL( TUgraphCreateFromEdgeList(tu, &graph, &edgeElements, NULL, instanceFile) );
  }
  fclose(instanceFile);

  if (graph)
  {
    /* Scan edges for (co)forest edges. */
    size_t numForestEdges = 0;
    size_t numCoforestEdges = 0;
    for (TU_GRAPH_ITER i = TUgraphEdgesFirst(graph); TUgraphEdgesValid(graph, i); i = TUgraphEdgesNext(graph, i))
    {
      TU_GRAPH_EDGE e = TUgraphEdgesEdge(graph, i);
      Element element = edgeElements[e];
      if (TUelementIsRow(element))
        numForestEdges++;
      else if (TUelementIsColumn(element))
        numCoforestEdges++;
    }

    /* Create list of (co)forest edges. */
    TU_GRAPH_EDGE* forestEdges = NULL;
    TU_CALL( TUallocBlockArray(tu, &forestEdges, numForestEdges) );
    for (size_t i = 0; i < numForestEdges; ++i)
      forestEdges[i] = -1;
    TU_GRAPH_EDGE* coforestEdges = NULL;
    TU_CALL( TUallocBlockArray(tu, &coforestEdges, numCoforestEdges) );
    for (size_t i = 0; i < numCoforestEdges; ++i)
      coforestEdges[i] = -1;

    for (TU_GRAPH_ITER i = TUgraphEdgesFirst(graph); TUgraphEdgesValid(graph, i); i = TUgraphEdgesNext(graph, i))
    {
      TU_GRAPH_EDGE e = TUgraphEdgesEdge(graph, i);
      Element element = edgeElements[e];
      if (TUelementIsRow(element))
      {
        size_t rowIndex = TUelementToRowIndex(element);
        if (rowIndex < numForestEdges)
          forestEdges[rowIndex] = e;
      }
      else if (TUelementIsColumn(element))
      {
        size_t columnIndex = TUelementToColumnIndex(element);
        if (columnIndex < numCoforestEdges)
          coforestEdges[columnIndex] = e;
      }
    }

    TU_CHRMAT* matrix = NULL;
    bool isCorrectForest = false;
    if (binary)
    {
      TU_CALL( TUcomputeGraphBinaryRepresentationMatrix(tu, graph, &matrix, NULL, numForestEdges, forestEdges,
        numCoforestEdges, coforestEdges, &isCorrectForest) );
    }
    else
    {
      TU_CALL( TUcomputeGraphTernaryRepresentationMatrix(tu, graph, &matrix, NULL, NULL, numForestEdges, forestEdges,
        numCoforestEdges, coforestEdges, &isCorrectForest) );
    }

    if (outputFormat == MATRIX_FORMAT_DENSE)
      TU_CALL( TUchrmatPrintDense(stdout, matrix, '0', false) );
    else if (outputFormat == MATRIX_FORMAT_SPARSE)
      TU_CALL( TUchrmatPrintSparse(stdout, matrix) );

    TU_CALL( TUchrmatFree(tu, &matrix) );

    free(coforestEdges);
    free(forestEdges);
    free(edgeElements);
    TUgraphFree(tu, &graph);
  }

  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
}

int main(int argc, char** argv)
{
  GraphFormat inputFormat = GRAPH_FORMAT_UNDEFINED;
  MatrixFormat outputFormat = MATRIX_FORMAT_DENSE;
  bool binary = false;
  char* instanceFileName = NULL;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-b"))
    {
      binary = true;
    }
    else if (!strcmp(argv[a], "-i") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "edgelist"))
        inputFormat = GRAPH_FORMAT_EDGELIST;
      else
      {
        printf("Error: unknown input graph format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-o") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        outputFormat = MATRIX_FORMAT_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        outputFormat = MATRIX_FORMAT_SPARSE;
      else
      {
        printf("Error: unknown output matrix format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!instanceFileName)
      instanceFileName = argv[a];
    else
    {
      printf("Error: Two input graph files <%s> and <%s> specified.\n\n", instanceFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (!instanceFileName)
  {
    puts("No input graph specified.\n");
    return printUsage(argv[0]);
  }

  TU_ERROR error = run(instanceFileName, inputFormat, outputFormat, binary);
  switch (error)
  {
  case TU_ERROR_INPUT:
    puts("Input error.");
    return EXIT_FAILURE;
  case TU_ERROR_MEMORY:
    puts("Memory error.");
    return EXIT_FAILURE;
  default:
    return EXIT_SUCCESS;
  }
}
