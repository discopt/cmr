#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <tu/matrix.h>
#include <tu/graphic.h>
#include <tu/graph.h>

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,
  FILEFORMAT_MATRIX_DENSE = 1,
  FILEFORMAT_MATRIX_SPARSE = 2,
  FILEFORMAT_GRAPH_EDGELIST = 3
} FileFormat;

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... FILE\n\n", program);
  puts("Converts graph to representation matrix or tests if matrix represents a graph, depending on input FILE.");
  puts("Options:");
  puts("  -i FORMAT  Format of input FILE; default: `dense'.");
  puts("  -o FORMAT  Format of output; default: `edgelist' if input is a matrix and `dense' if input is a graph.");
  puts("  -b         Consider binary representation matrices (default: ternary).");
  puts("Formats for matrices: dense, sparse");
  puts("Formats for graphs: edgelist");
  puts("If FILE is `-', then the input will be read from stdin.");
  return EXIT_FAILURE;
}

TU_ERROR matrixToGraph(const char* instanceFileName, FileFormat inputFormat, FileFormat outputFormat, bool binary)
{
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return TU_ERROR_INPUT;

  TU* tu = NULL;
  TU_CALL( TUcreateEnvironment(&tu) );

  /* Read matrix. */

  TU_CHRMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    TU_CALL( TUchrmatCreateFromDenseStream(tu, &matrix, instanceFile) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    TU_CALL( TUchrmatCreateFromSparseStream(tu, &matrix, instanceFile) );
  if (instanceFile != stdin)
    fclose(instanceFile);

  /* Transpose it. */

  TU_CHRMAT* transpose = NULL;
  TU_CALL( TUchrmatTranspose(tu, matrix, &transpose) );
  TU_CALL( TUchrmatFree(tu, &matrix) );

  /* Test for graphicness. */

  bool isGraphic;
  TU_GRAPH* graph = NULL;
  TU_GRAPH_EDGE* forestEdges = NULL;
  TU_GRAPH_EDGE* coforestEdges = NULL;
  bool* edgesReversed = NULL;
  if (binary)
    TU_CALL( TUtestBinaryGraphic(tu, transpose, &isGraphic, &graph, &forestEdges, &coforestEdges, NULL) );
  else
    TU_CALL( TUtestTernaryGraphic(tu, transpose, &isGraphic, &graph, &forestEdges, &coforestEdges, &edgesReversed, NULL) );

  fprintf(stderr, "%s input matrix is %sgraphic.\n", binary ? "Binary" : "Ternary", isGraphic ? "" : "NOT ");

  if (isGraphic)
  {
    if (outputFormat == FILEFORMAT_GRAPH_EDGELIST)
    {
      for (size_t row = 0; row < transpose->numColumns; ++row)
      {
        TU_GRAPH_EDGE e = forestEdges[row];
        TU_GRAPH_NODE u = TUgraphEdgeU(graph, e);
        TU_GRAPH_NODE v = TUgraphEdgeV(graph, e);
        if (edgesReversed && edgesReversed[e])
        {
          TU_GRAPH_NODE temp = u;
          u = v;
          v = temp;
        }
        printf("%d %d r%ld\n", u, v, row+1);
      }
      for (size_t column = 0; column < transpose->numRows; ++column)
      {
        TU_GRAPH_EDGE e = coforestEdges[column];
        TU_GRAPH_NODE u = TUgraphEdgeU(graph, e);
        TU_GRAPH_NODE v = TUgraphEdgeV(graph, e);
        if (edgesReversed && edgesReversed[e])
        {
          TU_GRAPH_NODE temp = u;
          u = v;
          v = temp;
        }
        printf("%d %d c%ld\n", u, v, column+1);
      }
    }

    if (edgesReversed)
      TU_CALL( TUfreeBlockArray(tu, &edgesReversed) );
    TU_CALL( TUfreeBlockArray(tu, &forestEdges) );
    TU_CALL( TUfreeBlockArray(tu, &coforestEdges) );
    TU_CALL( TUgraphFree(tu, &graph) );
  }

  /* Cleanup. */

  TU_CALL( TUchrmatFree(tu, &transpose) );
  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
}

TU_ERROR graphToMatrix(const char* instanceFileName, FileFormat inputFormat, FileFormat outputFormat, bool binary)
{
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return TU_ERROR_INPUT;

  TU* tu = NULL;
  TU_CALL( TUcreateEnvironment(&tu) );

  /* Read edge list. */

  TU_GRAPH* graph = NULL;
  Element* edgeElements = NULL;
  if (inputFormat == FILEFORMAT_GRAPH_EDGELIST)
  {
    TU_CALL( TUgraphCreateFromEdgeList(tu, &graph, &edgeElements, NULL, instanceFile) );
  }
  if (instanceFile != stdin)
    fclose(instanceFile);

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

  if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    TU_CALL( TUchrmatPrintDense(stdout, matrix, '0', false) );
  else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    TU_CALL( TUchrmatPrintSparse(stdout, matrix) );
  else
    assert(false);

  TU_CALL( TUchrmatFree(tu, &matrix) );

  free(coforestEdges);
  free(forestEdges);
  free(edgeElements);
  TUgraphFree(tu, &graph);

  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
}

int main(int argc, char** argv)
{
  FileFormat inputFormat = FILEFORMAT_UNDEFINED;
  FileFormat outputFormat = FILEFORMAT_UNDEFINED;
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
      if (!strcmp(argv[a+1], "dense"))
        inputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        inputFormat = FILEFORMAT_MATRIX_SPARSE;
      else if (!strcmp(argv[a+1], "edgelist"))
        inputFormat = FILEFORMAT_GRAPH_EDGELIST;
      else
      {
        printf("Error: unknown input file format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-o") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        outputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        outputFormat = FILEFORMAT_MATRIX_SPARSE;
      else if (!strcmp(argv[a+1], "edgelist"))
        outputFormat = FILEFORMAT_GRAPH_EDGELIST;
      else
      {
        printf("Error: unknown output format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!instanceFileName)
      instanceFileName = argv[a];
    else
    {
      printf("Error: Two input files <%s> and <%s> specified.\n\n", instanceFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (!instanceFileName)
  {
    puts("No input file specified.\n");
    return printUsage(argv[0]);
  }

  if (inputFormat == FILEFORMAT_UNDEFINED)
  {
    if (outputFormat == FILEFORMAT_UNDEFINED)
    {
      inputFormat = FILEFORMAT_MATRIX_DENSE;
      outputFormat = FILEFORMAT_GRAPH_EDGELIST;
    }
    else if (outputFormat == FILEFORMAT_MATRIX_DENSE || outputFormat == FILEFORMAT_MATRIX_SPARSE)
      inputFormat = FILEFORMAT_GRAPH_EDGELIST;
    else
      inputFormat = FILEFORMAT_MATRIX_DENSE;
  }
  else if (inputFormat == FILEFORMAT_MATRIX_DENSE || inputFormat == FILEFORMAT_MATRIX_SPARSE)
  {
    if (outputFormat == FILEFORMAT_UNDEFINED)
      outputFormat = FILEFORMAT_GRAPH_EDGELIST;
    else if (outputFormat == FILEFORMAT_MATRIX_DENSE || outputFormat == FILEFORMAT_MATRIX_SPARSE)
    {
      puts("Either input or output must be a graph.\n");
      return printUsage(argv[0]);
    }
  }
  else
  {
    /* Input format is graph format. */
    if (outputFormat == FILEFORMAT_UNDEFINED)
      outputFormat = FILEFORMAT_MATRIX_DENSE;
    else if (!(outputFormat == FILEFORMAT_MATRIX_DENSE || outputFormat == FILEFORMAT_MATRIX_SPARSE))
    {
      puts("Either input or output must be a matrix.\n");
      return printUsage(argv[0]);
    }
  }

  TU_ERROR error;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE || inputFormat == FILEFORMAT_MATRIX_SPARSE)
    error = matrixToGraph(instanceFileName, inputFormat, outputFormat, binary);
  else
    error = graphToMatrix(instanceFileName, inputFormat, outputFormat, binary);
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
