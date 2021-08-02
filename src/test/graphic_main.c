#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <cmr/matrix.h>
#include <cmr/graphic.h>
#include <cmr/graph.h>

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

CMR_ERROR matrixToGraph(const char* instanceFileName, FileFormat inputFormat, FileFormat outputFormat, bool binary)
{
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_CHRMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, &matrix, instanceFile) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, &matrix, instanceFile) );
  if (instanceFile != stdin)
    fclose(instanceFile);

  /* Transpose it. */

  CMR_CHRMAT* transpose = NULL;
  CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );
  CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  /* Test for graphicness. */

  bool isGraphic;
  CMR_GRAPH* graph = NULL;
  CMR_GRAPH_EDGE* forestEdges = NULL;
  CMR_GRAPH_EDGE* coforestEdges = NULL;
  bool* edgesReversed = NULL;

  clock_t startTime = clock();

  if (binary)
    CMR_CALL( CMRtestBinaryGraphic(cmr, transpose, &isGraphic, &graph, &forestEdges, &coforestEdges, NULL) );
  else
    CMR_CALL( CMRtestTernaryGraphic(cmr, transpose, &isGraphic, &graph, &forestEdges, &coforestEdges, &edgesReversed, NULL) );

  clock_t endTime = clock();
  fprintf(stderr, "Time: %f\n", (endTime - startTime) * 1.0 / CLOCKS_PER_SEC);

  fprintf(stderr, "%s input matrix is %sgraphic.\n", binary ? "Binary" : "Ternary", isGraphic ? "" : "NOT ");

  if (isGraphic)
  {
    if (outputFormat == FILEFORMAT_GRAPH_EDGELIST)
    {
      for (size_t row = 0; row < transpose->numColumns; ++row)
      {
        CMR_GRAPH_EDGE e = forestEdges[row];
        CMR_GRAPH_NODE u = CMRgraphEdgeU(graph, e);
        CMR_GRAPH_NODE v = CMRgraphEdgeV(graph, e);
        if (edgesReversed && edgesReversed[e])
        {
          CMR_GRAPH_NODE temp = u;
          u = v;
          v = temp;
        }
        printf("%d %d r%ld\n", u, v, row+1);
      }
      for (size_t column = 0; column < transpose->numRows; ++column)
      {
        CMR_GRAPH_EDGE e = coforestEdges[column];
        CMR_GRAPH_NODE u = CMRgraphEdgeU(graph, e);
        CMR_GRAPH_NODE v = CMRgraphEdgeV(graph, e);
        if (edgesReversed && edgesReversed[e])
        {
          CMR_GRAPH_NODE temp = u;
          u = v;
          v = temp;
        }
        printf("%d %d c%ld\n", u, v, column+1);
      }
    }

    if (edgesReversed)
      CMR_CALL( CMRfreeBlockArray(cmr, &edgesReversed) );
    CMR_CALL( CMRfreeBlockArray(cmr, &forestEdges) );
    CMR_CALL( CMRfreeBlockArray(cmr, &coforestEdges) );
    CMR_CALL( CMRgraphFree(cmr, &graph) );
  }

  /* Cleanup. */

  CMR_CALL( CMRchrmatFree(cmr, &transpose) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

CMR_ERROR graphToMatrix(const char* instanceFileName, FileFormat inputFormat, FileFormat outputFormat, bool binary)
{
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read edge list. */

  CMR_GRAPH* graph = NULL;
  CMR_ELEMENT* edgeElements = NULL;
  if (inputFormat == FILEFORMAT_GRAPH_EDGELIST)
  {
    CMR_CALL( CMRgraphCreateFromEdgeList(cmr, &graph, &edgeElements, NULL, instanceFile) );
  }
  if (instanceFile != stdin)
    fclose(instanceFile);

  /* Scan edges for (co)forest edges. */
  size_t numForestEdges = 0;
  size_t numCoforestEdges = 0;
  for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, i); i = CMRgraphEdgesNext(graph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(graph, i);
    CMR_ELEMENT element = edgeElements[e];
    if (CMRelementIsRow(element))
      numForestEdges++;
    else if (CMRelementIsColumn(element))
      numCoforestEdges++;
  }

  /* Create list of (co)forest edges. */
  CMR_GRAPH_EDGE* forestEdges = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &forestEdges, numForestEdges) );
  for (size_t i = 0; i < numForestEdges; ++i)
    forestEdges[i] = -1;
  CMR_GRAPH_EDGE* coforestEdges = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &coforestEdges, numCoforestEdges) );
  for (size_t i = 0; i < numCoforestEdges; ++i)
    coforestEdges[i] = -1;

  for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(graph); CMRgraphEdgesValid(graph, i); i = CMRgraphEdgesNext(graph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(graph, i);
    CMR_ELEMENT element = edgeElements[e];
    if (CMRelementIsRow(element))
    {
      size_t rowIndex = CMRelementToRowIndex(element);
      if (rowIndex < numForestEdges)
        forestEdges[rowIndex] = e;
    }
    else if (CMRelementIsColumn(element))
    {
      size_t columnIndex = CMRelementToColumnIndex(element);
      if (columnIndex < numCoforestEdges)
        coforestEdges[columnIndex] = e;
    }
  }

  CMR_CHRMAT* matrix = NULL;
  bool isCorrectForest = false;

  clock_t startTime = clock();

  if (binary)
  {
    CMR_CALL( CMRcomputeGraphBinaryRepresentationMatrix(cmr, graph, &matrix, NULL, numForestEdges, forestEdges,
      numCoforestEdges, coforestEdges, &isCorrectForest) );
  }
  else
  {
    CMR_CALL( CMRcomputeGraphTernaryRepresentationMatrix(cmr, graph, &matrix, NULL, NULL, numForestEdges, forestEdges,
      numCoforestEdges, coforestEdges, &isCorrectForest) );
  }

  clock_t endTime = clock();
  fprintf(stderr, "Time: %f\n", (endTime - startTime) * 1.0 / CLOCKS_PER_SEC);

  if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatPrintDense(cmr, stdout, matrix, '0', false) );
  else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatPrintSparse(stdout, matrix) );
  else
    assert(false);

  CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  free(coforestEdges);
  free(forestEdges);
  free(edgeElements);
  CMRgraphFree(cmr, &graph);

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
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

  CMR_ERROR error;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE || inputFormat == FILEFORMAT_MATRIX_SPARSE)
    error = matrixToGraph(instanceFileName, inputFormat, outputFormat, binary);
  else
    error = graphToMatrix(instanceFileName, inputFormat, outputFormat, binary);

  switch (error)
  {
  case CMR_ERROR_INPUT:
    puts("Input error.");
    return EXIT_FAILURE;
  case CMR_ERROR_MEMORY:
    puts("Memory error.");
    return EXIT_FAILURE;
  default:
    return EXIT_SUCCESS;
  }
}
