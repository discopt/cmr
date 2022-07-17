#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <cmr/matrix.h>
#include <cmr/graphic.h>
#include <cmr/graph.h>
#include <cmr/network.h>

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,       /**< Whether the file format of input/output was defined by the user. */
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
  FILEFORMAT_GRAPH_EDGELIST = 3,  /**< Edge list digraph format. */
  FILEFORMAT_GRAPH_DOT = 4,       /**< Dot digraph format. */
} FileFormat;

/**
 * \brief Prints the usage of the \p program to stdout.
 * 
 * \returns \c EXIT_FAILURE.
 */

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... FILE\n\n", program);
  puts("Converts digraph to (co)network matrix or tests if matrix is network, depending on input FILE.");
  puts("Options:");
  puts("  -i FORMAT  Format of input FILE; default: `dense'.");
  puts("  -o FORMAT  Format of output; default: `edgelist' if input is a matrix and `dense' if input is a digraph.");
  puts("  -t         Tests for being / converts to conetwork matrix.");
  puts("  -n         Output the elements of a minimal non-(co)network submatrix.");
  puts("  -N         Output a minimal non-(co)network submatrix.");
  puts("  -s         Print statistics about the computation to stderr.");
  puts("Formats for matrices: dense, sparse");
  puts("Formats for digraphs: edgelist, dot (output only)");
  puts("If FILE is `-', then the input will be read from stdin.");

  return EXIT_FAILURE;
}

/**
 * \brief Converts matrix from a file to a digraph if the former is (co)network.
 */

CMR_ERROR matrixToDigraph(
  const char* instanceFileName,   /**< File name containing the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,         /**< Format of the input matrix. */
  FileFormat outputFormat,        /**< Format of the output digraph. */
  bool conetwork,                 /**< Whether the input shall be checked for being conetwork instead of network. */
  bool outputNonNetworkElements,  /**< Whether to print the elements of a minimal non-(co)network submatrix. */
  bool outputNonNetworkMatrix,    /**< Whether to print a minimal non-(co)network submatrix. */
  bool printStats                 /**< Whether to print statistics to stderr. */
)
{
  clock_t readClock = clock();
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_CHRMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, instanceFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, instanceFile, &matrix) );
  if (instanceFile != stdin)
    fclose(instanceFile);
  fprintf(stderr, "Read %lux%lu matrix with %lu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Test for network. */

  bool isCoNetwork;
  CMR_GRAPH* digraph = NULL;
  CMR_GRAPH_EDGE* rowEdges = NULL;
  CMR_GRAPH_EDGE* columnEdges = NULL;
  bool* edgesReversed = NULL;

  CMR_NETWORK_STATISTICS stats;
  CMR_CALL( CMRstatsNetworkInit(&stats) );
  if (conetwork)
  {
    CMR_CALL( CMRtestConetworkMatrix(cmr, matrix, &isCoNetwork, &digraph, &rowEdges, &columnEdges, &edgesReversed,
      NULL, &stats) );
  }
  {
    CMR_CALL( CMRtestNetworkMatrix(cmr, matrix, &isCoNetwork, &digraph, &rowEdges, &columnEdges, &edgesReversed,
      NULL, &stats) );
  }

  fprintf(stderr, "Matrix %s%snetwork.\n", isCoNetwork ? "IS " : "IS NOT ", conetwork ? "co" : "");
  if (printStats)
    CMR_CALL( CMRstatsNetworkPrint(stderr, &stats, NULL) );

  if (isCoNetwork)
  {
    if (outputFormat == FILEFORMAT_GRAPH_EDGELIST)
    {
      if (conetwork)
      {
        for (size_t column = 0; column < matrix->numColumns; ++column)
        {
          CMR_GRAPH_EDGE e = columnEdges[column];
          CMR_GRAPH_NODE u = CMRgraphEdgeU(digraph, e);
          CMR_GRAPH_NODE v = CMRgraphEdgeV(digraph, e);
          if (edgesReversed && edgesReversed[e])
          {
            CMR_GRAPH_NODE temp = u;
            u = v;
            v = temp;
          }
          printf("%d %d c%ld\n", u, v, column+1);
        }
        for (size_t row = 0; row < matrix->numRows; ++row)
        {
          CMR_GRAPH_EDGE e = rowEdges[row];
          CMR_GRAPH_NODE u = CMRgraphEdgeU(digraph, e);
          CMR_GRAPH_NODE v = CMRgraphEdgeV(digraph, e);
          if (edgesReversed && edgesReversed[e])
          {
            CMR_GRAPH_NODE temp = u;
            u = v;
            v = temp;
          }
          printf("%d %d r%ld\n", u, v, row+1);
        }
      }
      else
      {
        for (size_t row = 0; row < matrix->numRows; ++row)
        {
          CMR_GRAPH_EDGE e = rowEdges[row];
          CMR_GRAPH_NODE u = CMRgraphEdgeU(digraph, e);
          CMR_GRAPH_NODE v = CMRgraphEdgeV(digraph, e);
          if (edgesReversed && edgesReversed[e])
          {
            CMR_GRAPH_NODE temp = u;
            u = v;
            v = temp;
          }
          printf("%d %d r%ld\n", u, v, row+1);
        }
        for (size_t column = 0; column < matrix->numColumns; ++column)
        {
          CMR_GRAPH_EDGE e = columnEdges[column];
          CMR_GRAPH_NODE u = CMRgraphEdgeU(digraph, e);
          CMR_GRAPH_NODE v = CMRgraphEdgeV(digraph, e);
          if (edgesReversed && edgesReversed[e])
          {
            CMR_GRAPH_NODE temp = u;
            u = v;
            v = temp;
          }
          printf("%d %d c%ld\n", u, v, column+1);
        }
      }
    }
    else if (outputFormat == FILEFORMAT_GRAPH_DOT)
    {
      char buffer[16];
      puts("digraph G {");
      for (size_t row = 0; row < matrix->numRows; ++row)
      {
        CMR_GRAPH_EDGE e = rowEdges[row];
        CMR_GRAPH_NODE u = CMRgraphEdgeU(digraph, e);
        CMR_GRAPH_NODE v = CMRgraphEdgeV(digraph, e);
        if (edgesReversed && edgesReversed[e])
        {
          CMR_GRAPH_NODE temp = u;
          u = v;
          v = temp;
        }
        printf(" v_%d -> v_%d [label=\"%s\",style=bold,color=red];\n", u, v,
          CMRelementString(CMRrowToElement(row), buffer));
      }
      for (size_t column = 0; column < matrix->numColumns; ++column)
      {
        CMR_GRAPH_EDGE e = columnEdges[column];
        CMR_GRAPH_NODE u = CMRgraphEdgeU(digraph, e);
        CMR_GRAPH_NODE v = CMRgraphEdgeV(digraph, e);
        if (edgesReversed && edgesReversed[e])
        {
          CMR_GRAPH_NODE temp = u;
          u = v;
          v = temp;
        }
        printf(" v_%d -> v_%d [label=\"%s\"];\n", u, v, CMRelementString(CMRcolumnToElement(column), buffer));
      }
      puts("}");
    }

    if (edgesReversed)
      CMR_CALL( CMRfreeBlockArray(cmr, &edgesReversed) );
    CMR_CALL( CMRfreeBlockArray(cmr, &rowEdges) );
    CMR_CALL( CMRfreeBlockArray(cmr, &columnEdges) );
    CMR_CALL( CMRgraphFree(cmr, &digraph) );
  }

  /* Cleanup. */

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

/**
 * \brief Converts the given digraph file to the corresponding (co)network matrix.
 */

CMR_ERROR digraphToMatrix(
  const char* instanceFileName, /**< File name containing the input digraph (may be `-' for stdin). */
  FileFormat inputFormat,       /**< Format of the input digraph. */
  FileFormat outputFormat,      /**< Format of the output matrix. */
  bool conetwork                /**< Whether the output shall be the conetwork matrix instead of the network matrix. */
)
{
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read edge list. */

  CMR_GRAPH* digraph = NULL;
  CMR_ELEMENT* edgeElements = NULL;
  if (inputFormat == FILEFORMAT_GRAPH_EDGELIST)
  {
    CMR_CALL( CMRgraphCreateFromEdgeList(cmr, &digraph, &edgeElements, NULL, instanceFile) );
  }
  if (instanceFile != stdin)
    fclose(instanceFile);

  /* Scan edges for (co)forest edges. */
  size_t numForestEdges = 0;
  size_t numCoforestEdges = 0;
  for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(digraph); CMRgraphEdgesValid(digraph, i); i = CMRgraphEdgesNext(digraph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(digraph, i);
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

  for (CMR_GRAPH_ITER i = CMRgraphEdgesFirst(digraph); CMRgraphEdgesValid(digraph, i); i = CMRgraphEdgesNext(digraph, i))
  {
    CMR_GRAPH_EDGE e = CMRgraphEdgesEdge(digraph, i);
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

  if (conetwork)
  {
    CMR_CALL( CMRcomputeNetworkMatrix(cmr, digraph, NULL, &matrix, NULL, numForestEdges, forestEdges, numCoforestEdges,
      coforestEdges, &isCorrectForest) );
  }
  else
  {
    CMR_CALL( CMRcomputeNetworkMatrix(cmr, digraph, &matrix, NULL, NULL, numForestEdges, forestEdges, numCoforestEdges,
      coforestEdges, &isCorrectForest) );
  }

  clock_t endTime = clock();
  fprintf(stderr, "Time: %f\n", (endTime - startTime) * 1.0 / CLOCKS_PER_SEC);

  if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );
  else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatPrintSparse(cmr, matrix, stdout) );
  else
    assert(false);

  CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  free(coforestEdges);
  free(forestEdges);
  free(edgeElements);
  CMR_CALL( CMRgraphFree(cmr, &digraph) );

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

int main(int argc, char** argv)
{
  FileFormat inputFormat = FILEFORMAT_UNDEFINED;
  FileFormat outputFormat = FILEFORMAT_UNDEFINED;
  bool transposed = false;
  char* instanceFileName = NULL;
  bool outputNonNetworkElements = false;
  bool outputNonNetworkMatrix = false;
  bool printStats = false;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-t"))
      transposed = true;
    else if (!strcmp(argv[a], "-n"))
      outputNonNetworkElements = true;
    else if (!strcmp(argv[a], "-N"))
      outputNonNetworkMatrix = true;
    else if (!strcmp(argv[a], "-s"))
      printStats = true;
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
      else if (!strcmp(argv[a+1], "dot"))
        outputFormat = FILEFORMAT_GRAPH_DOT;
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
  {
    error = matrixToDigraph(instanceFileName, inputFormat, outputFormat, transposed, outputNonNetworkElements,
      outputNonNetworkMatrix, printStats);
  }
  else
    error = digraphToMatrix(instanceFileName, inputFormat, outputFormat, transposed);

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
