#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include <cmr/matrix.h>
#include <cmr/graphic.h>
#include <cmr/graph.h>
#include <cmr/network.h>

typedef enum
{
  TASK_RECOGNIZE = 0, /**< Determine whether a matrix is graphic. */
  TASK_COMPUTE = 1    /**< Compute a graphic matrix from a graph (and a tree). */
} Task;

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,     /**< Whether the file format of input/output was defined by the user. */
  FILEFORMAT_MATRIX_DENSE = 1,  /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2, /**< Sparse matrix format. */
} FileFormat;

/**
 * \brief Converts matrix from a file to a digraph if the former is (co)network.
 */

CMR_ERROR recognizeNetwork(
  const char* inputFileName,            /**< File name of the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,               /**< Format of the input matrix. */
  bool conetwork,                       /**< Whether the input shall be checked for being conetwork instead of network. */
  const char* outputGraphFileName,      /**< File name of the output graph (may be NULL; may be `-' for stdout). */
  const char* outputTreeFileName,       /**< File name of the output tree (may be NULL; may be `-' for stdout). */
  const char* outputDotFileName,        /**< File name of the output dot file (may be NULL; may be `-' for stdout). */
  const char* outputSubmatrixFileName,  /**< File name of the output non-(co)network submatrix (may be NULL; may be `-' for stdout). */
  bool printStats,                      /**< Whether to print statistics to stderr. */
  double timeLimit                  /**< Time limit to impose. */
)
{
  clock_t readClock = clock();
  FILE* inputFile = strcmp(inputFileName, "-") ? fopen(inputFileName, "r") : stdin;
  if (!inputFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_CHRMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatCreateFromDenseStream(cmr, inputFile, &matrix) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatCreateFromSparseStream(cmr, inputFile, &matrix) );
  if (inputFile != stdin)
    fclose(inputFile);
  fprintf(stderr, "Read %lux%lu matrix with %lu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Test for network. */

  bool isCoNetwork;
  CMR_GRAPH* digraph = NULL;
  CMR_GRAPH_EDGE* rowEdges = NULL;
  CMR_GRAPH_EDGE* columnEdges = NULL;
  bool* edgesReversed = NULL;
  CMR_SUBMAT* submatrix = NULL;

  CMR_NETWORK_STATISTICS stats;
  CMR_CALL( CMRstatsNetworkInit(&stats) );
  if (conetwork)
  {
    CMR_CALL( CMRtestConetworkMatrix(cmr, matrix, &isCoNetwork, &digraph, &rowEdges, &columnEdges, &edgesReversed,
      outputSubmatrixFileName ? &submatrix : NULL, &stats, timeLimit) );
  }
  {
    CMR_CALL( CMRtestNetworkMatrix(cmr, matrix, &isCoNetwork, &digraph, &rowEdges, &columnEdges, &edgesReversed,
      outputSubmatrixFileName ? &submatrix : NULL, &stats, timeLimit) );
  }

  fprintf(stderr, "Matrix %s%snetwork.\n", isCoNetwork ? "IS " : "is NOT ", conetwork ? "co" : "");
  if (printStats)
    CMR_CALL( CMRstatsNetworkPrint(stderr, &stats, NULL) );

  if (isCoNetwork)
  {
    if (outputGraphFileName)
    {
      bool outputGraphToFile = strcmp(outputGraphFileName, "-");
      FILE* outputGraphFile = outputGraphToFile ? fopen(outputGraphFileName, "w") : stdout;
      fprintf(stderr, "Writing digraph to %s%s%s.\n", outputGraphToFile ? "file <" : "",
        outputGraphToFile ? outputGraphFileName : "stdout", outputGraphToFile ? ">" : "");

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
          fprintf(outputGraphFile, "%d %d c%ld\n", u, v, column+1);
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
          fprintf(outputGraphFile, "%d %d r%ld\n", u, v, row+1);
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
          fprintf(outputGraphFile, "%d %d r%ld\n", u, v, row+1);
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
          fprintf(outputGraphFile, "%d %d c%ld\n", u, v, column+1);
        }
      }
      
      if (outputGraphToFile)
        fclose(outputGraphFile);
    }

    if (outputTreeFileName)
    {
      // TODO: implement
      assert(!"NOT IMPLEMENTED");
      exit(EXIT_FAILURE);
    }

    if (outputDotFileName)
    {
      bool outputDotToFile = strcmp(outputDotFileName, "-");
      FILE* outputDotFile = outputDotToFile ? fopen(outputDotFileName, "w") : stdout;
      fprintf(stderr, "Writing final decomposition to %s%s%s.\n", outputDotToFile ? "file <" : "",
        outputDotToFile ? outputDotFileName : "stdout", outputDotToFile ? ">" : "");

      char buffer[16];
      fputs("digraph G {\n", outputDotFile);
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
        fprintf(outputDotFile, " v_%d -> v_%d [label=\"%s\",style=bold,color=red];\n", u, v,
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
        fprintf(outputDotFile, " v_%d -> v_%d [label=\"%s\"];\n", u, v, CMRelementString(CMRcolumnToElement(column), buffer));
      }
      fputs("}\n", outputDotFile);

      if (outputDotToFile)
        fclose(outputDotFile);
    }

    if (edgesReversed)
      CMR_CALL( CMRfreeBlockArray(cmr, &edgesReversed) );
    CMR_CALL( CMRfreeBlockArray(cmr, &rowEdges) );
    CMR_CALL( CMRfreeBlockArray(cmr, &columnEdges) );
    CMR_CALL( CMRgraphFree(cmr, &digraph) );
  }
  else
  {
    if (outputSubmatrixFileName)
    {
      bool outputSubmatrixToFile = strcmp(outputSubmatrixFileName, "-");
      fprintf(stderr, "Writing minimal non-%snetwork submatrix to %s%s%s.\n", conetwork ? "co" : "",
        outputSubmatrixToFile ? "file <" : "", outputSubmatrixToFile ? outputSubmatrixFileName : "stdout",
        outputSubmatrixToFile ? ">" : "");

      assert(submatrix);
      CMR_CALL( CMRsubmatWriteToFile(cmr, submatrix, matrix->numRows, matrix->numColumns, outputSubmatrixFileName) );
    }
  }

  CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
  CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

/**
 * \brief Converts the given graph file to the corresponding (co)network matrix.
 */

CMR_ERROR computeNetwork(
  const char* inputFileName,        /**< File name of input graph (may be `-' for stdin). */
  const char* outputMatrixFileName, /**< File name of output matrix (may be `-' for stdout). */
  FileFormat outputFormat,          /**< Matrix format for output. */
  bool conetwork,                   /**< Whether the output shall be the conetwork matrix instead of the network matrix. */
  const char* inputTreeFileName,    /**< File name of input tree (may be `-' for stdin). */
  bool printStats                   /**< Whether to print statistics to stderr. */
)
{
  FILE* inputFile = strcmp(inputFileName, "-") ? fopen(inputFileName, "r") : stdin;
  if (!inputFile)
    return CMR_ERROR_INPUT;

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read edge list. */

  CMR_GRAPH* digraph = NULL;
  CMR_ELEMENT* edgeElements = NULL;
  CMR_CALL( CMRgraphCreateFromEdgeList(cmr, &digraph, &edgeElements, NULL, inputFile) );
  if (inputFile != stdin)
    fclose(inputFile);

  if (inputTreeFileName)
  {
    // TODO: implement
    assert(!"NOT IMPLEMENTED");
    return EXIT_FAILURE;
  }

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

  if (printStats)
  {
    clock_t endTime = clock();
    fprintf(stderr, "Time: %f\n", (endTime - startTime) * 1.0 / CLOCKS_PER_SEC);
  }

  bool outputMatrixToFile = strcmp(outputMatrixFileName, "-");
  FILE* outputMatrixFile = outputMatrixToFile ? fopen(outputMatrixFileName, "w") : stdout;
  fprintf(stderr, "Writing %snetwork matrix to %s%s%s in %s format.\n", conetwork ? "co" : "",
    outputMatrixToFile ? "file <" : "", outputMatrixToFile ? outputMatrixFileName : "stdout",
    outputMatrixToFile ? ">" : "", outputFormat == FILEFORMAT_MATRIX_DENSE ? "dense" : "sparse");

  if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRchrmatPrintDense(cmr, matrix, outputMatrixFile, '0', false) );
  else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRchrmatPrintSparse(cmr, matrix, outputMatrixFile) );
  else
    assert(false);

  if (outputMatrixToFile)
    fclose(outputMatrixFile);

  CMR_CALL( CMRchrmatFree(cmr, &matrix) );

  free(coforestEdges);
  free(forestEdges);
  free(edgeElements);
  CMR_CALL( CMRgraphFree(cmr, &digraph) );

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

/**
 * \brief Prints the usage of the \p program to stdout.
 * 
 * \returns \c EXIT_FAILURE.
 */

int printUsage(const char* program)
{
  fputs("Usage:\n", stderr);
  fprintf(stderr, "%s IN-MAT [OPTION]...\n\n", program);
  fputs("  (1) determines whether the matrix given in file IN-MAT is (co)network.\n\n", stderr);
  fprintf(stderr, "%s -c IN-GRAPH OUT-MAT [OPTION]...\n\n", program);
  fputs("  (2) computes a (co)network matrix corresponding to the digraph from file IN-GRAPH and writes it to OUT-MAT.\n\n\n",
    stderr);
  fputs("Options specific to (1):\n", stderr);
  fputs("  -i FORMAT    Format of file IN-MAT, among `dense' and `sparse'; default: dense.\n", stderr);
  fputs("  -t           Test for being conetwork; default: test for being network.\n", stderr);
  fputs("  -G OUT-GRAPH Write a digraph to file OUT-GRAPH; default: skip computation.\n", stderr);
  fputs("  -T OUT-TREE  Write a directed spanning tree to file OUT-TREE; default: skip computation.\n", stderr);
  fputs("  -D OUT-DOT   Write a dot file OUT-DOT with the digraph and the directed spanning tree; default: skip computation.\n", stderr);
  fputs("  -N NON-SUB   Write a minimal non-(co)network submatrix to file NON-SUB; default: skip computation.\n\n", stderr);
  fputs("Options specific to (2):\n", stderr);
  fputs("  -o FORMAT    Format of file OUT-MAT, among `dense' and `sparse'; default: dense.\n", stderr);
  fputs("  -t           Return the transpose of the network matrix.\n", stderr);
  fputs("  -T IN-TREE   Read a directed tree from file IN-TREE; default: use first specified arcs as tree edges.\n\n", stderr);
  fputs("Common options:\n", stderr);
  fputs("  -s           Print statistics about the computation to stderr.\n\n", stderr);
  fputs("Advanced options:\n", stderr);
  fputs("  --time-limit LIMIT   Allow at most LIMIT seconds for the computation.\n\n", stderr);
  fputs("If IN-MAT, IN-GRAPH or IN-TREE is `-' then the matrix (resp. the digraph or directed tree) is read from stdin.\n", stderr);
  fputs("If OUT-GRAPH, OUT-TREE, OUT-DOT or NON-SUB is `-' then the digraph (resp. the directed tree, dot file or non-(co)network submatrix) is written to stdout.\n",
    stderr);

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  Task task = TASK_RECOGNIZE;
  FileFormat inputFormat = FILEFORMAT_UNDEFINED;
  FileFormat outputFormat = FILEFORMAT_UNDEFINED;
  bool transposed = false;
  bool printStats = false;
  char* inputFileName = NULL;
  char* treeFileName = NULL;
  char* outputFileName = NULL;
  char* outputGraphFileName = NULL;
  char* outputDotFileName = NULL;
  char* outputSubmatrixFileName = NULL;
  double timeLimit = DBL_MAX;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-c"))
      task = TASK_COMPUTE;
    else if (!strcmp(argv[a], "-t"))
      transposed = true;
    else if (!strcmp(argv[a], "-s"))
      printStats = true;
    else if (!strcmp(argv[a], "-i") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        inputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        inputFormat = FILEFORMAT_MATRIX_SPARSE;
      else
      {
        fprintf(stderr, "Error: Unknown input file format <%s>.\n\n", argv[a+1]);
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
      else
      {
        fprintf(stderr, "Error: Unknown output format <%s>.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!strcmp(argv[a], "-G") && a+1 < argc)
      outputGraphFileName = argv[++a];
    else if (!strcmp(argv[a], "-T") && a+1 < argc)
      treeFileName = argv[++a];
    else if (!strcmp(argv[a], "-D") && a+1 < argc)
      outputDotFileName = argv[++a];
    else if (!strcmp(argv[a], "-N") && a+1 < argc)
      outputSubmatrixFileName = argv[++a];
    else if (!strcmp(argv[a], "--time-limit") && (a+1 < argc))
    {
      if (sscanf(argv[a+1], "%lf", &timeLimit) == 0 || timeLimit <= 0)
      {
        fprintf(stderr, "Error: Invalid time limit <%s> specified.\n\n", argv[a+1]);
        return printUsage(argv[0]);
      }
      ++a;
    }
    else if (!inputFileName)
      inputFileName = argv[a];
    else if (!outputFileName)
      outputFileName = argv[a];
    else
    {
      fprintf(stderr,
        "Error: Three file names <%s>, <%s> and <%s> specified.\n\n", inputFileName, outputFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  CMR_ERROR error;
  if (!inputFileName)
  {
    fputs("Error: No input file specified.\n", stderr);
    return printUsage(argv[0]);
  }
  else if (task == TASK_RECOGNIZE)
  {
    if (outputFileName)
    {
      fprintf(stderr,
        "Error: Two file names <%s> and <%s> specified for recognition.\n\n", inputFileName, outputFileName);
      return printUsage(argv[0]);
    }
    if (outputFormat != FILEFORMAT_UNDEFINED)
    {
      fprintf(stderr,
        "Error: Option -o is invalid for recognition.\n\n");
      return printUsage(argv[0]);
    }
    if (inputFormat == FILEFORMAT_UNDEFINED)
      inputFormat = FILEFORMAT_MATRIX_DENSE;

    error = recognizeNetwork(inputFileName, inputFormat, transposed, outputGraphFileName, treeFileName, outputDotFileName,
      outputSubmatrixFileName, printStats, timeLimit);
  }
  else if (task == TASK_COMPUTE)
  {
    if (!outputFileName)
    {
      fputs("Error: No output file specified.\n", stderr);
      return printUsage(argv[0]);
    }
    if (inputFormat != FILEFORMAT_UNDEFINED)
    {
      fprintf(stderr,
        "Error: Option -i is invalid for computation.\n\n");
      return printUsage(argv[0]);
    }
    if (outputGraphFileName)
    {
      fprintf(stderr,
        "Error: Option -G is invalid for computation.\n\n");
      return printUsage(argv[0]);
    }
    if (outputDotFileName)
    {
      fprintf(stderr,
        "Error: Option -D is invalid for computation.\n\n");
      return printUsage(argv[0]);
    }
    if (outputSubmatrixFileName)
    {
      fprintf(stderr,
        "Error: Option -N is invalid for computation.\n\n");
      return printUsage(argv[0]);
    }
    if (outputFormat == FILEFORMAT_UNDEFINED)
      outputFormat = FILEFORMAT_MATRIX_DENSE;

    error = computeNetwork(inputFileName, outputFileName, outputFormat, transposed, treeFileName, printStats);
  }
  else
  {
    assert(false);
  }

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
