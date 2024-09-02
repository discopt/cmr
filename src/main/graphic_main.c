#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>

#include <cmr/env.h>
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
  FILEFORMAT_MATRIX_SPARSE = 2  /**< Sparse matrix format. */
} FileFormat;

/**
 * \brief Converts matrix from a file to a graph if the former is (co)graphic.
 */

CMR_ERROR recognizeGraphic(
  const char* inputMatrixFileName,      /**< File name of the input matrix (may be `-' for stdin). */
  FileFormat inputFormat,               /**< Format of the input matrix. */
  bool cographic,                       /**< Whether the input shall be checked for being cographic instead of
                                         **  graphic. */
  const char* outputGraphFileName,      /**< File name of the output graph (may be NULL; may be `-' for stdout). */
  const char* outputTreeFileName,       /**< File name of the output tree (may be NULL; may be `-' for stdout). */
  const char* outputDotFileName,        /**< File name of the output dot file (may be NULL; may be `-' for stdout). */
  const char* outputSubmatrixFileName,  /**< File name of the output non-(co)graphic submatrix (may be NULL; may be `-'
                                         **  for stdout). */
  bool printStats,                      /**< Whether to print statistics to stderr. */
  double timeLimit                      /**< Time limit to impose. */
)
{
  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read matrix. */

  CMR_CHRMAT* matrix = NULL;
  clock_t readClock = clock();
  CMR_ERROR error = CMR_OKAY;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    error = CMRchrmatCreateFromDenseFile(cmr, inputMatrixFileName, "-", &matrix);
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    error = CMRchrmatCreateFromSparseFile(cmr, inputMatrixFileName, "-", &matrix);
  else
    CMR_CALL(CMR_ERROR_INVALID);

  if (error)
  {
    fprintf(stderr, "Input error: %s\n", CMRgetErrorMessage(cmr));
    CMR_CALL( CMRfreeEnvironment(&cmr) );
    return CMR_ERROR_INPUT;
  }

  fprintf(stderr, "Read %zux%zu matrix with %zu nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (clock() - readClock) * 1.0 / CLOCKS_PER_SEC);

  /* Test for being binary first. */

  CMR_SUBMAT* submatrix = NULL;
  if (!CMRchrmatIsBinary(cmr, matrix, &submatrix))
  {
    CMR_CHRMAT* mat = NULL;
    CMR_CALL( CMRchrmatSlice(cmr, matrix, submatrix, &mat) );
    assert(mat->numRows == 1);
    assert(mat->numColumns == 1);
    assert(mat->numNonzeros == 1);
    fprintf(stderr, "Matrix is NOT %sgraphic since it is not binary: entry at row %zu, column %zu is %d.\n",
      cographic ? "co" : "", submatrix->rows[0] + 1, submatrix->columns[0] + 1, mat->entryValues[0]);

    CMR_CALL( CMRchrmatFree(cmr, &mat) );
    CMR_CALL( CMRsubmatFree(cmr, &submatrix) );
    CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    CMR_CALL( CMRfreeEnvironment(&cmr) );

    return CMR_OKAY;
  }

  /* Test for (co)graphicness. */

  bool isCoGraphic;
  CMR_GRAPH* graph = NULL;
  CMR_GRAPH_EDGE* rowEdges = NULL;
  CMR_GRAPH_EDGE* columnEdges = NULL;
  bool* edgesReversed = NULL;
  CMR_GRAPHIC_STATISTICS stats;
  CMR_CALL( CMRgraphicStatsInit(&stats) );
  if (cographic)
  {
    CMR_CALL( CMRgraphicTestTranspose(cmr, matrix, &isCoGraphic, &graph, &columnEdges, &rowEdges,
      outputSubmatrixFileName ? &submatrix : NULL, &stats, timeLimit) );
  }
  else
  {
    CMR_CALL( CMRgraphicTestMatrix(cmr, matrix, &isCoGraphic, &graph, &rowEdges, &columnEdges,
      outputSubmatrixFileName ? &submatrix : NULL, &stats, timeLimit) );
  }

  fprintf(stderr, "Matrix %s%sgraphic.\n", isCoGraphic ? "IS " : "is NOT ", cographic ? "co" : "");
  if (printStats)
    CMR_CALL( CMRgraphicStatsPrint(stderr, &stats, NULL) );

  if (isCoGraphic)
  {
    if (outputGraphFileName)
    {
      bool outputGraphToFile = strcmp(outputGraphFileName, "-");
      FILE* outputGraphFile = outputGraphToFile ? fopen(outputGraphFileName, "w") : stdout;
      fprintf(stderr, "Writing %sgraph to %s%s%s.\n", cographic ? "co" : "", outputGraphToFile ? "file <" : "",
        outputGraphToFile ? outputGraphFileName : "stdout", outputGraphToFile ? ">" : "");

      if (cographic)
      {
        for (size_t column = 0; column < matrix->numColumns; ++column)
        {
          CMR_GRAPH_EDGE e = columnEdges[column];
          CMR_GRAPH_NODE u = CMRgraphEdgeU(graph, e);
          CMR_GRAPH_NODE v = CMRgraphEdgeV(graph, e);
          if (edgesReversed && edgesReversed[e])
          {
            CMR_GRAPH_NODE temp = u;
            u = v;
            v = temp;
          }
          fprintf(outputGraphFile, "%d %d c%zu\n", u, v, column+1);
        }
        for (size_t row = 0; row < matrix->numRows; ++row)
        {
          CMR_GRAPH_EDGE e = rowEdges[row];
          CMR_GRAPH_NODE u = CMRgraphEdgeU(graph, e);
          CMR_GRAPH_NODE v = CMRgraphEdgeV(graph, e);
          if (edgesReversed && edgesReversed[e])
          {
            CMR_GRAPH_NODE temp = u;
            u = v;
            v = temp;
          }
          fprintf(outputGraphFile, "%d %d r%zu\n", u, v, row+1);
        }
      }
      else
      {
        for (size_t row = 0; row < matrix->numRows; ++row)
        {
          CMR_GRAPH_EDGE e = rowEdges[row];
          CMR_GRAPH_NODE u = CMRgraphEdgeU(graph, e);
          CMR_GRAPH_NODE v = CMRgraphEdgeV(graph, e);
          if (edgesReversed && edgesReversed[e])
          {
            CMR_GRAPH_NODE temp = u;
            u = v;
            v = temp;
          }
          fprintf(outputGraphFile, "%d %d r%zu\n", u, v, row+1);
        }
        for (size_t column = 0; column < matrix->numColumns; ++column)
        {
          CMR_GRAPH_EDGE e = columnEdges[column];
          CMR_GRAPH_NODE u = CMRgraphEdgeU(graph, e);
          CMR_GRAPH_NODE v = CMRgraphEdgeV(graph, e);
          if (edgesReversed && edgesReversed[e])
          {
            CMR_GRAPH_NODE temp = u;
            u = v;
            v = temp;
          }
          fprintf(outputGraphFile, "%d %d c%zu\n", u, v, column+1);
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
      fprintf(stderr, "Writing %sgraph to %s%s%s.\n", cographic ? "co" : "", outputDotToFile ? "file <" : "",
        outputDotToFile ? outputDotFileName : "stdout", outputDotToFile ? ">" : "");

      char buffer[16];
      fputs("graph G {\n", outputDotFile);
      for (size_t row = 0; row < matrix->numRows; ++row)
      {
        CMR_GRAPH_EDGE e = rowEdges[row];
        CMR_GRAPH_NODE u = CMRgraphEdgeU(graph, e);
        CMR_GRAPH_NODE v = CMRgraphEdgeV(graph, e);
        if (edgesReversed && edgesReversed[e])
        {
          CMR_GRAPH_NODE temp = u;
          u = v;
          v = temp;
        }
        const char* style = cographic ? "" : ",style=bold,color=red";
        fprintf(outputDotFile, " v_%d -- v_%d [label=\"%s\"%s];\n", u, v,
          CMRelementString(CMRrowToElement(row), buffer), style);
      }
      for (size_t column = 0; column < matrix->numColumns; ++column)
      {
        CMR_GRAPH_EDGE e = columnEdges[column];
        CMR_GRAPH_NODE u = CMRgraphEdgeU(graph, e);
        CMR_GRAPH_NODE v = CMRgraphEdgeV(graph, e);
        if (edgesReversed && edgesReversed[e])
        {
          CMR_GRAPH_NODE temp = u;
          u = v;
          v = temp;
        }
        const char* style = cographic ? ",style=bold,color=red" : "";
        fprintf(outputDotFile, " v_%d -- v_%d [label=\"%s\"%s];\n", u, v,
          CMRelementString(CMRcolumnToElement(column), buffer), style);
      }
      fputs("}\n", outputDotFile);

      if (outputDotToFile)
        fclose(outputDotFile);
    }

    if (edgesReversed)
      CMR_CALL( CMRfreeBlockArray(cmr, &edgesReversed) );
    CMR_CALL( CMRfreeBlockArray(cmr, &rowEdges) );
    CMR_CALL( CMRfreeBlockArray(cmr, &columnEdges) );
    CMR_CALL( CMRgraphFree(cmr, &graph) );
  }
  else
  {
    if (outputSubmatrixFileName)
    {
      bool outputSubmatrixToFile = strcmp(outputSubmatrixFileName, "-");
      fprintf(stderr, "Writing minimal non-%sgraphic submatrix to %s%s%s.\n", cographic ? "co" : "",
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
 * \brief Converts the given graph file to the corresponding (co)graphic matrix.
 */

CMR_ERROR computeGraphic(
  const char* inputFileName,        /**< File name of input graph (may be `-' for stdin). */
  const char* outputMatrixFileName, /**< File name of output matrix (may be `-' for stdout). */
  FileFormat outputFormat,          /**< Matrix format for output. */
  bool cographic,                   /**< Whether the output shall be the cographic matrix instead of the graphic matrix. */
  const char* inputTreeFileName,    /**< File name of input tree (may be `-' for stdin). */
  bool printStats                   /**< Whether to print statistics to stderr. */
)
{
  FILE* inputFile = strcmp(inputFileName, "-") ? fopen(inputFileName, "r") : stdin;
  if (!inputFile)
  {
    fprintf(stderr, "Input error: Could not open file <%s>.\n", inputFileName);
    return CMR_ERROR_INPUT;
  }

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  /* Read edge list. */

  CMR_GRAPH* graph = NULL;
  CMR_ELEMENT* edgeElements = NULL;
  CMR_CALL( CMRgraphCreateFromEdgeList(cmr, &graph, &edgeElements, NULL, inputFile) );
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

  if (cographic)
  {
    CMR_CALL( CMRgraphicComputeMatrix(cmr, graph, NULL, &matrix, numForestEdges, forestEdges, numCoforestEdges,
      coforestEdges, &isCorrectForest) );
  }
  else
  {
    CMR_CALL( CMRgraphicComputeMatrix(cmr, graph, &matrix, NULL, numForestEdges, forestEdges, numCoforestEdges,
      coforestEdges, &isCorrectForest) );
  }

  if (printStats)
  {
    clock_t endTime = clock();
    fprintf(stderr, "Time: %f\n", (endTime - startTime) * 1.0 / CLOCKS_PER_SEC);
  }

  bool outputMatrixToFile = strcmp(outputMatrixFileName, "-");
  FILE* outputMatrixFile = outputMatrixToFile ? fopen(outputMatrixFileName, "w") : stdout;
  fprintf(stderr, "Writing %sgraphic matrix to %s%s%s in %s format.\n", cographic ? "co" : "",
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
  CMR_CALL( CMRgraphFree(cmr, &graph) );

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

  fprintf(stderr, "%s IN-MAT [OPTION]...\n", program);
  fputs("  (1) determines whether the matrix given in file IN-MAT is (co)graphic.\n", stderr);
  fputs("\n", stderr);

  fprintf(stderr, "%s -c IN-GRAPH OUT-MAT [OPTION]...\n", program);
  fputs("  (2) computes a (co)graphic matrix corresponding to the graph from file IN-GRAPH and writes it to OUT-MAT.\n",
    stderr);
  fputs("\n", stderr);

  fputs("Options specific to (1):\n", stderr);
  fputs("  -i FORMAT    Format of file IN-MAT; default: dense.\n", stderr);
  fputs("  -t           Test for being cographic; default: test for being graphic.\n", stderr);
  fputs("  -G OUT-GRAPH Write a graph to file OUT-GRAPH; default: skip computation.\n", stderr);
  fputs("  -T OUT-TREE  Write a spanning tree to file OUT-TREE; default: skip computation.\n", stderr);
  fputs("  -D OUT-DOT   Write a dot file OUT-DOT with the graph and the spanning tree; default: skip computation.\n",
    stderr);
  fputs("  -N NON-SUB   Write a minimal non-(co)graphic submatrix to file NON-SUB; default: skip computation.\n",
    stderr);
  fputs("\n", stderr);

  fputs("Options specific to (2):\n", stderr);
  fputs("  -o FORMAT    Format of file OUT-MAT; default: dense.\n", stderr);
  fputs("  -t           Return the transpose of the graphic matrix.\n", stderr);
  fputs("  -T IN-TREE   Read a tree from file IN-TREE; default: use first specified edges as tree edges.\n", stderr);
  fputs("\n", stderr);

  fputs("Advanced options:\n", stderr);
  fputs("  --stats            Print statistics about the computation to stderr.\n", stderr);
  fputs("  --time-limit LIMIT Allow at most LIMIT seconds for the computation.\n", stderr);
  fputs("\n", stderr);

  fputs("Formats for matrices: dense, sparse\n", stderr);
  fputs("If IN-MAT, IN-GRAPH or IN-TREE is `-' then the matrix (resp. the graph or tree) is read from stdin.\n",
    stderr);
  fputs("If OUT-GRAPH, OUT-TREE, OUT-DOT or NON-SUB is `-' then the graph (resp. the tree, dot file or non-(co)graphic"
    " submatrix) is written to stdout.\n", stderr);

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
    else if (!strcmp(argv[a], "--stats"))
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

    error = recognizeGraphic(inputFileName, inputFormat, transposed, outputGraphFileName, treeFileName,
      outputDotFileName, outputSubmatrixFileName, printStats, timeLimit);
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

    error = computeGraphic(inputFileName, outputFileName, outputFormat, transposed, treeFileName, printStats);
  }
  else
  {
    assert(false);
  }

  switch (error)
  {
  case CMR_ERROR_INPUT:
    /* The actual function will have reported the details. */
    return EXIT_FAILURE;
  case CMR_ERROR_MEMORY:
    puts("Memory error.");
    return EXIT_FAILURE;
  default:
    return EXIT_SUCCESS;
  }
}
