#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>
#include <math.h>

#include <cmr/graphic.h>

static inline
size_t randRange(size_t first, size_t beyond)
{
  size_t N = beyond - first;
  size_t representatives = (RAND_MAX + 1u) / N;
  size_t firstInvalid = N * representatives;
  size_t x;
  do
  {
    x = rand();
  }
  while (x >= firstInvalid);
  return first + x / representatives;
}

int printUsage(const char* program)
{
  fprintf(stderr, "Usage: %s [OPTIONS] ROWS COLS\n\n", program);
  fputs("Creates a random ROWS-by-COLS 0/1 graphic matrix.\n", stderr);
  fputs("Options:\n", stderr);
  fputs("  -b NUM       Benchmarks the recognition algorithm for the created matrix with NUM repetitions.\n", stderr);
  fputs("Notes:", stderr);
  fputs("  -p and -s cannot be specified at the same time.\n", stderr);
  return EXIT_FAILURE;
}

int compare(const void* pa, const void* pb)
{
  size_t a = *((size_t*)(pa));
  size_t b = *((size_t*)(pb));
  return a < b ? -1 : (a > b);
}

CMR_ERROR genMatrixGraphic(
  size_t numRows,         /**< Number of rows of base matrix. */
  size_t numColumns,      /**< Number of columns of base matrix. */
  size_t benchmarkRepetitions /**< Whether to benchmark the recognition algorithm with the matrix instead of printing it. */
)
{
  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  size_t numBenchmarkNonzeros = 0;
  size_t numNodes = numRows + 1;
  size_t numEdges = numColumns;
  CMR_GRAPHIC_STATISTICS stats;
  CMR_CALL( CMRgraphicInitStatistics(&stats) );
  for (size_t benchmark = benchmarkRepetitions ? benchmarkRepetitions : 1; benchmark > 0; --benchmark)
  {
    clock_t startTime = clock();

    /* Init transpose of matrix. */
    size_t transposedMemNonzeros = 1;
    for (size_t x = numRows; x; x >>= 1)
      ++transposedMemNonzeros;
    transposedMemNonzeros *= numColumns;
    CMR_CHRMAT* transposed = NULL;
    CMR_CALL( CMRchrmatCreate(cmr, &transposed, numEdges, numNodes-1, transposedMemNonzeros) );
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

    size_t* columnNonzeros = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &columnNonzeros, numNodes - 1) );
    for (int e = 0; e < numEdges; ++e)
    {
      size_t numColumNonzeros = 0;
      int first = (int)(rand() * 1.0 * numNodes / RAND_MAX);
      int second = (int)(rand() * 1.0 * numNodes / RAND_MAX);
      while (treeDistance[first] > treeDistance[second])
      {
        columnNonzeros[numColumNonzeros++] = first-1;
        first = nextTreeNode[first];
      }
      while (treeDistance[second] > treeDistance[first])
      {
        columnNonzeros[numColumNonzeros++] = second-1;
        second = nextTreeNode[second];
      }
      while (first != second && first)
      {
        columnNonzeros[numColumNonzeros++] = first-1;
        first = nextTreeNode[first];
        columnNonzeros[numColumNonzeros++] = second-1;
        second = nextTreeNode[second];
      }
      transposed->rowSlice[e] = transposed->numNonzeros;

      qsort( columnNonzeros, numColumNonzeros, sizeof(size_t), &compare);

      for (size_t i = 0; i < numColumNonzeros; ++i)
      {
        size_t row = columnNonzeros[i];
        if (transposed->numNonzeros == transposedMemNonzeros)
        {
          transposedMemNonzeros *= 2;
          CMR_CALL( CMRreallocBlockArray(cmr, &transposed->entryColumns, transposedMemNonzeros) );
          CMR_CALL( CMRreallocBlockArray(cmr, &transposed->entryValues, transposedMemNonzeros) );
        }
        transposed->entryColumns[transposed->numNonzeros] = row;
        transposed->entryValues[transposed->numNonzeros] = 1;
        transposed->numNonzeros++;
      }
    }
    transposed->rowSlice[transposed->numRows] = transposed->numNonzeros;

    CMR_CALL( CMRfreeBlockArray(cmr, &columnNonzeros) );
    CMR_CALL( CMRfreeBlockArray(cmr, &nextTreeNode) );
    CMR_CALL( CMRfreeBlockArray(cmr, &treeDistance) );

    CMR_CHRMAT* matrix = NULL;
    CMR_CALL( CMRchrmatTranspose(cmr, transposed, &matrix) );
    CMR_CALL( CMRchrmatFree(cmr, &transposed) );

    double generationTime = (clock() - startTime) * 1.0 / CLOCKS_PER_SEC;
    fprintf(stderr, "Generated a %ldx%ld matrix with %ld nonzeros in %f seconds.\n", numRows, numColumns,
      matrix->numNonzeros, generationTime);

    if (benchmarkRepetitions)
    {
      /* Benchmark */      
      bool isGraphic;
      CMR_CALL( CMRtestGraphicMatrix(cmr, matrix, &isGraphic, NULL, NULL, NULL, NULL, &stats) );
    }
    else
    {
      /* Print matrix. */
      CMR_CALL( CMRchrmatPrintDense(cmr, matrix, stdout, '0', false) );
    }

    /* Cleanup. */
    CMR_CALL( CMRchrmatFree(cmr, &matrix) );
  }

  CMRgraphicPrintStatistics(stderr, &stats);

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

int main(int argc, char** argv)
{
  struct timeval curTime;
  gettimeofday(&curTime, NULL);
  srand(curTime.tv_usec);

  size_t numRows = SIZE_MAX;
  size_t numColumns = SIZE_MAX;
  size_t benchmarkRepetitions = 0;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-b") && a+1 < argc)
    {
      char* p;
      benchmarkRepetitions = strtoull(argv[a+1], &p, 10);
      if (*p != '\0' || benchmarkRepetitions == 0)
      {
        fprintf(stderr, "Error: invalid number of benchmark repetitions <%s>", argv[a+1]);
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a++;
    }
    else if (numRows == SIZE_MAX)
    {
      char* p = NULL;
      numRows = strtoull(argv[a], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
    }
    else if (numColumns == SIZE_MAX)
    {
      char* p = NULL;
      numColumns = strtoull(argv[a], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
    }
    else
    {
      fprintf(stderr, "Error: more than two size indicators specified: %ld %ld %s\n\n", numRows, numColumns, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (numRows == SIZE_MAX)
  {
    fputs("Error: no size indicator specified.\n", stderr);
    return printUsage(argv[0]);
  }
  else if (numColumns == SIZE_MAX)
  {
    fputs("Error: only one size indicator specified.\n", stderr);
    return printUsage(argv[0]);
  }
  else if (numRows <= 0 || numColumns <= 0)
  {
    fputs("Error: matrix must have at least 1 row and 1 column.\n", stderr);
    return printUsage(argv[0]);
  }

  CMR_ERROR error = genMatrixGraphic(numRows, numColumns, benchmarkRepetitions);
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
