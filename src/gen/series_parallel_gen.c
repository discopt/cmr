#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <sys/time.h>
#include <stdint.h>

#include <cmr/tu.h>
#include <cmr/series_parallel.h>

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
  fputs("Creates a random ROWS-by-COLS 0/1 or -1/0/+1 base matrix and augments it using series-parallel operations.\n", stderr);
  fputs("Options:\n", stderr);
  fputs("  -z ROWS COLS Add ROWS zero rows and COLS zero columns (default: 0 0).\n", stderr);
  fputs("  -u ROWS COLS Add ROWS unit rows and COLS unit columns (default: 0 0).\n", stderr);
  fputs("  -c ROWS COLS Add ROWS copied rows and COLS copied columns (default: 0 0).\n", stderr);
  fputs("  -t           Create a ternary matrix.\n", stderr);
  fputs("  -p PROB      Sets the probability of a nonzero in the base matrix (default: 0.5).\n", stderr);
  fputs("  -s SPARSITY  Sets the probability of a nonzero in the base matrix minimum such that each\n", stderr);
  fputs("               row/column has at least SPARSITY 1's in expectation.\n", stderr);
  fputs("  -r           Randomize matrix by permuting rows and columns afterwards (default: false).\n", stderr);
  fputs("  -b NUM       Benchmarks the recognition algorithm for the created matrix with NUM repetitions.\n", stderr);
  fputs("Notes:", stderr);
  fputs("  -p and -s cannot be specified at the same time.\n", stderr);
  return EXIT_FAILURE;
}

typedef struct _ListNonzero
{
  struct _ListNonzero* left;
  struct _ListNonzero* right;
  struct _ListNonzero* above;
  struct _ListNonzero* below;
  size_t row;
  size_t column;
  char value;
} ListNonzero;

static
CMR_ERROR addNonzero(
  CMR* cmr,
  ListNonzero* rowHeads,
  ListNonzero* columnHeads,
  size_t* pnumNonzeros,
  size_t row,
  size_t column,
  char value
)
{
  assert(cmr);

  ListNonzero* nz = NULL;
  CMR_CALL( CMRallocBlock(cmr, &nz) );
  nz->right = &rowHeads[row];
  nz->left = rowHeads[row].left;
  nz->left->right = nz;
  nz->right->left = nz;
  nz->below = &columnHeads[column];
  nz->above = columnHeads[column].above;
  nz->below->above = nz;
  nz->above->below = nz;
  nz->row = row;
  nz->column = column;
  nz->value = value;
  (*pnumNonzeros)++;
  
  return CMR_OKAY;
}

typedef struct
{
  size_t row;
  size_t column;
  char value;
} Nonzero;

int compareNonzeros(const void* a, const void* b)
{
  Nonzero* nza = (Nonzero*) a;
  Nonzero* nzb = (Nonzero*) b;
  if (nza->row < nzb->row)
    return -1;
  else if (nza->row > nzb->row)
    return +1;
  else
    return (int) nza->column - (int) nzb->column;
}

CMR_ERROR genMatrixSeriesParallel(
  size_t numBaseRows,         /**< Number of rows of base matrix. */
  size_t numBaseColumns,      /**< Number of columns of base matrix. */
  size_t numZeroRows,         /**< Number of added zero rows. */
  size_t numZeroColumns,      /**< Number of added zero columns. */
  size_t numUnitRows,         /**< Number of added unit rows. */
  size_t numUnitColumns,      /**< Number of added unit columns. */
  size_t numCopiedRows,       /**< Number of added copied rows. */
  size_t numCopiedColumns,    /**< Number of added copied columns. */
  bool ternary,               /**< Whether to create a ternary matrix. */
  double probability,         /**< Probability for each entry of base matrix to be a 1. */
  bool randomize,             /**< Whether to randomize afterwards via row/column permutations. */
  size_t benchmarkRepetitions /**< Whether to benchmark the recognition algorithm with the matrix instead of printing it. */
)
{
  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );

  size_t numTotalRows = numBaseRows + numZeroRows + numUnitRows + numCopiedRows;
  size_t numTotalColumns = numBaseColumns + numZeroColumns + numUnitColumns + numCopiedColumns;

  CMR_SP_STATISTICS stats;
  CMR_CALL( CMRstatsSeriesParallelInit(&stats) );
  size_t numBenchmarkNonzeros = 0;
  for (size_t benchmark = benchmarkRepetitions ? benchmarkRepetitions : 1; benchmark > 0; --benchmark)
  {
    size_t totalMemory = 0;
    ListNonzero* rowHeads = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &rowHeads, numTotalRows) );
    totalMemory += sizeof(ListNonzero) * numTotalRows;
    ListNonzero* columnHeads = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &columnHeads, numTotalColumns) );
    totalMemory += sizeof(ListNonzero) * numTotalColumns;

    for (size_t row = 0; row < numTotalRows; ++row)
    {
      rowHeads[row].left = &rowHeads[row];
      rowHeads[row].right = &rowHeads[row];
      rowHeads[row].above = NULL;
      rowHeads[row].below = NULL;
      rowHeads[row].row = row;
      rowHeads[row].column = SIZE_MAX;
    }
    for (size_t column = 0; column < numTotalColumns; ++column)
    {
      columnHeads[column].left = NULL;
      columnHeads[column].right = NULL;
      columnHeads[column].below = &columnHeads[column];
      columnHeads[column].above = &columnHeads[column];
      columnHeads[column].column = column;
      columnHeads[column].row = SIZE_MAX;
    }

    /* Create base matrix. */
    size_t numBaseNonzeros = 0;
    for (size_t row = 0; row < numBaseRows; ++row)
    {
      for (size_t column = 0; column < numBaseColumns; ++column)
      {
        if ((rand() * 1.0 / RAND_MAX) > probability)
          continue;

        char sign = (ternary && (rand() * 1.0 / RAND_MAX >= 0.5)) ? -1 : 1;
        CMR_CALL( addNonzero(cmr, rowHeads, columnHeads, &numBaseNonzeros, row, column, sign) );
        totalMemory += sizeof(ListNonzero);
      }
    }
    size_t numTotalNonzeros = numBaseNonzeros;

    /* Create a list of all operations. */
    char* operations = NULL;
    size_t numOperations = numZeroRows + numZeroColumns + numUnitRows + numUnitColumns + numCopiedRows + numCopiedColumns;
    CMR_CALL( CMRallocBlockArray(cmr, &operations, numOperations) );
    totalMemory += sizeof(char) * numOperations;
    char* op = operations;
    for (size_t i = 0; i < numZeroRows; ++i)
      (*op++) = 'z';
    for (size_t i = 0; i < numZeroColumns; ++i)
      (*op++) = 'Z';
    for (size_t i = 0; i < numUnitRows; ++i)
      (*op++) = 'u';
    for (size_t i = 0; i < numUnitColumns; ++i)
      (*op++) = 'U';
    for (size_t i = 0; i < numCopiedRows; ++i)
      (*op++) = 'c';
    for (size_t i = 0; i < numCopiedColumns; ++i)
      (*op++) = 'C';
    assert(op == &operations[numOperations]);

    /* Shuffle operations array. */
    for (size_t i = 0; i+1 < numOperations; ++i)
    {
      size_t j = randRange(i, numOperations);
      char tmp = operations[i];
      operations[i] = operations[j];
      operations[j] = tmp;
    }

    /* Start applying operations. */
    size_t numRows = numBaseRows;
    size_t numColumns = numBaseColumns;
    for (size_t i = 0; i < numOperations; ++i)
    {
      switch(operations[i])
      {
        case 'z':
        {
          ++numRows;
          break;
        }
        case 'Z':
        {
          ++numColumns;
          break;
        }
        case 'u':
        {
          size_t column = randRange(0, numColumns);
          char sign = (ternary && (rand() * 1.0 / RAND_MAX >= 0.5)) ? -1 : 1;
          CMR_CALL( addNonzero(cmr, rowHeads, columnHeads, &numTotalNonzeros, numRows, column, sign) );
          totalMemory += sizeof(ListNonzero);
          ++numRows;
          break;
        }
        case 'U':
        {
          size_t row = randRange(0, numRows);
          char sign = (ternary && (rand() * 1.0 / RAND_MAX >= 0.5)) ? -1 : 1;
          CMR_CALL( addNonzero(cmr, rowHeads, columnHeads, &numTotalNonzeros, row, numColumns, sign) );
          totalMemory += sizeof(ListNonzero);
          ++numColumns;
          break;
        }
        case 'c':
        {
          size_t row = randRange(0, numRows);
          char sign = (ternary && (rand() * 1.0 / RAND_MAX >= 0.5)) ? -1 : 1;
          for (ListNonzero* nz = rowHeads[row].right; nz->column != SIZE_MAX; nz = nz->right)
          {
            CMR_CALL( addNonzero(cmr, rowHeads, columnHeads, &numTotalNonzeros, numRows, nz->column, sign * nz->value) );
            totalMemory += sizeof(ListNonzero);
          }
          ++numRows;
          break;
        }
        case 'C':
        {
          size_t column = randRange(0, numColumns);
          char sign = (ternary && (rand() * 1.0 / RAND_MAX >= 0.5)) ? -1 : 1;
          for (ListNonzero* nz = columnHeads[column].below; nz->row != SIZE_MAX; nz = nz->below)
          {
            CMR_CALL( addNonzero(cmr, rowHeads, columnHeads, &numTotalNonzeros, nz->row, numColumns, sign * nz->value) );
            totalMemory += sizeof(ListNonzero);
          }
          ++numColumns;
          break;
        }
      }
    }

    /* Create permutations. */
    size_t* rowPermutation = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &rowPermutation, numTotalRows) );
    totalMemory += sizeof(size_t) * numTotalRows;
    for (size_t row = 0; row < numTotalRows; ++row)
      rowPermutation[row] = row;
    size_t* columnPermutation = NULL;
    CMR_CALL( CMRallocBlockArray(cmr, &columnPermutation, numTotalColumns) );
    totalMemory += sizeof(size_t) * numTotalColumns;
    for (size_t column = 0; column < numTotalColumns; ++column)
      columnPermutation[column] = column;
    if (randomize)
    {
      for (size_t row = 0; row < numTotalRows; ++row)
      {
        size_t r = randRange(row, numTotalRows);
        size_t tmp = rowPermutation[row];
        rowPermutation[row] = rowPermutation[r];
        rowPermutation[r] = tmp;
      }
      for (size_t column = 0; column < numTotalColumns; ++column)
      {
        size_t c = randRange(column, numTotalColumns);
        size_t tmp = columnPermutation[column];
        columnPermutation[column] = columnPermutation[c];
        columnPermutation[c] = tmp;
      }
    }

    fprintf(stderr, "Generated a %zux%zu matrix with %zu nonzeros.\n", numTotalRows, numTotalColumns, numTotalNonzeros);
    fprintf(stderr, "It contains a %zux%zu base matrix with %zu nonzeros, 1-entries generated with probability %g.\n",
      numBaseRows, numBaseColumns, numBaseNonzeros, probability);
    fprintf(stderr, "Series-parallel operations: %zux%zu zero, %zux%zu unit, %zux%zu copied\n", numZeroRows,
      numZeroColumns, numUnitRows, numUnitColumns, numCopiedRows, numCopiedColumns);
    if (randomize)
      fputs("Random row and column permutations were applied.\n", stderr);
    fprintf(stderr, "Amount of memory used: %.2g GB.\n", totalMemory / (1024.0*1024.0*1024.0));

    CMR_CALL( CMRfreeBlockArray(cmr, &operations) );

    if (benchmarkRepetitions)
    {
      /* Create array of all nonzeros. */
      
      Nonzero* nzs = NULL;
      CMR_CALL( CMRallocBlockArray(cmr, &nzs, numTotalNonzeros) );
      size_t* numRowNonzeros = NULL;
      CMR_CALL( CMRallocBlockArray(cmr, &numRowNonzeros, numRows) );
      for (size_t row = 0; row < numTotalRows; ++row)
        numRowNonzeros[row] = 0;
      size_t i = 0;
      for (size_t row = 0; row < numTotalRows; ++row)
      {
        for (ListNonzero* nz = rowHeads[row].right; nz->column != SIZE_MAX; nz = nz->right)
        {
          nzs[i].row = rowPermutation[nz->row];
          nzs[i].column = columnPermutation[nz->column];
          nzs[i].value = nz->value;
          numRowNonzeros[nzs[i].row]++;
          ++i;
        }
      }

      qsort(nzs, numTotalNonzeros, sizeof(Nonzero), compareNonzeros);

      /* Create matrix. */

      CMR_CHRMAT* matrix = NULL;
      CMR_CALL( CMRchrmatCreate(cmr, &matrix, numTotalRows, numTotalColumns, numTotalNonzeros) );

      size_t row = 0;
      for (size_t i = 0; i < numTotalNonzeros; ++i)
      {
        while (row <= nzs[i].row)
          matrix->rowSlice[row++] = i;

        matrix->entryColumns[i] = nzs[i].column;
        matrix->entryValues[i] = nzs[i].value;
      }
      while (row <= matrix->numRows)
        matrix->rowSlice[row++] = matrix->numNonzeros;

      CMR_CALL( CMRfreeBlockArray(cmr, &numRowNonzeros) );
      CMR_CALL( CMRfreeBlockArray(cmr, &nzs) );

      CMR_SP_REDUCTION* reductions = NULL;
      CMR_CALL( CMRallocBlockArray(cmr, &reductions, matrix->numRows + matrix->numColumns) );
      size_t numReductions;

      /* Benchmark */

      CMR_SUBMAT* wheelMatrix = NULL;
      if (ternary)
      {
        CMR_CALL( CMRtestTernarySeriesParallel(cmr, matrix, NULL, reductions, &numReductions, NULL, &wheelMatrix,
          &stats, DBL_MAX) );
      }
      else
      {
        CMR_CALL( CMRtestBinarySeriesParallel(cmr, matrix, NULL, reductions, &numReductions, NULL, &wheelMatrix,
          &stats, DBL_MAX) );
      }
      numBenchmarkNonzeros += numTotalNonzeros;

      CMR_CALL( CMRsubmatFree(cmr, &wheelMatrix) );
      CMR_CALL( CMRfreeBlockArray(cmr, &reductions) );
      CMR_CALL( CMRchrmatFree(cmr, &matrix) );
    }
    else
    {
      /* Print matrix. */

      printf("%zu %zu %zu\n\n", numTotalRows, numTotalColumns, numTotalNonzeros);
      for (size_t row = 0; row < numTotalRows; ++row)
      {
        for (ListNonzero* nz = rowHeads[row].right; nz->column != SIZE_MAX; nz = nz->right)
        {
          printf("%zu %zu %d\n", rowPermutation[nz->row]+1, columnPermutation[nz->column]+1, nz->value);
        }
      }
    }

    /* Cleanup. */

    CMR_CALL( CMRfreeBlockArray(cmr, &columnPermutation) );
    CMR_CALL( CMRfreeBlockArray(cmr, &rowPermutation) );
    for (size_t row = 0; row < numTotalRows; ++row)
    {
      for (ListNonzero* nz = rowHeads[row].right; nz != &rowHeads[row]; )
      {
        ListNonzero* current = nz;
        nz = nz->right;
        CMR_CALL( CMRfreeBlock(cmr, &current) );
      }
    }

    CMR_CALL( CMRfreeBlockArray(cmr, &columnHeads) );
    CMR_CALL( CMRfreeBlockArray(cmr, &rowHeads) );
  }

  if (benchmarkRepetitions > 0)
  {
    CMR_CALL( CMRstatsSeriesParallelPrint(stderr, &stats, NULL) );
    printf("Average number of nonzeros:     %f\n", (double)(numBenchmarkNonzeros) * 1.0 / benchmarkRepetitions);
  }

  CMR_CALL( CMRfreeEnvironment(&cmr) );

  return CMR_OKAY;
}

int main(int argc, char** argv)
{
  struct timeval curTime;
  gettimeofday(&curTime, NULL);
  srand(curTime.tv_usec);

  size_t numBaseRows = SIZE_MAX;
  size_t numBaseColumns = SIZE_MAX;
  size_t numZeroRows = 0;
  size_t numZeroColumns = 0;
  size_t numUnitRows = 0;
  size_t numUnitColumns = 0;
  size_t numCopiedRows = 0;
  size_t numCopiedColumns = 0;
  bool ternary = false;
  double probability = -1.0;
  double sparsity = -1.0;
  bool randomize = false;
  size_t benchmarkRepetitions = 0;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-z") && a+2 < argc)
    {
      char* p;
      numZeroRows = strtoull(argv[a+1], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      numZeroColumns = strtoull(argv[a+2], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a+= 2;
    }
    else if (!strcmp(argv[a], "-u") && a+2 < argc)
    {
      char* p;
      numUnitRows = strtoull(argv[a+1], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      numUnitColumns = strtoull(argv[a+2], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a+= 2;
    }
    else if (!strcmp(argv[a], "-c") && a+2 < argc)
    {
      char* p;
      numCopiedRows = strtoull(argv[a+1], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      numCopiedColumns = strtoull(argv[a+2], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a+= 2;
    }
    else if (!strcmp(argv[a], "-t"))
      ternary = true;
    else if (!strcmp(argv[a], "-p") && a+1 < argc)
    {
      char* p;
      probability = strtod(argv[a+1], &p);
      if (*p != '\0' || probability < 0.0 || probability > 1.0)
      {
        fprintf(stderr, "Error: invalid probablity <%s>", argv[a+1]);
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a++;
    }
    else if (!strcmp(argv[a], "-s") && a+1 < argc)
    {
      char* p;
      sparsity = strtod(argv[a+1], &p);
      if (*p != '\0' || sparsity < 0.0)
      {
        fprintf(stderr, "Error: invalid sparsity <%s>", argv[a+1]);
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
      a++;
    }
    else if (!strcmp(argv[a], "-r"))
      randomize = true;
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
    else if (numBaseRows == SIZE_MAX)
    {
      char* p = NULL;
      numBaseRows = strtoull(argv[a], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
    }
    else if (numBaseColumns == SIZE_MAX)
    {
      char* p = NULL;
      numBaseColumns = strtoull(argv[a], &p, 10);
      if (*p != '\0')
      {
        printUsage(argv[0]);
        return EXIT_FAILURE;
      }
    }
    else
    {
      fprintf(stderr, "Error: more than two size indicators specified: %zu %zu %s\n\n", numBaseRows, numBaseColumns,
        argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (numBaseRows == SIZE_MAX)
  {
    fputs("Error: no size indicator specified.\n", stderr);
    return printUsage(argv[0]);
  }
  else if (numBaseColumns == SIZE_MAX)
  {
    fputs("Error: only one size indicator specified.\n", stderr);
    return printUsage(argv[0]);
  }
  else if (numBaseRows <= 0 || numBaseColumns <= 0)
  {
    fputs("Error: base matrix must have at least 1 row and 1 column.\n", stderr);
    return printUsage(argv[0]);
  }
  else if (sparsity >= 0.0 && probability >= 0.0)
  {
    fputs("Error: at most one of -p and -s can be specified.\n", stderr);
    return printUsage(argv[0]);
  }

  if (sparsity >= 0.0)
  {
    size_t smallerSize = numBaseRows < numBaseColumns ? numBaseRows : numBaseColumns;
    probability = sparsity * 1.0 / smallerSize;
  }
  if (probability < 0.0)
    probability = 0.5;
  if (probability > 1.0)
    probability = 1.0;

  CMR_ERROR error = genMatrixSeriesParallel(numBaseRows, numBaseColumns, numZeroRows, numZeroColumns, numUnitRows,
    numUnitColumns, numCopiedRows, numCopiedColumns, ternary, probability, randomize, benchmarkRepetitions);
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
