#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>

#include <tu/tu.h>
#include <tu/series_parallel.h>

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
  fputs("Creates a random ROWS-by-COLS base matrix and augments it using series-parallel operations.\n", stderr);
  fputs("Options:\n", stderr);
  fputs("  -z ROWS COLS Add ROWS zero rows and COLS zero columns (default: 0 0).\n", stderr);
  fputs("  -u ROWS COLS Add ROWS unit rows and COLS unit columns (default: 0 0).\n", stderr);
  fputs("  -c ROWS COLS Add ROWS copied rows and COLS copied columns (default: 0 0).\n", stderr);
  fputs("  -p PROB      Sets the probability of a '1' in the base matrix (default: 0.5).\n", stderr);
  fputs("  -s SPARSITY  Sets the probability of a '1' in the base matrix minimum such that each\n", stderr);
  fputs("               row/column has at least SPARSITY 1's in expectation.\n", stderr);
  fputs("  -r           Randomize matrix by permuting rows and columns afterwards (default: false).\n", stderr);
  fputs("Notes:", stderr);
  fputs("  -p and -s cannot be specified at the same time.\n", stderr);
  return EXIT_FAILURE;
}

typedef struct _Nonzero
{
  struct _Nonzero* left;
  struct _Nonzero* right;
  struct _Nonzero* above;
  struct _Nonzero* below;
  size_t row;
  size_t column;
} Nonzero;

static
TU_ERROR addNonzero(
  TU* tu,
  Nonzero* rowHeads,
  Nonzero* columnHeads,
  size_t* pnumNonzeros,
  size_t row,
  size_t column
)
{
  assert(tu);

  Nonzero* nz = NULL;
  TU_CALL( TUallocBlock(tu, &nz) );
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
  (*pnumNonzeros)++;
  
  return TU_OKAY;
}

TU_ERROR genMatrixSeriesParallel(
  size_t numBaseRows,       /**< Number of rows of base matrix. */
  size_t numBaseColumns,    /**< Number of columns of base matrix. */
  size_t numZeroRows,       /**< Number of added zero rows. */
  size_t numZeroColumns,    /**< Number of added zero columns. */
  size_t numUnitRows,       /**< Number of added unit rows. */
  size_t numUnitColumns,    /**< Number of added unit columns. */
  size_t numCopiedRows,     /**< Number of added copied rows. */
  size_t numCopiedColumns,  /**< Number of added copied columns. */
  double probability,       /**< Probability for each entry of base matrix to be a 1. */
  bool randomize            /**< Whether to randomize afterwards via row/column permutations. */
)
{
  TU* tu = NULL;
  TU_CALL( TUcreateEnvironment(&tu) );

  size_t numTotalRows = numBaseRows + numZeroRows + numUnitRows + numCopiedRows;
  size_t numTotalColumns = numBaseColumns + numZeroColumns + numUnitColumns + numCopiedColumns;

  Nonzero* rowHeads = NULL;
  TU_CALL( TUallocBlockArray(tu, &rowHeads, numTotalRows) );
  for (size_t row = 0; row < numTotalRows; ++row)
  {
    rowHeads[row].left = &rowHeads[row];
    rowHeads[row].right = &rowHeads[row];
    rowHeads[row].above = NULL;
    rowHeads[row].below = NULL;
    rowHeads[row].row = row;
    rowHeads[row].column = SIZE_MAX;
  }

  Nonzero* columnHeads = NULL;
  TU_CALL( TUallocBlockArray(tu, &columnHeads, numTotalRows) );
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

      TU_CALL( addNonzero(tu, rowHeads, columnHeads, &numBaseNonzeros, row, column) );
    }
  }
  size_t numTotalNonzeros = numBaseNonzeros;

  /* Create a list of all operations. */
  char* operations = NULL;
  size_t numOperations = numZeroRows + numZeroColumns + numUnitRows + numUnitColumns + numCopiedRows + numCopiedColumns;
  TU_CALL( TUallocBlockArray(tu, &operations, numOperations) );
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
        TU_CALL( addNonzero(tu, rowHeads, columnHeads, &numTotalNonzeros, numRows, column) );
        ++numRows;
        break;
      }
      case 'U':
      {
        size_t row = randRange(0, numRows);
        TU_CALL( addNonzero(tu, rowHeads, columnHeads, &numTotalNonzeros, row, numColumns) );
        ++numColumns;
        break;
      }
      case 'c':
      {
        size_t row = randRange(0, numRows);
        for (Nonzero* nz = rowHeads[row].right; nz->column != SIZE_MAX; nz = nz->right)
          TU_CALL( addNonzero(tu, rowHeads, columnHeads, &numTotalNonzeros, numRows, nz->column) );
        ++numRows;
        break;
      }
      case 'C':
      {
        size_t column = randRange(0, numColumns);
        for (Nonzero* nz = columnHeads[column].below; nz->row != SIZE_MAX; nz = nz->below)
          TU_CALL( addNonzero(tu, rowHeads, columnHeads, &numTotalNonzeros, nz->row, numColumns) );
        ++numColumns;
        break;
      }
    }
  }

  fprintf(stderr, "Generated a %ldx%ld matrix with %ld nonzeros.\n", numTotalRows, numTotalColumns, numTotalNonzeros);
  fprintf(stderr, "It contains a %ldx%ld base matrix with %ld nonzeros, 1-entries generated with probability %g.\n",
    numBaseRows, numBaseColumns, numBaseNonzeros, probability);
  fprintf(stderr, "Series-parallel operations: %ldx%ld zero, %ldx%ld unit, %ldx%ld copied\n", numZeroRows,
    numZeroColumns, numUnitRows, numUnitColumns, numCopiedRows, numCopiedColumns);
  if (randomize)
    fputs("Random row and column permutations were applied.\n", stderr);

  /* Create permutations. */
  size_t* rowPermutation = NULL;
  TU_CALL( TUallocBlockArray(tu, &rowPermutation, numTotalRows) );
  for (size_t row = 0; row < numTotalRows; ++row)
    rowPermutation[row] = row;
  size_t* columnPermutation = NULL;
  TU_CALL( TUallocBlockArray(tu, &columnPermutation, numTotalColumns) );
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

  /* Print matrix. */

  printf("%ld %ld %ld\n\n", numTotalRows, numTotalColumns, numTotalNonzeros);
  for (size_t row = 0; row < numTotalRows; ++row)
  {
    for (Nonzero* nz = rowHeads[row].right; nz->column != SIZE_MAX; nz = nz->right)
    {
      printf("%ld %ld 1\n", rowPermutation[nz->row]+1, columnPermutation[nz->column]+1);
    }
  }

  /* Cleanup. */

  TU_CALL( TUfreeBlockArray(tu, &columnPermutation) );
  TU_CALL( TUfreeBlockArray(tu, &rowPermutation) );
  TU_CALL( TUfreeBlockArray(tu, &operations) );
  for (size_t row = 0; row < numTotalRows; ++row)
  {
    for (Nonzero* nz = rowHeads[row].right; nz != &rowHeads[row]; )
    {
      Nonzero* current = nz;
      nz = nz->right;
      TU_CALL( TUfreeBlock(tu, &current) );
    }
  }
  TU_CALL( TUfreeBlockArray(tu, &columnHeads) );
  TU_CALL( TUfreeBlockArray(tu, &rowHeads) );

  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
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
  double probability = -1.0;
  double sparsity = -1.0;
  bool randomize = false;
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
      fprintf(stderr, "Error: more than two size indicators specified: %ld %ld %s\n\n", numBaseRows, numBaseColumns,
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
  else if (probability < 0.0)
    probability = 0.5;

  TU_ERROR error = genMatrixSeriesParallel(numBaseRows, numBaseColumns, numZeroRows, numZeroColumns, numUnitRows,
    numUnitColumns, numCopiedRows, numCopiedColumns, probability, randomize);
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
