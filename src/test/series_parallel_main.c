#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <tu/matrix.h>
#include <tu/series_parallel.h>

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,
  FILEFORMAT_MATRIX_DENSE = 1,
  FILEFORMAT_MATRIX_SPARSE = 2
} FileFormat;

int printUsage(const char* program)
{
  printf("Usage: %s [OPTION]... FILE\n\n", program);
  puts("Detects all series/parallel 2-sums of the ternary or binary matrix and outputs the remaining submatrix.");
  puts("Options:");
  puts("  -i FORMAT  Format of input FILE; default: `dense'.");
  puts("  -o FORMAT  Format of output; default: `dense'.");
  puts("Formats: dense, sparse");
  puts("If FILE is `-', then the input will be read from stdin.");
  return EXIT_FAILURE;
}

TU_ERROR matrixSeriesParallel2Sums(const char* instanceFileName, FileFormat inputFormat, FileFormat outputFormat)
{
  clock_t startTime, endTime;
  FILE* instanceFile = strcmp(instanceFileName, "-") ? fopen(instanceFileName, "r") : stdin;
  if (!instanceFile)
    return TU_ERROR_INPUT;

  TU* tu = NULL;
  TU_CALL( TUcreateEnvironment(&tu) );

  /* Read matrix. */

  startTime = clock();
  TU_CHRMAT* matrix = NULL;
  if (inputFormat == FILEFORMAT_MATRIX_DENSE)
    TU_CALL( TUchrmatCreateFromDenseStream(tu, &matrix, instanceFile) );
  else if (inputFormat == FILEFORMAT_MATRIX_SPARSE)
    TU_CALL( TUchrmatCreateFromSparseStream(tu, &matrix, instanceFile) );
  if (instanceFile != stdin)
    fclose(instanceFile);
  endTime = clock();
  fprintf(stderr, "Read %dx%d matrix with %d nonzeros: %fs\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (endTime - startTime) * 1.0 / CLOCKS_PER_SEC);

  /* Run the search. */

  TU_SERIES_PARALLEL* operations = NULL;
  size_t numRemoved = 0;
  TU_CALL( TUallocBlockArray(tu, &operations, matrix->numRows + matrix->numColumns) );
  TU_SUBMAT* remainingSubmatrix = NULL;
//   TU_SUBMAT* wheelSubmatrix = NULL;

  startTime = clock();
  TU_CALL( TUfindSeriesParallel(tu, matrix, operations, &numRemoved, &remainingSubmatrix, /*&wheelSubmatrix, */true) );
  endTime = clock();

  fprintf(stderr, "Removed %ld series/parallel/zero elements: %fs\n", numRemoved,
    (endTime - startTime) * 1.0 / CLOCKS_PER_SEC);

  startTime = clock();
  TU_CHRMAT* remainingMatrix = NULL;
  if (remainingSubmatrix)
    TU_CALL( TUchrmatFilterSubmat(tu, matrix, remainingSubmatrix, &remainingMatrix) );
  endTime = clock();
  fprintf(stderr, "Extracted %dx%d matrix with %d nonzeros: %fs\n", remainingMatrix->numRows,
    remainingMatrix->numColumns, remainingMatrix->numNonzeros, (endTime - startTime) * 1.0 / CLOCKS_PER_SEC);

  /* Cleanup. */

  TU_CALL( TUchrmatFree(tu, &remainingMatrix) );
//   TU_CALL( TUsubmatFree(tu, &wheelSubmatrix) );
  TU_CALL( TUsubmatFree(tu, &remainingSubmatrix) );
  TU_CALL( TUfreeBlockArray(tu, &operations) );
  TU_CALL( TUchrmatFree(tu, &matrix) );
  TU_CALL( TUfreeEnvironment(&tu) );

  return TU_OKAY;
}

int main(int argc, char** argv)
{
  FileFormat inputFormat = FILEFORMAT_MATRIX_DENSE;
  FileFormat outputFormat = FILEFORMAT_MATRIX_DENSE;
  char* instanceFileName = NULL;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }

    else if (!strcmp(argv[a], "-i") && a+1 < argc)
    {
      if (!strcmp(argv[a+1], "dense"))
        inputFormat = FILEFORMAT_MATRIX_DENSE;
      else if (!strcmp(argv[a+1], "sparse"))
        inputFormat = FILEFORMAT_MATRIX_SPARSE;
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

  TU_ERROR error = matrixSeriesParallel2Sums(instanceFileName, inputFormat, outputFormat);
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
