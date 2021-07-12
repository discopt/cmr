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
  puts("Applies all possible series-parallel reductions to the ternary or binary matrix in FILE.");
  puts("Options:");
  puts("  -i FORMAT  Format of input FILE; default: `dense'.");
  puts("  -o FORMAT  Format of output matrices; default: `dense'.");
  puts("  -sp        Output the list of series-parallel reductions.");
  puts("  -r         Output the elements of the reduced matrix.");
  puts("  -R         Output the reduced matrix.");
  puts("  -w         Output the elements of a wheel matrix if not series-parallel.");
  puts("  -W         Output a wheel matrix if not series-parallel.");
  puts("Matrix formats: dense, sparse");
  puts("If FILE is `-', then the input will be read from stdin.");
  return EXIT_FAILURE;
}

TU_ERROR matrixSeriesParallel2Sums(
  const char* instanceFileName, /**< File name of instance. */
  FileFormat inputFormat,       /**< Format of input matrix. */
  FileFormat outputFormat,      /**< Format of output matrices. */
  bool outputReductions,        /**< Whether to output the list of series-parallel reductions. */
  bool outputReducedElements,   /**< Whether to output the elements of the reduced matrix. */
  bool outputReducedMatrix,     /**< Whether to output the reduced matrix. */
  bool outputWheelElements,     /**< Whether to output the elements of a wheel matrix if not series-parallel. */
  bool outputWheelMatrix        /**< Whether to output a wheel matrix if not series-parallel. */
)
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
  fprintf(stderr, "Read %dx%d matrix with %d nonzeros in %f seconds.\n", matrix->numRows, matrix->numColumns,
    matrix->numNonzeros, (endTime - startTime) * 1.0 / CLOCKS_PER_SEC);

  /* Run the search. */

  TU_SP* operations = NULL;
  size_t numOperations = 0;
  TU_CALL( TUallocBlockArray(tu, &operations, matrix->numRows + matrix->numColumns) );
  TU_SUBMAT* reducedSubmatrix = NULL;
  TU_SUBMAT* wheelSubmatrix = NULL;

  startTime = clock();
  TU_CALL( TUfindSeriesParallel(tu, matrix, operations, &numOperations,
    (outputReducedElements || outputReducedMatrix) ? &reducedSubmatrix : NULL,
    (outputWheelElements || outputWheelMatrix) ? &wheelSubmatrix : NULL, NULL, NULL, true) );
  endTime = clock();

  fprintf(stderr, "Recognition done in %fs seconds. Matrix %s series-parallel.\n",
    (endTime - startTime) * 1.0 / CLOCKS_PER_SEC,
    numOperations == matrix->numRows + matrix->numColumns ? "IS" : "is NOT");

  if (outputReductions)
  {
    fprintf(stderr, "Printing %ld series-parallel reductions.\n", numOperations);
    printf("%ld\n", numOperations);
    for (size_t i = 0; i < numOperations; ++i)
      printf("%s\n", TUspString(operations[i], NULL));
  }

  if (outputReducedElements)
  {
    fprintf(stderr, "\nReduced submatrix consists of these elements:\n");
    printf("%ld rows:", reducedSubmatrix->numRows);
    for (size_t r = 0; r < reducedSubmatrix->numRows; ++r)
      printf(" %ld", reducedSubmatrix->rows[r]+1);
    printf("\n%ld columns: ", reducedSubmatrix->numColumns);
    for (size_t c = 0; c < reducedSubmatrix->numColumns; ++c)
      printf(" %ld", reducedSubmatrix->columns[c]+1);
    printf("\n");
  }

  if (outputReducedMatrix)
  {
    startTime = clock();
    TU_CHRMAT* reducedMatrix = NULL;
    TU_CALL( TUchrmatFilterSubmat(tu, matrix, reducedSubmatrix, &reducedMatrix) );
    endTime = clock();
    fprintf(stderr, "\nExtracted reduced %dx%d matrix with %d nonzeros in %f seconds.\n", reducedMatrix->numRows,
      reducedMatrix->numColumns, reducedMatrix->numNonzeros, (endTime - startTime) * 1.0 / CLOCKS_PER_SEC );
    if (outputFormat == FILEFORMAT_MATRIX_DENSE)
      TU_CALL( TUchrmatPrintDense(tu, stdout, reducedMatrix, '0', false) );
    else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
      TU_CALL( TUchrmatPrintSparse(stdout, reducedMatrix) );
    TU_CALL( TUchrmatFree(tu, &reducedMatrix) );
  }

  if (wheelSubmatrix && outputWheelElements)
  {
    fprintf(stderr, "\nWheel submatrix of order %ld consists of these elements of the input matrix:\n",
      wheelSubmatrix->numRows);
    printf("%ld rows:", wheelSubmatrix->numRows);
    for (size_t r = 0; r < wheelSubmatrix->numRows; ++r)
      printf(" %ld", wheelSubmatrix->rows[r]+1);
    printf("\n%ld columns: ", wheelSubmatrix->numColumns);
    for (size_t c = 0; c < wheelSubmatrix->numColumns; ++c)
      printf(" %ld", wheelSubmatrix->columns[c]+1);
    printf("\n");
  }

  if (wheelSubmatrix && outputWheelMatrix)
  {
    startTime = clock();
    TU_CHRMAT* wheelMatrix = NULL;
    TU_CALL( TUchrmatFilterSubmat(tu, matrix, wheelSubmatrix, &wheelMatrix) );
    endTime = clock();
    fprintf(stderr, "\nExtracted %dx%d wheel matrix with %d nonzeros in %f seconds.\n", wheelMatrix->numRows,
      wheelMatrix->numColumns, wheelMatrix->numNonzeros, (endTime - startTime) * 1.0 / CLOCKS_PER_SEC );
    if (outputFormat == FILEFORMAT_MATRIX_DENSE)
      TU_CALL( TUchrmatPrintDense(tu, stdout, wheelMatrix, '0', false) );
    else if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
      TU_CALL( TUchrmatPrintSparse(stdout, wheelMatrix) );
    TU_CALL( TUchrmatFree(tu, &wheelMatrix) );
  }
  
  /* Cleanup. */

  TU_CALL( TUsubmatFree(tu, &wheelSubmatrix) );
  TU_CALL( TUsubmatFree(tu, &reducedSubmatrix) );
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
  bool outputReductions = false;
  bool outputReducedElements = false;
  bool outputReducedMatrix = false;
  bool outputWheelElements = false;
  bool outputWheelMatrix = false;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-sp"))
      outputReductions = true;
    else if (!strcmp(argv[a], "-r"))
      outputReducedElements = true;
    else if (!strcmp(argv[a], "-R"))
      outputReducedMatrix = true;
    else if (!strcmp(argv[a], "-w"))
      outputWheelElements = true;
    else if (!strcmp(argv[a], "-W"))
      outputWheelMatrix = true;
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

  TU_ERROR error = matrixSeriesParallel2Sums(instanceFileName, inputFormat, outputFormat, outputReductions,
    outputReducedElements, outputReducedMatrix, outputWheelElements, outputWheelMatrix);
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
