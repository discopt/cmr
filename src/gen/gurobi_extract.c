#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <stdint.h>
#include <math.h>

#include <cmr/graphic.h>

#include <gurobi_c.h>

typedef enum
{
  FILEFORMAT_UNDEFINED = 0,       /**< Whether the file format of output was defined by the user. */
  FILEFORMAT_MATRIX_DENSE = 1,    /**< Dense matrix format. */
  FILEFORMAT_MATRIX_SPARSE = 2,   /**< Sparse matrix format. */
} FileFormat;

#define GRB_CALL(x) \
  do \
  { \
    int _grb_error = x; \
    if (_grb_error) \
    { \
      fprintf(stderr, "[Gurobi error <%d>: %s]", _grb_error, GRBgeterrormsg(env)); \
    } \
  } while(0)

static
CMR_ERROR extractMatrixGurobi(
  char* mipFileName,        /**< File name of MIP. */
  FileFormat outputFormat   /**< Output file format. */
  )
{
  GRBenv* env = NULL;
  GRBemptyenv(&env);
  GRB_CALL( GRBsetintparam(env, "outputflag", 0) );
  GRBstartenv(env);

  GRBmodel* model = NULL;
  GRB_CALL( GRBreadmodel(env, mipFileName, &model) );

  int numRows, numColumns, numNonzeros;
  GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_NUMCONSTRS, &numRows) );
  GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_NUMVARS, &numColumns) );
  GRB_CALL( GRBgetintattr(model, GRB_INT_ATTR_NUMNZS, &numNonzeros) );

  CMR* cmr = NULL;
  CMR_CALL( CMRcreateEnvironment(&cmr) );
  
  CMR_DBLMAT* matrix = NULL;
  CMR_CALL( CMRdblmatCreate(cmr, &matrix, numRows, numColumns, numNonzeros) );
  size_t entry = 0;
  
  int* begin = NULL;
  int* indices = NULL;
  double* values = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &begin, numRows+1) );
  CMR_CALL( CMRallocBlockArray(cmr, &indices, numNonzeros) );
  CMR_CALL( CMRallocBlockArray(cmr, &values, numNonzeros) );
  
  GRB_CALL( GRBgetconstrs(model, &numNonzeros, begin, indices, values, 0, numRows) );
  for (size_t row = 0; row < numRows; ++row)
  {
    matrix->rowSlice[row] = entry;
    size_t beyond = (row+1 < numRows) ? begin[row+1] : numNonzeros;
    for (size_t i = begin[row]; i < beyond; ++i)
    {
      matrix->entryColumns[entry] = indices[i];
      matrix->entryValues[entry] = values[i];
      ++entry;
    }
  }
  matrix->rowSlice[numRows] = entry;
  assert(entry == numNonzeros);
  
  if (outputFormat == FILEFORMAT_MATRIX_SPARSE)
    CMR_CALL( CMRdblmatPrintSparse(cmr, matrix, stdout) );
  else if (outputFormat == FILEFORMAT_MATRIX_DENSE)
    CMR_CALL( CMRdblmatPrintDense(cmr, matrix, stdout, '0', false) );

  CMR_CALL( CMRfreeBlockArray(cmr, &values) );
  CMR_CALL( CMRfreeBlockArray(cmr, &indices) );
  CMR_CALL( CMRfreeBlockArray(cmr, &begin) );
  CMR_CALL( CMRdblmatFree(cmr, &matrix) );
  CMR_CALL( CMRfreeEnvironment(&cmr) );
  
  GRBfreemodel(model);

  GRBfreeenv(env);
  
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

  fprintf(stderr, "%s MIPFILE [OPTION]...\n", program);
  fputs("  extracts the coefficient matrix of a mixed-integer program in MIPFILE and writes it to stdout.\n", stderr);
  fputs("\n", stderr);

  fputs("Options:\n", stderr);
  fputs("  -o FORMAT  Format of output matrix; default: dense.\n", stderr);
  fputs("\n", stderr);

  fputs("Formats for matrices: dense, sparse\n", stderr);
  fputs("MIPFILE must refer to a file that Gurobi can read.\n", stderr);

  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  FileFormat outputFormat = FILEFORMAT_MATRIX_DENSE;
  char* mipFileName = NULL;
  for (int a = 1; a < argc; ++a)
  {
    if (!strcmp(argv[a], "-h"))
    {
      printUsage(argv[0]);
      return EXIT_SUCCESS;
    }
    else if (!strcmp(argv[a], "-o") && (a+1 < argc))
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
    else if (!mipFileName)
    {
      mipFileName = argv[a];
    }
    else
    {
      printf("Error: more than two MIP file names: %s %s\n\n", mipFileName, argv[a]);
      return printUsage(argv[0]);
    }
  }

  if (!mipFileName)
  {
    puts("Error: no MIP file name specified.\n");
    return printUsage(argv[0]);
  }

  CMR_ERROR error = extractMatrixGurobi(mipFileName, outputFormat);
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
