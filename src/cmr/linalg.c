// #define CMR_DEBUG /* Uncomment to debug this file. */

#include "linalg.h"
#include "listmatrix.h"
#include "sort.h"

#include "env_internal.h"
#include <assert.h>
#include <stdint.h>
#include <limits.h>

/**
 * TODO: Implement (transposed) HNF according to
 *
 * `Polynomial Algorithms for Computing the Smith and Hermite Normal Forms of an Integer Matrix'
 *
 **/

/**
 * \brief Computes the gcd of \p a and \p b along with the Bezout coefficients.
 *
 * Runs the extended Euclidean algorithm. For the gcd \f$ g \f$, we will have \f$ g = s a + t b. \f$.
 * We are guaranteed to have \f$ t \neq 0 \f$.
 *
 */

static
int64_t gcdExt(
  int64_t a,    /**< First number . */
  int64_t b,    /**< Second number. */
  int64_t* ps,  /**< Pointer for storing the first Bezout coefficient. */
  int64_t* pt   /**< Pointer for storing the second Bezout coefficient. */
)
{
  int64_t s = 0;
  int64_t old_s = 1;
  int64_t t = 1;
  int64_t old_t = 0;
  int64_t r = b;
  int64_t old_r = a;

  while(r != 0)
  {
    int64_t q = old_r / r;
    int64_t temp = r;
    r = old_r - q * r;
    old_r = temp;

    temp = s;
    s = old_s - q * s;
    old_s = temp;

    temp = t;
    t = old_t - q * t;
    old_t = temp;
  }

  if (old_r > 0)
  {
    *ps = old_s;
    *pt = old_t;
    return old_r;
  }
  else
  {
    *ps = -old_s;
    *pt = -old_t;
    return -old_r;
  }
}

typedef struct _RowInfo64
{
  size_t row;
  int64_t value;
  long priority;
} RowInfo64;

static
int compareOtherRows64(const void* first, const void* second)
{
  const RowInfo64* row1 = (const RowInfo64*) first;
  const RowInfo64* row2 = (const RowInfo64*) second;
  return row1->priority - row2->priority;
}


#if defined(CMR_WITH_GMP)

typedef struct _RowInfoGMP
{
  size_t row;
  mpz_t value;
  long priority;
} RowInfoGMP;

static
int compareOtherRowsGMP(const void* first, const void* second)
{
  const RowInfoGMP* row1 = (const RowInfoGMP*) first;
  const RowInfoGMP* row2 = (const RowInfoGMP*) second;
  return row1->priority - row2->priority;
}

static CMR_ERROR CMRintmatComputeUpperDiagonalGMP(CMR* cmr, CMR_INTMAT* matrix, bool invert, size_t* prank,
  CMR_SUBMAT* permutations, CMR_INTMAT** presult, CMR_INTMAT** ptranspose)
{
  assert(cmr);
  assert(matrix);
  assert(permutations);
  assert(prank);

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "CMRintmatComputeUpperDiagonalGMP for matrix.\n");
  CMRintmatPrintDense(cmr, matrix, stdout, '0', true);
#endif /* CMR_DEBUG */

  for (size_t row = 0; row < matrix->numRows; ++row)
    permutations->rows[row] = row;
  for (size_t column = 0; column < matrix->numColumns; ++column)
    permutations->columns[column] = column;

  ListMatGMP* listmatrix = NULL;
  CMR_CALL( CMRlistmatGMPAlloc(cmr, matrix->numRows, matrix->numColumns,
    2 * matrix->numNonzeros + matrix->numRows + matrix->numColumns, &listmatrix) );
  CMR_CALL( CMRlistmatGMPInitializeFromIntMatrix(cmr, listmatrix, matrix) );

  size_t maxRank = matrix->numRows < matrix->numColumns ? matrix->numRows : matrix->numColumns;
  size_t maxNumRowsColumns = matrix->numRows > matrix->numColumns ? matrix->numRows : matrix->numColumns;
  mpz_t* densePivot = NULL; /* Dense representation of the pivot row / column. */
  CMR_CALL( CMRallocStackArray(cmr, &densePivot, maxNumRowsColumns) );
  bool* denseProcessed = NULL; /* Indicator whether a row/column appeared in pivot and other. */
  CMR_CALL( CMRallocStackArray(cmr, &denseProcessed, maxNumRowsColumns) );
  size_t numOtherRows;
  RowInfoGMP* otherRowInfos = NULL; /* Array with rows that have a nonzero in the pivot column. */
  CMR_CALL( CMRallocStackArray(cmr, &otherRowInfos, matrix->numRows) );
  size_t* originalRowsToPermutedRows = NULL; /* Array that maps permuted rows to rows. */
  CMR_CALL( CMRallocStackArray(cmr, &originalRowsToPermutedRows, matrix->numRows) );
  size_t* originalColumnsToPermutedColumns = NULL; /* Array that maps permuted columns to columns. */
  CMR_CALL( CMRallocStackArray(cmr, &originalColumnsToPermutedColumns, matrix->numColumns) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    originalRowsToPermutedRows[row] = row;
  for (size_t column = 0; column < matrix->numColumns; ++column)
    originalColumnsToPermutedColumns[column] = column;

  for (size_t e = 0; e < maxNumRowsColumns; ++e)
  {
    mpz_init(densePivot[e]);
    denseProcessed[e] = false;
  }

  *prank = 0;
  mpz_t minEntryValue;
  mpz_t minEntryAbs;
  mpz_t pivotValue;
  mpz_init(minEntryValue);
  mpz_init(minEntryAbs);
  mpz_init(pivotValue);
  while (*prank < maxRank)
  {
    /* TODO: maintain smallest elements instead of searching for them. */

    mpz_set_si(minEntryAbs, 0);
    size_t minEntryRow = SIZE_MAX;
    size_t minEntryColumn = SIZE_MAX;
    size_t minEntryArea = SIZE_MAX;
    for (size_t permutedRow = *prank; permutedRow < matrix->numRows; ++permutedRow)
    {
      size_t row = permutations->rows[permutedRow];
      CMRdbgMsg(4, "Scanning permuted row %ld (row %ld) for good pivot entries.\n", permutedRow, row);
      for (ListMatGMPNonzero* nz = listmatrix->rowElements[row].head.right; nz != &listmatrix->rowElements[row].head;
        nz = nz->right)
      {
        mpz_t x;
        mpz_init(x);
        mpz_abs(x, nz->value);
        int comparison = mpz_cmp(x, minEntryAbs);
        if (mpz_sgn(minEntryAbs) > 0 && comparison > 0)
        {
          mpz_clear(x);
          continue;
        }

        size_t area = (listmatrix->rowElements[row].numNonzeros - 1)
          * (listmatrix->columnElements[nz->column].numNonzeros - 1);
        if (mpz_sgn(minEntryAbs) <= 0 || comparison < 0 || (comparison == 0 && area < minEntryArea))
        {
          minEntryRow = row;
          minEntryColumn = nz->column;
          minEntryArea = area;
          mpz_set(minEntryValue, nz->value);
          mpz_abs(minEntryAbs, minEntryValue);
        }
        mpz_clear(x);
      }
    }

#if defined(CMR_DEBUG)
    char buffer[1024];
    CMRdbgMsg(2, "Pivot row %d, column %d has min nonzero %s and fill-in area %ld.\n", minEntryRow, minEntryColumn,
      mpz_get_str(buffer, 10, minEntryValue), minEntryArea, minEntryArea);
    CMR_CALL( CMRlistmatGMPPrintDense(cmr, listmatrix, stdout) );
#endif /* CMR_DEBUG */

    if (minEntryRow == SIZE_MAX)
      break;

    /* Find the permuted row and column. */
    mpz_set(pivotValue, minEntryValue);
    size_t pivotRow = minEntryRow;
    size_t pivotColumn = minEntryColumn;
    size_t pivotPermutedRow = originalRowsToPermutedRows[pivotRow];
    size_t pivotPermutedColumn = originalColumnsToPermutedColumns[pivotColumn];

    /* Update permutations such that rowPermutation[*prank],columnPermutation[*prank] denotes the pivot entry. */
    CMRdbgMsg(2, "Pivoting at %ld,%ld\nResulting permutation is:\n", pivotRow, pivotColumn);
    size_t swap;
    swap = permutations->rows[*prank];
    permutations->rows[*prank] = permutations->rows[pivotPermutedRow];
    permutations->rows[pivotPermutedRow] = swap;
    originalRowsToPermutedRows[permutations->rows[*prank]] = *prank;
    originalRowsToPermutedRows[permutations->rows[pivotPermutedRow]] = pivotPermutedRow;

    swap = permutations->columns[*prank];
    permutations->columns[*prank] = permutations->columns[pivotPermutedColumn];
    permutations->columns[pivotPermutedColumn] = swap;
    originalColumnsToPermutedColumns[permutations->columns[*prank]] = *prank;
    originalColumnsToPermutedColumns[permutations->columns[pivotPermutedColumn]] = pivotPermutedColumn;

#if defined(CMR_DEBUG)
    CMRsubmatWriteToStream(cmr, permutations, matrix->numRows, matrix->numColumns, stdout);
#endif /* CMR_DEBUG */

    /* Go through the nonzeros in the pivot column and sort them to prioritize. */
    numOtherRows = 0;
    bool scalePivotRow = mpz_sgn(pivotValue) < 0;
    for (ListMatGMPNonzero* nz = listmatrix->columnElements[pivotColumn].head.below;
      nz != &listmatrix->columnElements[pivotColumn].head; nz = nz->below)
    {
      if (nz->row == pivotRow)
      {
        assert(originalRowsToPermutedRows[nz->row] == *prank);
      }
      else if (originalRowsToPermutedRows[nz->row] > *prank)
      {
        mpz_t s, t;
        mpz_init(s);
        mpz_init(t);
        mpz_gcdext(NULL, s, t, nz->value, pivotValue);
        if (mpz_sgn(s) == 0)
          otherRowInfos[numOtherRows].priority = 0; /* Highest priority since divisible by pivot value. */
        else
          otherRowInfos[numOtherRows].priority = listmatrix->rowElements[nz->row].numNonzeros;
        otherRowInfos[numOtherRows].row = nz->row;
        mpz_init_set(otherRowInfos[numOtherRows].value, nz->value);
        ++numOtherRows;
        scalePivotRow = false;
        mpz_clear(t);
        mpz_clear(s);
      }
      else if (invert)
      {
        otherRowInfos[numOtherRows].row = nz->row;
        mpz_init_set(otherRowInfos[numOtherRows].value, nz->value);
        otherRowInfos[numOtherRows].priority = INT32_MAX; /* Lowest priority for top rows. */
        ++numOtherRows;
      }
    }

    /* We can now increase the rank. */
    ++(*prank);

    /* If the pivot entry is negative and have no other rows to combine with (except for upper-diagonal part) then
     * we just scale the entire row. */
    if (scalePivotRow)
    {
      for (ListMatGMPNonzero* nz = listmatrix->rowElements[pivotRow].head.right;
        nz != &listmatrix->rowElements[pivotRow].head; nz = nz->right)
      {
        mpz_neg(nz->value, nz->value);
      }
    }

    /* Continue if no updates are needed. */
    if (numOtherRows == 0)
    {
      continue;
    }

    /* Sort other rows. */
    CMR_CALL( CMRsort(cmr, numOtherRows, otherRowInfos, sizeof(RowInfoGMP), compareOtherRowsGMP) );

    /* Copy pivot row to densePivot array. */
    for (ListMatGMPNonzero* nz = listmatrix->rowElements[pivotRow].head.right;
      nz != &listmatrix->rowElements[pivotRow].head; nz = nz->right)
    {
      mpz_set(densePivot[nz->column], nz->value);
      CMRdbgMsg(4, "Copying %ld into densePivot[%ld]\n", nz->value, nz->column);
    }

    CMRdbgMsg(2, "Processing all other rows.\n");

    /* Process every other row. */
    for (size_t i = 0; i < numOtherRows; ++i)
    {
      mpz_set(pivotValue, densePivot[pivotColumn]);
      size_t otherRow = otherRowInfos[i].row;
      char buffer[1024];
      CMRdbgMsg(4, "Other row %ld has value %s, priority %lld and %ld nonzeros.\n", otherRowInfos[i].row,
        mpz_get_str(buffer, 10, otherRowInfos[i].value), otherRowInfos[i].priority, listmatrix->rowElements[otherRowInfos[i].row].numNonzeros);

      mpz_t gcd;
      mpz_t U_11;
      mpz_t U_12;
      mpz_t U_21;
      mpz_t U_22;
      mpz_init(gcd);
      mpz_init(U_11);
      mpz_init(U_12);
      mpz_init(U_21);
      mpz_init(U_22);
      if (otherRowInfos[i].priority < INT32_MAX)
      {
        /* Rows below the diagonal will have a zero in the pivot column. */
        mpz_gcdext(gcd, U_12, U_11, otherRowInfos[i].value, pivotValue);

#if defined(CMR_DEBUG)
        char buffer1[1024];
        char buffer2[1024];
        char buffer3[1024];
        char buffer4[1024];
        char buffer5[1024];
        char buffer6[1024];
        char buffer7[1024];
        char buffer8[1024];
        CMRdbgMsg(6, "gcd is %s, s = %s, t = %s. Bezout: %s = %s * %s + %s * %s\n",
          mpz_get_str(buffer1, 10, gcd), mpz_get_str(buffer2, 10, U_11), mpz_get_str(buffer3, 10, U_12),
          mpz_get_str(buffer4, 10, gcd), mpz_get_str(buffer5, 10, U_11), mpz_get_str(buffer6, 10, pivotValue),
          mpz_get_str(buffer7, 10, U_12), mpz_get_str(buffer8, 10, otherRowInfos[i].value));
#endif /* CMR_DEBUG */

        mpz_div(U_21, otherRowInfos[i].value, gcd);
        mpz_neg(U_21, U_21);
        mpz_div(U_22, pivotValue, gcd);
      }
      else
      {
        /* For rows above the diagonal we cannot risk to modify the pivot row, but we still subtract an integer
         * multiple. */
        mpz_set_si(U_11, 1);
        mpz_set_si(U_12, 0);
        mpz_fdiv_q(U_21, otherRowInfos[i].value, pivotValue);
        mpz_neg(U_21, U_21);
        mpz_set_si(U_22, 1);
      }


#if defined(CMR_DEBUG)
      char buffer1[1024];
      char buffer2[1024];
      char buffer3[1024];
      char buffer4[1024];
      CMRdbgMsg(6, "U = [[%s, %s], [%s, %s]]\n", mpz_get_str(buffer1, 10, U_11), mpz_get_str(buffer2, 10, U_12),
        mpz_get_str(buffer3, 10, U_21), mpz_get_str(buffer4, 10, U_22));
#endif /* CMR_DEBUG */

      /*
       * Unimodular matrix U applied to the rows (from the left):
       *
       *   [   s  ,  t  ]
       *   [ -b/g , a/g ]
       *
       * Its determinant is s*a/g - t*(-b)/g = (s*a + t*b)/g = g/g = 1.
       *
       * The pivot column will then have the updated entries g (pivot row) and (-b)/g*a + a/g*b = 0 (other row).
       */

      /* Go through this other row. We add new nonzeros to the linked lists, but do not remove new zeros.
       * Also the values in the linked list are not updated. */
      for (ListMatGMPNonzero* iter = listmatrix->rowElements[otherRow].head.right;
        iter != &listmatrix->rowElements[otherRow].head; )
      {
        ListMatGMPNonzero* nz = iter;
        iter = iter->right;
        size_t column = nz->column;
        denseProcessed[column] = true;

        /* Apply unimodular transformation to both entries. */
        mpz_t p_old;
        mpz_init_set(p_old, densePivot[column]);
        /* o_old is nz->value */
        mpz_t temp, p_new, o_new;
        mpz_init(temp);
        mpz_init(p_new);
        mpz_init(o_new);
        mpz_mul(p_new, U_11, p_old);
        mpz_mul(temp, U_12, nz->value);
        mpz_add(p_new, p_new, temp);
        mpz_mul(o_new, U_21, p_old);
        mpz_mul(temp, U_22, nz->value);
        mpz_add(o_new, o_new, temp);

#if defined(CMR_DEBUG)
        CMRdbgMsg(8, "Other row's nonzero %s vs. pivot %s -> %s and %s; all in column %ld.\n",
          mpz_get_str(buffer1, 10, nz->value), mpz_get_str(buffer2, 10, p_old), mpz_get_str(buffer3, 10, o_new),
          mpz_get_str(buffer4, 10, p_new), nz->column);
#endif /* CMR_DEBUG */

        mpz_set(densePivot[column], p_new);
        mpz_set(nz->value, o_new);
        if (mpz_sgn(p_old) == 0 && mpz_sgn(p_new) != 0)
        {
          /* A new nonzero in the pivot row is added as a 1 because we copy back from the dense vector anyhow. */
          ptrdiff_t memoryShift = 0;
          CMRdbgMsg(8, "Inserting a dummy 1 into pivot row %d in column %d.\n", pivotRow, column);
          mpz_t temp;
          mpz_init_set_si(temp, 1);
          CMR_CALL( CMRlistmatGMPInsert(cmr, listmatrix, pivotRow, column, temp, 0, &memoryShift) );
          mpz_clear(temp);
          if (memoryShift)
            iter += memoryShift;
        }

        mpz_clear(p_old);
        mpz_clear(temp);
        mpz_clear(p_new);
        mpz_clear(o_new);
      }

      /* Go through pivot row. */
      for (ListMatGMPNonzero* nz = listmatrix->rowElements[pivotRow].head.right;
        nz != &listmatrix->rowElements[pivotRow].head; nz = nz->right)
      {
        size_t column = nz->column;
        if (denseProcessed[column])
          continue;

        /* Apply unimodular transformation to both entries. */
        mpz_t p_old, p_new, o_new;
        mpz_init(p_old);
        mpz_init(p_new);
        mpz_init(o_new);
        mpz_set(p_old, densePivot[column]);
        mpz_mul(p_new, U_11, p_old);
        mpz_mul(o_new, U_21, p_old);

        mpz_set(densePivot[column], p_new);
        ptrdiff_t memoryShift;
        CMRdbgMsg(8, "Inserting into other row %ld in column %ld. Value is %ld\n", otherRow, column, o_new);
        CMR_CALL( CMRlistmatGMPInsert(cmr, listmatrix, otherRow, column, o_new, 0, &memoryShift) );
        if (memoryShift)
          nz += memoryShift;

        mpz_clear(p_old);
        mpz_clear(p_new);
        mpz_clear(o_new);
      }

      /* Go through other row again to remove nonzeros with value 0. */
      for (ListMatGMPNonzero* iter = listmatrix->rowElements[otherRow].head.right;
        iter != &listmatrix->rowElements[otherRow].head; )
      {
        ListMatGMPNonzero* nz = iter;
        iter = iter->right;
        denseProcessed[nz->column] = false;
        CMRdbgMsg(8, "Deselecting nonzero at %ld,%ld\n", nz->row, nz->column);
        if (mpz_sgn(nz->value) == 0)
        {
          CMRdbgMsg(8, "Removing nonzero at %ld,%ld\n", nz->row, nz->column);
          CMR_CALL( CMRlistmatGMPDelete(cmr, listmatrix, nz) );
        }
      }

      mpz_clear(gcd);
      mpz_clear(U_11);
      mpz_clear(U_12);
      mpz_clear(U_21);
      mpz_clear(U_22);
    }
    for (size_t i = 0; i < numOtherRows; ++i)
      mpz_clear(otherRowInfos[i].value);

    /* Go through pivot row and update the linked list data. */
    for (ListMatGMPNonzero* iter = listmatrix->rowElements[pivotRow].head.right;
      iter != &listmatrix->rowElements[pivotRow].head; )
    {
      ListMatGMPNonzero* nz = iter;
      iter = iter->right;
      mpz_t entry;
      mpz_init(entry);
      mpz_swap(entry, densePivot[nz->column]);

      char buffer[1024];
      CMRdbgMsg(6, "Updating entry in pivot row %ld, column %ld, value: %s\n\n", pivotRow, nz->column,
        mpz_get_str(buffer, 10, entry));
      if (mpz_sgn(entry) != 0)
        mpz_swap(nz->value, entry);
      else
        CMR_CALL( CMRlistmatGMPDelete(cmr, listmatrix, nz) );
      mpz_clear(entry);
    }
  }
  mpz_clear(pivotValue);
  mpz_clear(minEntryAbs);
  mpz_clear(minEntryValue);

  /* If requested, write the resulting matrix back into an int matrix. */
  CMR_ERROR error = CMR_OKAY;
  if (presult)
  {
    CMR_INTMAT* result = NULL;
    CMR_CALL( CMRintmatCreate(cmr, &result, listmatrix->numRows, listmatrix->numColumns, listmatrix->numNonzeros) );
    result->rowSlice[0] = 0;
    result->numNonzeros = 0;
    for (size_t permRow = 0; (permRow < matrix->numRows) && (error == CMR_OKAY); ++permRow)
    {
      size_t row = permutations->rows[permRow];
      for (ListMatGMPNonzero* nz = listmatrix->rowElements[row].head.right;
        nz != &listmatrix->rowElements[row].head; nz = nz->right)
      {
        size_t permColumn = originalColumnsToPermutedColumns[nz->column];
        if (!mpz_fits_slong_p(nz->value))
        {
          error = CMR_ERROR_OVERFLOW;
          break;
        }
        long value = mpz_get_si(nz->value);
        if (value > INT_MAX || value < INT_MIN)
        {
          error = CMR_ERROR_OVERFLOW;
          break;
        }
        result->entryValues[result->numNonzeros] = (int) value;
        result->entryColumns[result->numNonzeros] = permColumn;
        result->numNonzeros++;
      }
      result->rowSlice[permRow+1] = result->numNonzeros;
    }
    if (error == CMR_OKAY)
      CMR_CALL( CMRintmatSortNonzeros(cmr, result) );

    if (ptranspose)
      CMR_CALL( CMRintmatTranspose(cmr, result, ptranspose) );

    if (presult)
      *presult = NULL;
    if (presult && (error == CMR_OKAY))
      *presult = result;
    else
      CMR_CALL( CMRintmatFree(cmr, &result) );
  }

  CMR_CALL( CMRfreeStackArray(cmr, &originalColumnsToPermutedColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &originalRowsToPermutedRows) );
  CMR_CALL( CMRfreeStackArray(cmr, &otherRowInfos) );
  CMR_CALL( CMRfreeStackArray(cmr, &denseProcessed) );
  CMR_CALL( CMRfreeStackArray(cmr, &densePivot) );

  CMR_CALL( CMRlistmatGMPFree(cmr, &listmatrix) );

  return error;
}

#endif /* CMR_WITH_GMP */


CMR_ERROR CMRintmatComputeUpperDiagonal(CMR* cmr, CMR_INTMAT* matrix, bool invert, size_t* prank,
  CMR_SUBMAT** ppermutations, CMR_INTMAT** presult, CMR_INTMAT** ptranspose)
{
  assert(cmr);
  assert(matrix);
  assert(ppermutations);
  assert(prank);

  bool isIntTooSmall = false;

#if defined(CMR_DEBUG)
  CMRdbgMsg(0, "CMRintmatComputeUpperDiagonal for matrix.\n");
  CMRintmatPrintDense(cmr, matrix, stdout, '0', true);
#endif /* CMR_DEBUG */

  CMR_CALL( CMRsubmatCreate(cmr, matrix->numRows, matrix->numColumns, ppermutations) );
  CMR_SUBMAT* permutations = *ppermutations;
  assert(permutations);

  for (size_t row = 0; row < matrix->numRows; ++row)
    permutations->rows[row] = row;
  for (size_t column = 0; column < matrix->numColumns; ++column)
    permutations->columns[column] = column;

  ListMat64* listmatrix = NULL;
  CMR_CALL( CMRlistmat64Alloc(cmr, matrix->numRows, matrix->numColumns,
    2 * matrix->numNonzeros + matrix->numRows + matrix->numColumns, &listmatrix) );
  CMR_CALL( CMRlistmat64InitializeFromIntMatrix(cmr, listmatrix, matrix) );

  size_t maxRank = matrix->numRows < matrix->numColumns ? matrix->numRows : matrix->numColumns;
  size_t maxNumRowsColumns = matrix->numRows > matrix->numColumns ? matrix->numRows : matrix->numColumns;
  int64_t* densePivot = NULL; /* Dense representation of the pivot row / column. */
  CMR_CALL( CMRallocStackArray(cmr, &densePivot, maxNumRowsColumns) );
  bool* denseProcessed = NULL; /* Indicator whether a row/column appeared in pivot and other. */
  CMR_CALL( CMRallocStackArray(cmr, &denseProcessed, maxNumRowsColumns) );
  size_t numOtherRows;
    RowInfo64* otherRowInfos = NULL; /* Array with rows that have a nonzero in the pivot column. */
  CMR_CALL( CMRallocStackArray(cmr, &otherRowInfos, matrix->numRows) );
  size_t* originalRowsToPermutedRows = NULL; /* Array that maps permuted rows to rows. */
  CMR_CALL( CMRallocStackArray(cmr, &originalRowsToPermutedRows, matrix->numRows) );
  size_t* originalColumnsToPermutedColumns = NULL; /* Array that maps permuted columns to columns. */
  CMR_CALL( CMRallocStackArray(cmr, &originalColumnsToPermutedColumns, matrix->numColumns) );
  for (size_t row = 0; row < matrix->numRows; ++row)
    originalRowsToPermutedRows[row] = row;
  for (size_t column = 0; column < matrix->numColumns; ++column)
    originalColumnsToPermutedColumns[column] = column;

  for (size_t e = 0; e < maxNumRowsColumns; ++e)
  {
    densePivot[e] = 0;
    denseProcessed[e] = false;
  }

  *prank = 0;
  while (*prank < maxRank)
  {
    /* TODO: maintain smallest elements instead of searching for them. */

    int64_t minEntryValue = INT64_MAX;
    size_t minEntryRow = SIZE_MAX;
    size_t minEntryColumn = SIZE_MAX;
    size_t minEntryArea = SIZE_MAX;
    for (size_t permutedRow = *prank; permutedRow < matrix->numRows; ++permutedRow)
    {
      size_t row = permutations->rows[permutedRow];
      CMRdbgMsg(4, "Permuted row %ld is row %ld\n", permutedRow, row);
      for (ListMat64Nonzero* nz = listmatrix->rowElements[row].head.right; nz != &listmatrix->rowElements[row].head;
        nz = nz->right)
      {
        int64_t x = llabs(nz->value);
        if (x > llabs(minEntryValue))
          continue;

        size_t area = (listmatrix->rowElements[row].numNonzeros - 1)
          * (listmatrix->columnElements[nz->column].numNonzeros - 1);
        if (x < llabs(minEntryValue) || (x == llabs(minEntryValue) && area < minEntryArea))
        {
          minEntryRow = row;
          minEntryColumn = nz->column;
          minEntryArea = area;
          minEntryValue = nz->value;
        }
      }
    }

#if defined(CMR_DEBUG)
    CMRdbgMsg(2, "Pivot row %d, column %d has min nonzero %ld and fill-in area %ld.\n", minEntryRow, minEntryColumn,
      minEntryValue, minEntryArea, minEntryArea);
    CMR_CALL( CMRlistmat64PrintDense(cmr, listmatrix, stdout) );
#endif /* CMR_DEBUG */

    if (minEntryRow == SIZE_MAX)
      break;

    /* Find the permuted row and column. */
    int64_t pivotValue = minEntryValue;
    size_t pivotRow = minEntryRow;
    size_t pivotColumn = minEntryColumn;
    size_t pivotPermutedRow = originalRowsToPermutedRows[pivotRow];
    size_t pivotPermutedColumn = originalColumnsToPermutedColumns[pivotColumn];

    /* Update permutations such that rowPermutation[*prank],columnPermutation[*prank] denotes the pivot entry. */
    CMRdbgMsg(2, "Pivoting at %ld,%ld\nResulting permutation is:\n", pivotRow, pivotColumn);
    CMRdbgMsg(2, "pivotPermutedRow = %ld\n", pivotPermutedRow);
    CMRdbgMsg(2, "pivotPermutedColumn = %ld\n", pivotPermutedColumn);
    CMRdbgMsg(2, "rank = %ld\n", *prank);

    size_t swap;
    swap = permutations->rows[*prank];
    permutations->rows[*prank] = permutations->rows[pivotPermutedRow];
    permutations->rows[pivotPermutedRow] = swap;
    originalRowsToPermutedRows[permutations->rows[*prank]] = *prank;
    originalRowsToPermutedRows[permutations->rows[pivotPermutedRow]] = pivotPermutedRow;

    swap = permutations->columns[*prank];
    permutations->columns[*prank] = permutations->columns[pivotPermutedColumn];
    permutations->columns[pivotPermutedColumn] = swap;
    originalColumnsToPermutedColumns[permutations->columns[*prank]] = *prank;
    originalColumnsToPermutedColumns[permutations->columns[pivotPermutedColumn]] = pivotPermutedColumn;

#if defined(CMR_DEBUG)
    CMRsubmatWriteToStream(cmr, permutations, matrix->numRows, matrix->numColumns, stdout);
#endif /* CMR_DEBUG */

    /* Go through the nonzeros in the pivot column and sort them to prioritize. */
    numOtherRows = 0;
    bool scalePivotRow = pivotValue < 0;
    for (ListMat64Nonzero* nz = listmatrix->columnElements[pivotColumn].head.below;
      nz != &listmatrix->columnElements[pivotColumn].head; nz = nz->below)
    {
      if (nz->row == pivotRow)
      {
        assert(originalRowsToPermutedRows[nz->row] == *prank);
      }
      else if (originalRowsToPermutedRows[nz->row] > *prank)
      {
        int64_t s, t;
        gcdExt(nz->value, pivotValue, &s, &t);
        if (s == 0)
          otherRowInfos[numOtherRows].priority = 0; /* Highest priority since divisible by pivot value. */
        else
          otherRowInfos[numOtherRows].priority = listmatrix->rowElements[nz->row].numNonzeros;
        otherRowInfos[numOtherRows].row = nz->row;
        otherRowInfos[numOtherRows].value = nz->value;
        ++numOtherRows;
        scalePivotRow = false;
      }
      else if (invert)
      {
        otherRowInfos[numOtherRows].row = nz->row;
        otherRowInfos[numOtherRows].value = nz->value;
        otherRowInfos[numOtherRows].priority = INT32_MAX; /* Lowest priority for top rows. */
        ++numOtherRows;
      }
    }

    /* We can now increase the rank. */
    ++(*prank);

    /* If the pivot entry is negative and have no other rows to combine with (except for upper-diagonal part) then
     * we just scale the entire row. */
    if (scalePivotRow)
    {
      for (ListMat64Nonzero* nz = listmatrix->rowElements[pivotRow].head.right;
        nz != &listmatrix->rowElements[pivotRow].head; nz = nz->right)
      {
        nz->value *= -1;
      }
    }

    /* Continue if no updates are needed. */
    if (numOtherRows == 0)
    {
      continue;
    }

    /* Sort other rows. */
    CMR_CALL( CMRsort(cmr, numOtherRows, otherRowInfos, sizeof(RowInfo64), compareOtherRows64) );

    /* Copy pivot row to densePivot array. */
    for (ListMat64Nonzero* nz = listmatrix->rowElements[pivotRow].head.right;
      nz != &listmatrix->rowElements[pivotRow].head; nz = nz->right)
    {
      densePivot[nz->column] = nz->value;
      CMRdbgMsg(4, "Copying %ld into densePivot[%ld]\n", nz->value, nz->column);
    }

    CMRdbgMsg(2, "Processing all other rows.\n");

    /* Process every other row. */
    for (size_t i = 0; i < numOtherRows && !isIntTooSmall; ++i)
    {
      pivotValue = densePivot[pivotColumn];
      size_t otherRow = otherRowInfos[i].row;
      CMRdbgMsg(4, "Other row %ld has value %lld, priority %lld and %ld nonzeros.\n", otherRowInfos[i].row,
        otherRowInfos[i].value, otherRowInfos[i].priority, listmatrix->rowElements[otherRowInfos[i].row].numNonzeros);

      int64_t U_11, U_12, U_21, U_22;
      if (otherRowInfos[i].priority < INT32_MAX)
      {
        /* Rows below the diagonal will have a zero in the pivot column. */
        int64_t gcd = gcdExt(otherRowInfos[i].value, pivotValue, &U_12, &U_11);
        CMRdbgMsg(6, "gcd is %d, s = %d, t = %d. Bezout: %d = %d * %d + %d * %d\n",
          gcd, U_11, U_12, gcd, U_11, pivotValue, U_12, otherRowInfos[i].value);
        U_21 = -otherRowInfos[i].value / gcd;
        U_22 = pivotValue / gcd;
      }
      else
      {
        /* For rows above the diagonal we cannot risk to modify the pivot row, but we still subtract an integer
         * multiple. */
        U_11 = 1;
        U_12 = 0;
        U_21 = - (otherRowInfos[i].value / pivotValue);
        if (otherRowInfos[i].value % pivotValue < 0)
          ++U_21;
        U_22 = 1;
      }

      CMRdbgMsg(6, "U = [[%ld, %ld], [%ld, %ld]]\n", U_11, U_12, U_21, U_22);

      int64_t normU = llabs(U_11) + llabs(U_12) + llabs(U_21) + llabs(U_22);
      if (normU <= 0 || normU > INT32_MAX)
      {
        isIntTooSmall = true;
        break;
      }

      /*
       * Unimodular matrix U applied to the rows (from the left):
       *
       *   [   s  ,  t  ]
       *   [ -b/g , a/g ]
       *
       * Its determinant is s*a/g - t*(-b)/g = (s*a + t*b)/g = g/g = 1.
       *
       * The pivot column will then have the updated entries g (pivot row) and (-b)/g*a + a/g*b = 0 (other row).
       */

      /* Go through this other row. We add new nonzeros to the linked lists, but do not remove new zeros.
       * Also the values in the linked list are not updated. */
      for (ListMat64Nonzero* iter = listmatrix->rowElements[otherRow].head.right;
        iter != &listmatrix->rowElements[otherRow].head; )
      {
        ListMat64Nonzero* nz = iter;
        iter = iter->right;
        size_t column = nz->column;
        denseProcessed[column] = true;

        /* Apply unimodular transformation to both entries. */
        int64_t p_old = densePivot[column];
        int64_t o_old = nz->value;
        int64_t p_new = U_11 * p_old + U_12 * o_old;
        int64_t o_new = U_21 * p_old + U_22 * o_old;
        CMRdbgMsg(8, "Other row's nonzero %ld vs. pivot %ld -> %ld and %ld; all in column %ld.\n", o_old, p_old, o_new,
          p_new, nz->column);

        if (llabs(p_new) > INT32_MAX || llabs(o_new) > INT32_MAX)
        {
          isIntTooSmall = true;
          break;
        }

        densePivot[column] = p_new;
        nz->value = o_new;
        if (p_old == 0 && p_new != 0)
        {
          /* A new nonzero in the pivot row is added as a 1 because we copy back from the dense vector anyhow. */
          ptrdiff_t memoryShift = 0;
          CMRdbgMsg(8, "Inserting a dummy 1 into pivot row %d in column %d.\n", pivotRow, column);
          CMR_CALL( CMRlistmat64Insert(cmr, listmatrix, pivotRow, column, 1, 0, &memoryShift) );
          if (memoryShift)
            iter += memoryShift;
        }
      }

      /* Go through pivot row. */
      for (ListMat64Nonzero* nz = listmatrix->rowElements[pivotRow].head.right;
        nz != &listmatrix->rowElements[pivotRow].head; nz = nz->right)
      {
        size_t column = nz->column;
        if (denseProcessed[column])
          continue;

        /* Apply unimodular transformation to both entries. */
        int64_t p_old = densePivot[column];
        int64_t p_new = U_11 * p_old;
        int64_t o_new = U_21 * p_old;

        densePivot[column] = p_new;
        ptrdiff_t memoryShift;
        CMRdbgMsg(8, "Inserting into other row %ld in column %ld. Value is %ld\n", otherRow, column, o_new);
        CMR_CALL( CMRlistmat64Insert(cmr, listmatrix, otherRow, column, o_new, 0, &memoryShift) );
        if (memoryShift)
          nz += memoryShift;
      }

      /* Go through other row again to remove nonzeros with value 0. */
      for (ListMat64Nonzero* iter = listmatrix->rowElements[otherRow].head.right;
        iter != &listmatrix->rowElements[otherRow].head; )
      {
        ListMat64Nonzero* nz = iter;
        iter = iter->right;
        denseProcessed[nz->column] = false;
        CMRdbgMsg(8, "Deselecting nonzero at %ld,%ld\n", nz->row, nz->column);
        if (nz->value == 0)
        {
          CMRdbgMsg(8, "Removing nonzero at %ld,%ld\n", nz->row, nz->column);
          CMR_CALL( CMRlistmat64Delete(cmr, listmatrix, nz) );
        }
      }
    }

    /* Go through pivot row and update the linked list data. */
    for (ListMat64Nonzero* iter = listmatrix->rowElements[pivotRow].head.right;
      iter != &listmatrix->rowElements[pivotRow].head; )
    {
      ListMat64Nonzero* nz = iter;
      iter = iter->right;
      int64_t entry = densePivot[nz->column];
      densePivot[nz->column] = 0;

      CMRdbgMsg(6, "Updating entry in pivot row %ld, column %ld, value: %ld\n\n", pivotRow, nz->column, entry);
      if (entry != 0)
        nz->value = entry;
      else
        CMR_CALL( CMRlistmat64Delete(cmr, listmatrix, nz) );
    }
  }

  /* If requested, write the resulting (transposed) matrix back into an int matrix. */
  if (!isIntTooSmall && (presult || ptranspose))
  {
    CMR_INTMAT* result = NULL;
    CMR_CALL( CMRintmatCreate(cmr, &result, listmatrix->numRows, listmatrix->numColumns, listmatrix->numNonzeros) );
    result->rowSlice[0] = 0;
    result->numNonzeros = 0;
    for (size_t permRow = 0; permRow < matrix->numRows; ++permRow)
    {
      size_t row = permutations->rows[permRow];
      for (ListMat64Nonzero* nz = listmatrix->rowElements[row].head.right;
        nz != &listmatrix->rowElements[row].head; nz = nz->right)
      {
        size_t permColumn = originalColumnsToPermutedColumns[nz->column];
        result->entryValues[result->numNonzeros] = nz->value;
        result->entryColumns[result->numNonzeros] = permColumn;
        result->numNonzeros++;
      }
      result->rowSlice[permRow+1] = result->numNonzeros;
    }
    CMR_CALL( CMRintmatSortNonzeros(cmr, result) );

    if (ptranspose)
      CMR_CALL( CMRintmatTranspose(cmr, result, ptranspose) );

    if (presult)
      *presult = result;
    else
      CMR_CALL( CMRintmatFree(cmr, &result) );
  }

  CMR_CALL( CMRfreeStackArray(cmr, &originalColumnsToPermutedColumns) );
  CMR_CALL( CMRfreeStackArray(cmr, &originalRowsToPermutedRows) );
  CMR_CALL( CMRfreeStackArray(cmr, &otherRowInfos) );
  CMR_CALL( CMRfreeStackArray(cmr, &denseProcessed) );
  CMR_CALL( CMRfreeStackArray(cmr, &densePivot) );

  CMR_CALL( CMRlistmat64Free(cmr, &listmatrix) );

#if defined(CMR_WITH_GMP)

  if (isIntTooSmall)
  {
    CMR_CALL( CMRintmatComputeUpperDiagonalGMP(cmr, matrix, invert, prank, *ppermutations, presult, ptranspose) );

    return CMR_OKAY;
  }

#endif /* CMR_WITH_GMP */

  return isIntTooSmall ? CMR_ERROR_OVERFLOW : CMR_OKAY;
}
