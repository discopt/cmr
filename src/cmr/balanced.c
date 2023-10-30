//#define CMR_DEBUG /* Uncomment to debug */

#include <cmr/balanced.h>

#include <cmr/camion.h>

#include <cmr/series_parallel.h>
#include "env_internal.h"

#include "matrix_internal.h"
#include "camion_internal.h"

#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>

/* Start selecting columns until maxDepth many are selected, then do final balancedness test on resulting submatrix. */ 
CMR_ERROR enumerateCols(CMR* cmr, size_t depth, size_t maxDepth, CMR_CHRMAT* matrix, size_t* selectedRows, size_t* NNZ2Columns, size_t NNZ2ColCount, size_t* selectedCols, size_t* NNZperRow, bool* pisBalanced, CMR_SUBMAT** psubmatrix)
{
  // NOTE THAT selectedCols CONTAINS THE ENTRIES OF NNZ2Columns THAT ARE SELECTED!!
  /* We can stop iterating if we found that the matrix is unbalanced. */
  if (!*pisBalanced) return CMR_OKAY;

  if (depth == maxDepth)
  {
    /* The final step, we have selected maxDepth columns. */
    CMRdbgMsg(2, "Selected columns: ");
    for (size_t i = 0; i < maxDepth; ++i)
    {
      CMRdbgMsg(0, "%d ", NNZ2Columns[selectedCols[i]]);
    }
    CMRdbgMsg(0, "\n");

    /* Count the sum of the entries in the selected rows and columns. */
    size_t sum = 0;
    for (size_t i = 0; i < maxDepth; ++i)
    {
      size_t row = selectedRows[i];
      /* Counts how many NZ entries have been found in the row, so that we can stop once we found 2. */
      int NNZfound = 0;
      for (size_t j = 0; j < maxDepth; ++j) 
      {
        size_t col = NNZ2Columns[selectedCols[j]];
        size_t entry = SIZE_MAX;
        CMR_CALL( CMRchrmatFindEntry(matrix, row, col, &entry) );
        if (entry != SIZE_MAX) {
          NNZfound++;
          sum += matrix->entryValues[entry];
        }
        CMRdbgMsg(3, "value at %ld %ld: %ld\n", row, col, matrix->entryValues[entry]);
      }
    }
    CMRdbgMsg(4, "Sum of entries: %ld\n", sum);
    if (sum % 4 != 0) 
    {
      /* This submatrix is unbalanced. */
      CMRdbgMsg(4, "Unbalanced submatrix found!\n");
      *pisBalanced = false;

      /* Output the unbalanced submatrix. */
      CMRdbgMsg(4, "Creating submatrix...\n");
      CMRdbgMsg(5, "TEST0\n");
      CMR_CALL( CMRsubmatCreate(cmr, maxDepth, maxDepth, psubmatrix) );
      CMRdbgMsg(5, "TEST1\n");
      CMR_SUBMAT* submatrix = *psubmatrix;
      CMRdbgMsg(5, "TEST2\n");
      for (size_t i = 0; i < maxDepth; ++i)
      {
        CMRdbgMsg(5, "i: %ld\n", i);
        submatrix->rows[i] = selectedRows[i];
        submatrix->columns[i] = NNZ2Columns[selectedCols[i]];
      }
      CMRdbgMsg(4, "Submatrix created!\n");
    }
  } 
  else 
  {
    /* The recursion step, we pick one column and enumerate over the other columns. */
    for (size_t col = (depth == 0) ? 0 : selectedCols[depth - 1] + 1; col < NNZ2ColCount - maxDepth + depth + 1; ++col)
    {
      // size_t col = NNZ2Columns[i];
      /* First check if adding the column would satisfy the following rules: 
        1. No entry in NNZperRow should become >2.
        2. Adding the column should make at least one entry in NNZperRow change from 0 to 1. */
      /* Set to true if adding the column would add a 3rd NZ to one of the rows. */
      bool rule1isBroken = false;
      /* Set to false once a NZ is added to a row that previously contained 0 NZ. */
      bool rule2isBroken = true;
      /* Contains a bool for each row. Each element is true if adding the column adds a NZ to the row. */
      bool* NZadditions = NULL;
      CMR_CALL( CMRallocStackArray(cmr, &NZadditions, maxDepth) );
      /* Counts how many NZ entries have been found in the column, so that we can stop once we found 2. */
      int NNZfound = 0;
      for (size_t row = 0; NNZfound < 2; ++row) {
        size_t entry;
        // TODO Is there a better way to get the indices of the NZ entries in a column of a matrix in the sparse format?
        CMR_CALL( CMRchrmatFindEntry(matrix, selectedRows[row], NNZ2Columns[col], &entry) );
        // CMRdbgMsg(3, "Value at %ld %ld: %ld\n", selectedRows[row], NNZ2Columns[col], matrix->entryValues[entry]);
        if (entry != SIZE_MAX) {
          NNZfound++;
          NZadditions[row] = true;
          if (NNZperRow[row] == 2) {
            /* There are already two NZ entries in this row, and adding this column would make that 3. This breaks rule 1. */
            rule1isBroken = true;
            break;
          } else if (NNZperRow[row] == 0) {
            /* One of the rows goes from NNZ = 0 to NNZ = 1. */
            rule2isBroken = false;
          }
        } else {
          NZadditions[row] = false;
        }
      }
      /* If the addition of this column breaks rule 1, it can be skipped. */
      if (rule1isBroken) {
        // CMRdbgMsg(2, "Rule 1 broken.\n");
        CMR_CALL( CMRfreeStackArray(cmr, &NZadditions) ); 
        continue;
      }
      /* If the addition of this column breaks rule 2, it can be skipped, except when we are choosing the final column, then rule 2 is always broken. */
      if (rule2isBroken && (depth != maxDepth - 1)) {
        // CMRdbgMsg(2, "Rule 2 broken. Current depth: %ld. maxDepth: %ld. Debug: %s\n", depth, maxDepth, depth != maxDepth - 1 ? "true" : "false");
        CMR_CALL( CMRfreeStackArray(cmr, &NZadditions) ); 
        continue;
      }
      
      /* Add to NNZperRow counter. */
      for (size_t row = 0; row < maxDepth; ++row) {
        if (NZadditions[row]) {NNZperRow[row]++;}
      }
      selectedCols[depth] = col;

      enumerateCols(cmr, depth + 1, maxDepth, matrix, selectedRows, NNZ2Columns, NNZ2ColCount, selectedCols, NNZperRow, pisBalanced, psubmatrix);
      
      /* Remove from NNZperRow counter. */
      for (size_t row = 0; row < maxDepth; ++row) {
        if (NZadditions[row]) {NNZperRow[row]--;}
      }
      
      CMR_CALL( CMRfreeStackArray(cmr, &NZadditions) );
    }
  }

  return CMR_OKAY;
}

/* Start selecting rows until maxDepth many are selected, then start selecting columns. */
CMR_ERROR enumerateRows(CMR* cmr, size_t depth, size_t maxDepth, CMR_CHRMAT* matrix, size_t* selectedRows, bool* pisBalanced, CMR_SUBMAT** psubmatrix, double timeLimit)
{
  /* We can stop iterating if we found that the matrix is unbalanced. */
  if (!*pisBalanced) return CMR_OKAY;

  if (depth == maxDepth)
  {
    /* The final step, we have selected maxDepth rows. */
    CMRdbgMsg(0, "Selected rows: ");
    for (size_t i = 0; i < maxDepth; ++i)
    {
      CMRdbgMsg(0, "%d ", selectedRows[i]);
    }
    CMRdbgMsg(0, "\n");

    /* Here we keep track of NNZ in each column of the selected rows. */
    size_t* NNZPerCol = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &NNZPerCol, matrix->numColumns) );
    for (size_t col = 0; col < matrix->numColumns; ++col) { NNZPerCol[col] = 0; }

    for (size_t i = 0; i < maxDepth; ++i)
    {
      size_t row = selectedRows[i];
      for (size_t entry = matrix->rowSlice[row]; entry < matrix->rowSlice[row + 1]; ++entry) {
        NNZPerCol[matrix->entryColumns[entry]]++;
      }
    }
    /* Counts the number of columns with NNZ = 2 in the selected rows. */
    size_t NNZ2ColCount = 0;
    for (size_t col = 0; col < matrix->numColumns; ++col) {
      if (NNZPerCol[col] == 2) {
        NNZ2ColCount++;
      }
    }    
    CMRdbgMsg(1, "Number of cols with NNZ=2: %ld\n", NNZ2ColCount);
    if (NNZ2ColCount < maxDepth) {
      /* There are not enough columns with NNZ=2 to make a square submatrix */
      CMRdbgMsg(1, "Not enough cols to make a square submatrix\n");
      return CMR_OKAY;
    }
    /* Here we store the indices of the columns with NNZ = 2 in the selected rows. */
    size_t* NNZ2Columns = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &NNZ2Columns, NNZ2ColCount) );
    size_t NZ2ColsAllocated = 0;
    for (size_t col = 0; NZ2ColsAllocated < NNZ2ColCount; ++col) {
      CMRdbgMsg(1, "NNZ in col %ld: %ld\n", col, NNZPerCol[col]);
      if (NNZPerCol[col] == 2) {
        NNZ2Columns[NZ2ColsAllocated] = col;
        NZ2ColsAllocated++;
      }
    }

    CMRdbgMsg(1, "Columns with NNZ=2: ");
    for (size_t i = 0; i < NNZ2ColCount; i++) {
      size_t col = NNZ2Columns[i];
      CMRdbgMsg(0, "%ld ", col);
    }
    CMRdbgMsg(0, "\n");

    /* The columns with NNZ=2 are found. */

    // /* Make new reduced submatrix, with NNZ=2 cols and selected rows, then enumerateColumns over that reduced submatrix. */
    // /* First make the new submatrix. */
    // CMR_SUBMAT* submatrix = NULL;
    // CMR_CALL( CMRsubmatCreate(cmr, maxDepth, NZ2Count, &submatrix) );
    // submatrix->columns = NZ2Columns;
    // submatrix->rows = selectedRows;
    // CMR_CHRMAT* colNNZ2Matrix = NULL;
    // CMR_CALL( CMRchrmatZoomSubmat(cmr, matrix, submatrix, &colNNZ2Matrix) );
    // CMRdbgMsg(1, "Cols have NNZ=2 only matrix:\n");
    // CMRchrmatPrintDense(cmr, colNNZ2Matrix, stdout, '0', true);
    // fflush(stdout); 
    //
    // /* Now reduce the submatrix colNNZ2Matrix further, to remove duplicate columns. */
    // bool isSeriesParallel = true;
    // CMR_SP_REDUCTION* reductions = NULL;
    // CMR_CALL( CMRallocStackArray(cmr, &reductions, colNNZ2Matrix->numRows + colNNZ2Matrix->numColumns) );
    // size_t numReductions;
    // CMR_SUBMAT* submatrix2 = NULL;
    // CMR_SUBMAT* violatorSubmatrix = NULL;
    // CMR_SP_STATISTICS* stats = NULL;
    //
    // CMR_CALL( CMRtestTernarySeriesParallel(cmr, colNNZ2Matrix, &isSeriesParallel, reductions, &numReductions, 
    //   &submatrix2, &violatorSubmatrix, stats, timeLimit) );
    //
    // CMR_CALL( CMRfreeStackArray(cmr, &reductions) );
    //
    // CMR_CHRMAT* reducedMatrix = NULL;
    // CMR_CALL( CMRchrmatZoomSubmat(cmr, colNNZ2Matrix, submatrix2, &reducedMatrix) );
    // // CMRdbgMsg(0, "New matrix size: %zdx%zd\n", reducedMatrix->numRows, reducedMatrix->numColumns);
    //
    // CMRdbgMsg(1, "Reduced matrix:\n");
    // CMRchrmatPrintDense(cmr, reducedMatrix, stdout, '0', true);
    // fflush(stdout); 
    //
    // /* Check if reducedMatrix has enough columns to make a square submatrix. */
    // if (reducedMatrix->numColumns < maxDepth) {
    //   CMRdbgMsg(1, "Not enough columns left over to make square submatrix.\n");
    //   return CMR_OKAY;
    // }

    // QUESTION: Is this point only reachable by ending up with a square submatrix? 
    // Let's say that reducedMatrix is n by n+1, then there is a column c that isn't a linear combination of the other columns.
    // c has NNZ=2, so must there be another column c2 such that the NZ values of c and c2 are a 2x2 square? 
    // Then if this 2x2 square is unbalanced it must have been found when iterating over square submatrices of size 2.
    // So if we ever have numCols > numRows, there must be a 2x2 unbalanced hole in the reducedMatrix.
    // So can we assume that numCols == numRows?
    // Then recudedMatrix must be either a hole, or multiple disconnected holes?
    // If it is multiple holes, we would have already iterated over each hole before, as it is smaller than the current size.
    // If we end up at this point, each of these smaller holes must have been balanced.
    // Then can we just count all the nonzero entries, in either case?
    // In the case that it is multiple smaller holes then it is balanced as the sum must still be divisible by 4.
    // In the case that it is a single hole we can check if it is balanced by summing and checking if it is divisible by 4.
    // ANSWER: I misunderstood how exactly the ternaryparallelseriesreduction worked. There are cases where my assumption doesn't hold.

    /* Start enumerating over the columns. */
    /* Holds the selected rows used while iterating over all submatrices of size submatrixSize. */
    size_t* selectedCols = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &selectedCols, maxDepth) );
    /* Array to keep count of the NNZ per row with the current selection of columns. */
    size_t* NNZperRow = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &NNZperRow, maxDepth) );
    for (int i = 0; i < maxDepth; ++i) {
      NNZperRow[i] = 0;
    }

    enumerateCols(cmr, 0, maxDepth, matrix, selectedRows, NNZ2Columns, NNZ2ColCount, selectedCols, NNZperRow, pisBalanced, psubmatrix);
    CMR_CALL( CMRfreeStackArray(cmr, &NNZperRow) );
    CMR_CALL( CMRfreeStackArray(cmr, &selectedCols) );

    CMR_CALL( CMRfreeStackArray(cmr, &NNZ2Columns) );
    CMR_CALL( CMRfreeStackArray(cmr, &NNZPerCol) );
  }
  else
  {
    /* The recursion step, we pick one row and enumerate over the other rows. */
    for (size_t i = (depth == 0) ? 0 : selectedRows[depth - 1] + 1; i < matrix->numRows - maxDepth + depth + 1; ++i)
    {
      selectedRows[depth] = i;
      enumerateRows(cmr, depth + 1, maxDepth, matrix, selectedRows, pisBalanced, psubmatrix, timeLimit);
    }
  }
  return CMR_OKAY;
}

CMR_ERROR CMRtestBalancedEnumeration(CMR* cmr, CMR_CHRMAT* matrix, bool* pisBalanced, CMR_SUBMAT** psubmatrix,
  double timeLimit)
{
  assert(cmr);
  assert(matrix);

  /* clock_t totalClock = clock(); */
  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  printf("DEBUG: CMRtestBalancedEnumeration called!\n");
  // printf("Hello world!\n");

  // printf("Pointer: %p\n", matrix); /* Prints pointer adress of matrix */
  // printf("Number of rows: %zd\n", matrix->numRows);
  // printf("Number of columns: %zd\n", matrix->numColumns);
  // size_t entry = SIZE_MAX;
  // CMR_CALL( CMRchrmatFindEntry(matrix, 1, 1, &entry) ); /* Returns the index of the entry at row = 1, col = 1 */
  // printf("%d\n", matrix->entryValues[entry]); /* Try using CMRchrmatFindEntry. When entry is SIZE_MAX it returns 0? Figure out why */

  /* Test for displaying the matrix */
  // for (size_t row = 0; row < matrix->numRows; ++row)
  // {
  //   for (size_t col = 0; col < matrix->numColumns; ++col)
  //   {
  //     size_t entry = SIZE_MAX;
  //     CMR_CALL( CMRchrmatFindEntry(matrix, row, col, &entry) ); /* Returns the index of the entry */
  //     printf("%d ", matrix->entryValues[entry]);
  //   }
  //   printf("\n");
  // }

  /* We can also display the matrix using: */
  CMRdbgMsg(0, "Input matrix:\n");
  CMRchrmatPrintDense(cmr, matrix, stdout, '0', true);
  // CMRchrmatPrintSparse(cmr, matrix, stdout);
  fflush(stdout); 

  /* We first filter out all the rows with NNZ >= 2 */
  // for (size_t row = 0; row < matrix->numRows; ++row)
  // {
  //   if (matrix->rowSlice[row+1] - matrix->rowSlice[row] >= 2)
  //     printf("filteredRow: %zd\n", row);
  // }
  /* Use https://discopt.github.io/cmr/series__parallel_8h.html#a0cd0d3c42bd90cd3510c67d233524b64 instead */
  bool isSeriesParallel = true;
  CMR_SP_REDUCTION* reductions = NULL;
  CMR_CALL( CMRallocStackArray(cmr, &reductions, matrix->numRows + matrix->numColumns) );
  size_t numReductions;

  // CMR_SUBMAT* submatrix = NULL;
  // CMR_CALL( CMRsubmatCreate(cmr, matrix->numRows, matrix->numColumns, &submatrix));
  // printf("Print the submatrix:\n");
  // CMRsubmatWriteToStream(cmr, submatrix, matrix->numRows, matrix->numColumns, stdout);

  CMR_SUBMAT* submatrix = NULL;
  CMR_SUBMAT* violatorSubmatrix = NULL;
  // CMR_SEPA* separation = NULL;
  CMR_SP_STATISTICS* stats = NULL;

  CMR_CALL( CMRtestTernarySeriesParallel(cmr, matrix, &isSeriesParallel, reductions, &numReductions, 
    &submatrix, &violatorSubmatrix, stats, timeLimit) );

  /* How to print the reductions of the matrix. 
     These are all the steps that are taken when constructing the reduced submatrix. */
  // printf("Now print the reductions\n");
  // for (size_t reduction = 0; reduction < numReductions; ++reduction)
  // {
  //   // printf("%d ", reductions[reduction]);
  //   printf("%s\n", CMRspReductionString(reductions[reduction], NULL));
  // }

  CMR_CALL( CMRfreeStackArray(cmr, &reductions) );
  
  /* How to print the submatrix to the terminal: */
  // printf("Print the submatrix:\n");
  // CMRsubmatWriteToStream(cmr, submatrix, matrix->numRows, matrix->numColumns, stdout);
  // fflush(stdout); 

  // printf("Reduced submatrix size: %zdx%zd\n", submatrix->numRows, submatrix->numColumns);

  /* Matrix object that is a reduced submatrix of the original matrix. */
  CMR_CHRMAT* reducedMatrix = NULL;
  CMR_CALL( CMRchrmatZoomSubmat(cmr, matrix, submatrix, &reducedMatrix) );
  // CMRdbgMsg(0, "New matrix size: %zdx%zd\n", reducedMatrix->numRows, reducedMatrix->numColumns);

  CMRdbgMsg(0, "Reduced matrix:\n");
  CMRchrmatPrintDense(cmr, reducedMatrix, stdout, '0', true);
  fflush(stdout); 


  /* Here comes the rest of the exponential algorithm. */

  /* Start with the assumption that the matrix is balanced, until we find that it is not. Then stop searching. */
  bool isBalancedTemp;
  isBalancedTemp = true;
  /* The maximum size for a square submatrix of the reduced matrix. */
  size_t maxSubmatrixSize;
  maxSubmatrixSize = (reducedMatrix->numRows < reducedMatrix->numColumns) ? reducedMatrix->numRows : reducedMatrix->numColumns;
  

  CMRdbgMsg(0, "Max submatrix size: %d\n", maxSubmatrixSize);

  /* The size of the submatrices that we are now iterating over. */
  for (size_t submatrixSize = 2; submatrixSize <= maxSubmatrixSize; ++submatrixSize)
  {
    CMRdbgMsg(0, "\nCurrently evaluating submatrices of size: %d\n", submatrixSize);
    /* Holds the selected rows used while iterating over all submatrices of size submatrixSize. */
    size_t* selectedRows = NULL;
    CMR_CALL( CMRallocStackArray(cmr, &selectedRows, submatrixSize) );
    enumerateRows(cmr, 0, submatrixSize, reducedMatrix, selectedRows, &isBalancedTemp, psubmatrix, timeLimit);
    CMR_CALL( CMRfreeStackArray(cmr, &selectedRows) );
    
    if (!isBalancedTemp) 
    {
      break;
    }
  }
  
  *pisBalanced = isBalancedTemp;

  return CMR_OKAY;
}

CMR_ERROR CMRtestBalancedGraph(CMR* cmr, CMR_CHRMAT* matrix, bool* pisBalanced, CMR_SUBMAT** psubmatrix,
  double timeLimit)
{
  assert(cmr);
  assert(matrix);

  /* clock_t totalClock = clock(); */
  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  printf("DEBUG: CMRtestBalancedGraph called!\n");

  return CMR_OKAY;
}
CMR_ERROR CMRtestBalanced(CMR* cmr, CMR_CHRMAT* matrix, bool* pisBalanced, CMR_SUBMAT** psubmatrix,
  double timeLimit)
{
  assert(cmr);
  assert(matrix);

  /* clock_t totalClock = clock(); */
  if (!CMRchrmatIsTernary(cmr, matrix, psubmatrix))
    return CMR_OKAY;

  printf("DEBUG: CMRtestBalanced called!\n");
  CMR_CALL( CMRtestBalancedEnumeration(cmr, matrix, pisBalanced, psubmatrix, timeLimit) );

  return CMR_OKAY;
}
