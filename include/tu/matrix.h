#ifndef TU_MATRIX_H
#define TU_MATRIX_H

#include <tu/env.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Row-wise representation of sparse double matrix
 */

typedef struct
{
  /**
   * \brief Number of rows
   */
  int numRows;

  /**
   * \brief Number of columns
   */
  int numColumns;

  /**
   * \brief Number of nonzeros
   */
  int numNonzeros;

  /**
   * \brief Array mapping each row to its first entry
   */
  int* rowStarts;

  /**
   * \brief Array mapping each entry to its column
   */
  int* entryColumns;

  /**
   * \brief Array mapping each entry to its value
   */
  double* entryValues;
} TU_DBLMAT;

/**
 * \brief Creates a double matrix of size \p numRows times \p numColumns with \p numNonzeros
 *        nonzeros. The row starts and entries are allocated but not initialized.
 */

TU_EXPORT
void TUdblmatCreate(
  TU* tu,                     /**< TU environment. */
  TU_DBLMAT** matrix,  /**< Pointer for storing the created matrix. */
  int numRows,                /**< Number of rows. */
  int numColumns,             /**< Number of columns. */
  int numNonzeros             /**< Number of nonzeros. */
);

/**
 * \brief Frees the memory of a double matrix.
 */

TU_EXPORT
void TUdblmatFree(
  TU* tu,                   /**< TU environment. */
  TU_DBLMAT** matrix /**< Double matrix */
);

/**
 * \brief Copy a double matrix.
 */
TU_EXPORT
void TUdblmatCopy(
  TU* tu,                   /**< TU environment. */
  TU_DBLMAT* matrix, /**< Given matrix. */
  TU_DBLMAT** result /**< Pointer to created copy of matrix. */
);

/**
 * \brief Create the transposed double matrix.
 */
TU_EXPORT
void TUdblmatTranspose(
  TU* tu,                   /**< TU environment. */
  TU_DBLMAT* matrix, /**< Given matrix. */
  TU_DBLMAT** result /**< Pointer to created transpose of matrix. */
);

/**
 * \brief Row-wise representation of sparse int matrix.
 */

typedef struct
{
  /**
   * \brief Number of rows
   */
  int numRows;

  /**
   * \brief Number of columns
   */
  int numColumns;

  /**
   * \brief Number of nonzeros
   */
  int numNonzeros;

  /**
   * \brief Array mapping each row to its first entry
   */
  int* rowStarts;

  /**
   * \brief Array mapping each entry to its column
   */
  int* entryColumns;

  /**
   * \brief Array mapping each entry to its value
   */
  int* entryValues;
} TU_INTMAT;


/**
 * \brief Creates an int matrix of size \p numRows times \p numColumns with \p numNonzeros
 *        nonzeros. The row starts and entries are allocated but not initialized.
 */

TU_EXPORT
void TUintmatCreate(
  TU* tu,                 /**< TU environment. */
  TU_INTMAT** matrix, /**< Pointer for storing the created matrix. */
  int numRows,            /**< Number of rows. */
  int numColumns,         /**< Number of columns. */
  int numNonzeros         /**< Number of nonzeros. */
);

/**
 * \brief Frees the memory of an int matrix.
 */

TU_EXPORT
void TUintmatFree(
  TU* tu,                 /**< TU environment. */
  TU_INTMAT** matrix  /**< Double matrix */
);

/**
 * \brief Copy an int matrix.
 */
TU_EXPORT
void TUintmatCopy(
  TU* tu,                 /**< TU environment. */
  TU_INTMAT* matrix,  /**< Given matrix. */
  TU_INTMAT** result  /**< Pointer to created copy of matrix. */
);

/**
 * \brief Create the transposed int matrix.
 */
TU_EXPORT
void TUintmatTranspose(
  TU* tu,                   /**< TU environment. */
  TU_INTMAT* matrix, /**< Given matrix. */
  TU_INTMAT** result /**< Pointer to created transpose of matrix. */
);

/**
 * \brief Row-wise representation of sparse char matrix
 */

typedef struct
{
  /**
   * \brief Number of rows
   */
  int numRows;

  /**
   * \brief Number of columns
   */
  int numColumns;

  /**
   * \brief Number of nonzeros
   */
  int numNonzeros;

  /**
   * \brief Array mapping each row to its first entry
   */
  int* rowStarts;

  /**
   * \brief Array mapping each entry to its column
   */
  int* entryColumns;

  /**
   * \brief Array mapping each entry to its value
   */
  char* entryValues;
} TU_CHRMAT;



/**
 * \brief Creates a char matrix of size \p numRows times \p numColumns with \p numNonzeros
 *        nonzeros. The row starts and entries are allocated but not initialized.
 */

TU_EXPORT
void TUchrmatCreate(
  TU* tu,                   /**< TU environment. */
  TU_CHRMAT** matrix,  /**< Pointer for storing the created matrix. */
  int numRows,              /**< Number of rows. */
  int numColumns,           /**< Number of columns. */
  int numNonzeros           /**< Number of nonzeros. */
);

/**
 * \brief Frees the memory of a char matrix.
 */

TU_EXPORT
void TUchrmatFree(
  TU* tu,                 /**< TU environment. */
  TU_CHRMAT** matrix /**< Char matrix */
);

/**
 * \brief Copy a char matrix.
 */
TU_EXPORT
void TUchrmatCopy(
  TU* tu,                 /**< TU environment. */
  TU_CHRMAT* matrix, /**< Given matrix. */
  TU_CHRMAT** result /**< Pointer to created copy of matrix. */
);

/**
 * \brief Create the transposed char matrix.
 */
TU_EXPORT
void TUchrmatTranspose(
  TU* tu,                   /**< TU environment. */
  TU_CHRMAT* matrix, /**< Given matrix. */
  TU_CHRMAT** result /**< Pointer to created transpose of matrix. */
);

/**
 * \brief Prints a double matrix.
 */

TU_EXPORT
void TUdblmatPrintDense(
  FILE* stream,             /**< File stream to print to */
  TU_DBLMAT* matrix, /**< Double matrix */
  char zeroChar,            /**< Character to print for a zero */
  bool header               /**< Whether to print row and column indices */
);

/**
 * \brief Prints an int matrix.
 */

TU_EXPORT
void TUintmatPrintDense(
  FILE* stream,           /**< File stream to print to */
  TU_INTMAT* matrix,  /**< Int matrix */
  char zeroChar,          /**< Character to print for a zero */
  bool header             /**< Whether to print row and column indices */
);

/**
 * \brief Prints a char matrix.
 */

TU_EXPORT
void TUchrmatPrintDense(
  FILE* stream,           /**< File stream to print to */
  TU_CHRMAT* matrix, /**< Char matrix */
  char zeroChar,          /**< Character to print for a zero */
  bool header             /**< Whether to print row and column indices */
);

/**
 * \brief Checks whether two double matrices are equal.
 */

TU_EXPORT
bool TUdblmatCheckEqual(
  TU_DBLMAT* matrix1,  /**< First matrix */
  TU_DBLMAT* matrix2   /**< Second matrix */
);

/**
 * \brief Checks whether two int matrices are equal.
 */

TU_EXPORT
bool TUintmatCheckEqual(
  TU_INTMAT* matrix1, /**< First matrix */
  TU_INTMAT* matrix2  /**< Second matrix */
);


/**
 * \brief Checks whether two char matrices are equal.
 */

TU_EXPORT
bool TUchrmatCheckEqual(
  TU_CHRMAT* matrix1,  /**< First matrix */
  TU_CHRMAT* matrix2   /**< Second matrix */
);

/**
 * \brief Checks whether two double matrices are transposes of each other.
 */

TU_EXPORT
bool TUdblmatCheckTranspose(
  TU_DBLMAT* matrix1,  /**< First matrix */
  TU_DBLMAT* matrix2   /**< Second matrix */
);


/**
 * \brief Checks whether two int matrices are transposes of each other.
 */

TU_EXPORT
bool TUintmatCheckTranspose(
  TU_INTMAT* matrix1, /**< First matrix */
  TU_INTMAT* matrix2  /**< Second matrix */
);

/**
 * \brief Checks whether two char matrices are transposes of each other.
 */

TU_EXPORT
bool TUchrmatCheckTranspose(
  TU_CHRMAT* matrix1, /**< First matrix */
  TU_CHRMAT* matrix2 /**< Second matrix */
);

/**
 * \brief Checks whether double matrix has each row sorted by minor.
 */

TU_EXPORT
bool TUdblmatCheckSorted(
  TU_DBLMAT* matrix /** Double matrix */
);

/**
 * \brief Checks whether int matrix has each row sorted by minor.
 */

TU_EXPORT
bool TUintmatCheckSorted(
  TU_INTMAT* matrix /** Int matrix */
);

/**
 * \brief Checks whether char matrix has each row sorted by minor.
 */

TU_EXPORT
bool TUchrmatCheckSorted(
  TU_CHRMAT* matrix /** Char matrix */
);


/**
 * \brief Row and column indices for a submatrix
 *
 * Does not contain information about the matrix it refers to.
 */
typedef struct
{
  /**
   * \brief Number of rows
   */
  int numRows;
  /**
   * \brief Array with row indices
   */
  int* rows;
  /**
   * \brief Number of columns
   */
  int numColumns;
  /**
   * \brief Array with column indices
   */
  int* columns;
} TU_SUBMATRIX;

/**
 * \brief Creates a submatrix of given size.
 *
 * Only allocates the memory. Use \ref TUgetSubmatrixRows and \ref TUgetSubmatrixColumns to modify
 * the row and column indices, respectively.
 */
TU_EXPORT
void TUcreateSubmatrix(
  TU* tu,                   /**< TU environment. */
  TU_SUBMATRIX** submatrix, /** Pointer to where the submatrix is to be stored. */
  int numRows,              /**< Number of rows */
  int numColumns            /**< Number of columns */
);

/**
 * \brief Creates a 1x1 submatrix.
 */

TU_EXPORT
void TUcreateSubmatrix1x1(
  TU* tu,                   /**< TU environment. */
  TU_SUBMATRIX** submatrix, /**< Pointer to submatrix */
  int row,                  /**< Row of entry */
  int column                /**< Column of entry */
);

/**
 * \brief Frees a submatrix.
 */
TU_EXPORT
void TUfreeSubmatrix(
  TU* tu,                   /**< TU environment */
  TU_SUBMATRIX** submatrix  /**< Pointer to submatrix */
);

/**
 * \brief Creates a submatrix of a char matrix explicitly.
 */
TU_EXPORT
void TUfilterCharSubmatrix(
  TU* tu,                   /**< TU environment. */
  TU_CHRMAT* matrix,   /**< Given matrix */
  TU_SUBMATRIX* submatrix,  /**< Specified submatrix */
  TU_CHRMAT** result   /**< Resulting submatrix as a char matrix. */
);


/**
 * \brief Checks if double matrix has only entries in {0, 1} with absolute error
 * tolerance \p epsilon.
 */

TU_EXPORT
bool TUisBinaryDouble(
  TU* tu,                   /**< TU environment. */
  TU_DBLMAT* matrix, /**< Char matrix */
  double epsilon,           /**< Absolute error tolerance */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

/**
 * \brief Checks if int matrix has only entries in {0, 1}.
 */

TU_EXPORT
bool TUisBinaryInt(
  TU* tu,                   /**< TU environment. */
  TU_INTMAT* matrix,    /**< Char matrix */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

/**
 * \brief Checks if matrix has only entries in {0, 1}.
 */

TU_EXPORT
bool TUisBinaryChar(
  TU* tu,                   /**< TU environment. */
  TU_CHRMAT* matrix,   /**< Char matrix */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

/**
 * \brief Checks if double matrix has only entries in {-1, 0, +1} with absolute error
 * tolerance \p epsilon.
 */

TU_EXPORT
bool TUisTernaryDouble(
  TU* tu,                   /**< TU environment. */
  TU_DBLMAT* matrix, /**< Char matrix */
  double epsilon,           /**< Absolute error tolerance */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

/**
 * \brief Checks if int matrix has only entries in {-1, 0, +1}.
 */

TU_EXPORT
bool TUisTernaryInt(
  TU* tu,                   /**< TU environment. */
  TU_INTMAT* matrix,    /**< Char matrix */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

/**
 * \brief Checks if char matrix has only entries in {-1, 0, +1}.
 */

TU_EXPORT
bool TUisTernaryChar(
  TU* tu,                   /**< TU environment. */
  TU_CHRMAT* matrix,   /**< Char matrix */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_MATRIX_H */
