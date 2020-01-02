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
} TU_MATRIX_DOUBLE;

/**
 * \brief Clears the arrays in an int matrix.
 */

TU_EXPORT
void TUclearMatrixDouble(
  TU_MATRIX_DOUBLE* matrix /**< Double matrix */
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
} TU_MATRIX_INT;

/**
 * \brief Clears the arrays in an int matrix.
 */

TU_EXPORT
void TUclearMatrixInt(
  TU_MATRIX_INT* matrix /**< Int matrix */
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
} TU_MATRIX_CHAR;

/**
 * \brief Clears the arrays in a char matrix.
 */

TU_EXPORT
void TUclearMatrixChar(
  TU_MATRIX_CHAR* matrix /**< Char matrix */
);

/**
 * \brief Prints a double matrix.
 */

TU_EXPORT
void TUprintMatrixDenseDouble(
  FILE* stream,             /**< File stream to print to */
  TU_MATRIX_DOUBLE* matrix, /**< Double matrix */
  char zeroChar,            /**< Character to print for a zero */
  bool header               /**< Whether to print row and column indices */
);

/**
 * \brief Prints an int matrix.
 */

TU_EXPORT
void TUprintMatrixDenseInt(
  FILE* stream,           /**< File stream to print to */
  TU_MATRIX_INT* matrix,  /**< Int matrix */
  char zeroChar,          /**< Character to print for a zero */
  bool header             /**< Whether to print row and column indices */
);

/**
 * \brief Prints a char matrix.
 */

TU_EXPORT
void TUprintMatrixDenseChar(
  FILE* stream,           /**< File stream to print to */
  TU_MATRIX_CHAR* matrix, /**< Char matrix */
  char zeroChar,          /**< Character to print for a zero */
  bool header             /**< Whether to print row and column indices */
);

/**
 * \brief Checks whether two double matrices are equal.
 */

TU_EXPORT
bool TUcheckMatrixEqualDouble(
  TU_MATRIX_DOUBLE* matrix1,  /**< First matrix */
  TU_MATRIX_DOUBLE* matrix2   /**< Second matrix */
);

/**
 * \brief Checks whether two int matrices are equal.
 */

TU_EXPORT
bool TUcheckMatrixEqualInt(
  TU_MATRIX_INT* matrix1, /**< First matrix */
  TU_MATRIX_INT* matrix2  /**< Second matrix */
);


/**
 * \brief Checks whether two char matrices are equal.
 */

TU_EXPORT
bool TUcheckMatrixEqualChar(
  TU_MATRIX_CHAR* matrix1,  /**< First matrix */
  TU_MATRIX_CHAR* matrix2   /**< Second matrix */
);

/**
 * \brief Checks whether two double matrices are transposes of each other.
 */

TU_EXPORT
bool TUcheckMatrixTransposeDouble(
  TU_MATRIX_DOUBLE* matrix1, /**< First matrix */
  TU_MATRIX_DOUBLE* matrix2 /**< Second matrix */
);


/**
 * \brief Checks whether two int matrices are transposes of each other.
 */

TU_EXPORT
bool TUcheckMatrixTransposeInt(
  TU_MATRIX_INT* matrix1, /**< First matrix */
  TU_MATRIX_INT* matrix2 /**< Second matrix */
);


/**
 * \brief Checks whether two char matrices are transposes of each other.
 */

TU_EXPORT
bool TUcheckMatrixTransposeChar(
  TU_MATRIX_CHAR* matrix1, /**< First matrix */
  TU_MATRIX_CHAR* matrix2 /**< Second matrix */
);

/**
 * \brief Checks whether double matrix has each row sorted by minor.
 */

TU_EXPORT
bool TUcheckMatrixSortedDouble(
  TU_MATRIX_DOUBLE* matrix /** Double matrix */
);

/**
 * \brief Checks whether int matrix has each row sorted by minor.
 */

TU_EXPORT
bool TUcheckMatrixSortedInt(
  TU_MATRIX_INT* matrix /** Int matrix */
);

/**
 * \brief Checks whether char matrix has each row sorted by minor.
 */

TU_EXPORT
bool TUcheckMatrixSortedChar(
  TU_MATRIX_CHAR* matrix /** Char matrix */
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
  TU_SUBMATRIX** submatrix, /** Pointer to submatrix */
  int numRows, /**< Number of rows */
  int numColumns /**< Number of columns */
);

/**
 * \brief Creates a 1x1 submatrix.
 */

TU_EXPORT
void TUcreateSubmatrix1x1(
  TU_SUBMATRIX** submatrix, /**< Pointer to submatrix */
  int row, /**< Row of entry */
  int column /**< Column of entry */
);

/**
 * \brief Frees a submatrix.
 */
TU_EXPORT
void TUfreeSubmatrix(
  TU_SUBMATRIX** submatrix /**< Pointer to submatrix */
);


/**
 * \brief Creates a submatrix of a char matrix explicitly.
 */
TU_EXPORT
void TUfilterSubmatrixChar(
  TU_MATRIX_CHAR* matrix,   /**< Given matrix */
  TU_SUBMATRIX* submatrix,  /**< Specified submatrix */
  TU_MATRIX_CHAR* result    /**< Resulting submatrix as a char matrix. */
);

/**
 * \brief Checks if double matrix has only entries in {-1, 0, +1} with absolute error
 * tolerance \p epsilon.
 */

TU_EXPORT
bool TUisTernaryDouble(
  TU_MATRIX_DOUBLE* matrix, /**< Char matrix */
  double epsilon,           /**< Absolute error tolerance */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

/**
 * \brief Checks if int matrix has only entries in {-1, 0, +1}.
 */

TU_EXPORT
bool TUisTernaryInt(
  TU_MATRIX_INT* matrix,    /**< Char matrix */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

/**
 * \brief Checks if char matrix has only entries in {-1, 0, +1}.
 */

TU_EXPORT
bool TUisTernaryChar(
  TU_MATRIX_CHAR* matrix,   /**< Char matrix */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_MATRIX_H */
