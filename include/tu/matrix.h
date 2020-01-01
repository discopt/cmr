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
} TU_SPARSE_DOUBLE;

/**
 * \brief Clears the arrays in a sparse int matrix.
 */

TU_EXPORT
void TUclearSparseDouble(
  TU_SPARSE_DOUBLE* sparse /**< Sparse double matrix */
);
  
/**
 * \brief Row-wise representation of sparse int matrix
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
} TU_SPARSE_INT;

/**
 * \brief Clears the arrays in a sparse int matrix.
 */

TU_EXPORT
void TUclearSparseInt(
  TU_SPARSE_INT* sparse /**< Sparse int matrix */
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
} TU_SPARSE_CHAR;

/**
 * \brief Clears the arrays in a sparse char matrix.
 */

TU_EXPORT
void TUclearSparseChar(
  TU_SPARSE_CHAR* sparse /**< Sparse char matrix */
);

/**
 * \brief Prints a sparse double matrix.
 */

TU_EXPORT
void TUprintSparseAsDenseDouble(
  FILE* stream, /**< File stream to print to */
  TU_SPARSE_DOUBLE* sparse, /**< Sparse double matrix */
  char zeroChar, /**< Character to print for a zero */
  bool header /**< Whether to print row and column indices */
);

/**
 * \brief Prints a sparse int matrix.
 */

TU_EXPORT
void TUprintSparseAsDenseInt(
  FILE* stream, /**< File stream to print to */
  TU_SPARSE_INT* sparse, /**< Sparse int matrix */
  char zeroChar, /**< Character to print for a zero */
  bool header /**< Whether to print row and column indices */
);

/**
 * \brief Prints a sparse char matrix.
 */

TU_EXPORT
void TUprintSparseAsDenseChar(
  FILE* stream, /**< File stream to print to */
  TU_SPARSE_CHAR* sparse, /**< Sparse char matrix */
  char zeroChar, /**< Character to print for a zero */
  bool header /**< Whether to print row and column indices */
);

/**
 * \brief Checks whether two sparse double matrices are equal.
 */

TU_EXPORT
bool TUcheckSparseEqualDouble(
  TU_SPARSE_DOUBLE* matrix1, /**< First matrix */
  TU_SPARSE_DOUBLE* matrix2 /**< Second matrix */
);

/**
 * \brief Checks whether two sparse int matrices are equal.
 */

TU_EXPORT
bool TUcheckSparseEqualInt(
  TU_SPARSE_INT* matrix1, /**< First matrix */
  TU_SPARSE_INT* matrix2 /**< Second matrix */
);


/**
 * \brief Checks whether two sparse char matrices are equal.
 */

TU_EXPORT
bool TUcheckSparseEqualChar(
  TU_SPARSE_CHAR* matrix1, /**< First matrix */
  TU_SPARSE_CHAR* matrix2 /**< Second matrix */
);

/**
 * \brief Checks whether two sparse double matrices are transposes of each other.
 */

TU_EXPORT
bool TUcheckSparseTransposeDouble(
  TU_SPARSE_DOUBLE* matrix1, /**< First matrix */
  TU_SPARSE_DOUBLE* matrix2 /**< Second matrix */
);


/**
 * \brief Checks whether two sparse int matrices are transposes of each other.
 */

TU_EXPORT
bool TUcheckSparseTransposeInt(
  TU_SPARSE_INT* matrix1, /**< First matrix */
  TU_SPARSE_INT* matrix2 /**< Second matrix */
);


/**
 * \brief Checks whether two sparse char matrices are transposes of each other.
 */

TU_EXPORT
bool TUcheckSparseTransposeChar(
  TU_SPARSE_CHAR* matrix1, /**< First matrix */
  TU_SPARSE_CHAR* matrix2 /**< Second matrix */
);

/**
 * \brief Checks whether sparse double matrix has each row sorted by minor.
 */

TU_EXPORT
bool TUcheckSparseSortedDouble(
  TU_SPARSE_DOUBLE* sparse /** Sparse double matrix */
);

/**
 * \brief Checks whether sparse int matrix has each row sorted by minor.
 */

TU_EXPORT
bool TUcheckSparseSortedInt(
  TU_SPARSE_INT* sparse /** Sparse int matrix */
);

/**
 * \brief Checks whether sparse char matrix has each row sorted by minor.
 */

TU_EXPORT
bool TUcheckSparseSortedChar(
  TU_SPARSE_CHAR* sparse /** Sparse char matrix */
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
 * \brief Creates a submatrix of a sparse char matrix explicitly.
 */
TU_EXPORT
void TUfilterSubmatrixChar(
  TU_SPARSE_CHAR* matrix,   /**< Given matrix */
  TU_SUBMATRIX* submatrix,  /**< Specified submatrix */
  TU_SPARSE_CHAR* result    /**< Resulting submatrix as a sparse char matrix. */
);

/**
 * \brief Checks if sparse double matrix has only entries in {-1, 0, +1} with absolute error
 * tolerance \p epsilon.
 */

TU_EXPORT
bool TUisTernaryDouble(
  TU_SPARSE_DOUBLE* sparse, /**< Sparse char matrix */
  double epsilon,           /**< Absolute error tolerance */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

/**
 * \brief Checks if sparse int matrix has only entries in {-1, 0, +1}.
 */

TU_EXPORT
bool TUisTernaryInt(
  TU_SPARSE_INT* sparse,    /**< Sparse char matrix */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

/**
 * \brief Checks if sparse char matrix has only entries in {-1, 0, +1}.
 */

TU_EXPORT
bool TUisTernaryChar(
  TU_SPARSE_CHAR* sparse,   /**< Sparse char matrix */
  TU_SUBMATRIX** submatrix  /**< If not \c NULL, a non-ternary entry is stored in \c *submatrix. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_MATRIX_H */
