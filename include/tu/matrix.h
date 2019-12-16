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


#ifdef __cplusplus
}
#endif

#endif /* TU_MATRIX_H */
