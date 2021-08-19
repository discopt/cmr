#ifndef CMR_MATRIX_H
#define CMR_MATRIX_H

#include <cmr/env.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
 * \brief Row and column indices for a submatrix
 *
 * Does not contain information about the matrix it refers to.
 */
typedef struct
{
  size_t numRows;     /**< \brief Number of rows. */
  size_t* rows;       /**< \brief Array with row indices. */
  size_t numColumns;  /** \brief Number of columns. */
  size_t* columns;    /** \brief Array with column indices. */
} CMR_SUBMAT;


/**
 * \brief Creates a submatrix of given size.
 *
 * Only allocates the memory. Use rows and columns attributes of *\p psubmatrix to actually set the row and column
 * indices, respectively.
 */
CMR_EXPORT
CMR_ERROR CMRsubmatCreate(
  CMR* cmr,               /**< \ref CMR environment. */
  size_t numRows,         /**< Number of rows */
  size_t numColumns,      /**< Number of columns */
  CMR_SUBMAT** psubmatrix /**< Pointer to where the submatrix is to be stored. */
);

/**
 * \brief Creates a 1 \f$ \times \f$ 1 submatrix.
 */

CMR_EXPORT
CMR_ERROR CMRsubmatCreate1x1(
  CMR* cmr,               /**< \ref CMR environment. */
  size_t row,             /**< Row of entry */
  size_t column,          /**< Column of entry */
  CMR_SUBMAT** psubmatrix /**< Pointer to submatrix */
);

/**
 * \brief Frees a submatrix.
 */

CMR_EXPORT
CMR_ERROR CMRsubmatFree(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SUBMAT** psubmatrix /**< Pointer to submatrix. */
);


/**
 * \brief Row-wise representation of sparse double matrix.
 * 
 * The nonzeros of row \c r of the matrix are stored at positions \ref rowSlice[r] up to, but not including,
 * \ref rowSlice[r+1], of the arrays \ref entryColumns and \ref entryValues.
 * The former contains the respective columns while the latter contains the actual matrix entries of the nonzeros.
 * The nonzeros of each row must be ordered by increasing column.
 * Moreover, no duplicates are allowed and all stored nonzero entries must indeed have a value different from 0.
 */

typedef struct
{
  size_t numRows;       /**< \brief Number of rows. */
  size_t numColumns;    /**< \brief Number of columns. */
  size_t numNonzeros;   /**< \brief Number of and memory allocated for nonzeros. */
  size_t* rowSlice;     /**< \brief Array mapping each row to the index of its first entry. */
  size_t* entryColumns; /**< \brief Array mapping each entry to its column.*/
  double* entryValues;  /**< \brief Array mapping each entry to its value. */
} CMR_DBLMAT;

/**
 * \brief Row-wise representation of sparse int matrix.
 * 
 * The nonzeros of the matrix are stored in the arrays \ref entryColumns and \ref \entryValues, each of length
 * \ref numNonzeros.
 * The nonzeros of row \c r are stored at positions \f$ p \in \{a, a+1, \dotsc, b-2, b-1 \}\f$, where
 * \f$ a := \f$ \ref rowSlice[r] and \f$ b := \f$ \ref rowSlice[r+1].
 * The column is \ref entryColumns[p] while the actual matrix entry is \ref entryValues[p].
 * The nonzeros of each row must be sorted by column in ascending order.
 * Moreover, no duplicates are allowed and all stored entries must indeed be nonzero.
 */

typedef struct
{
  size_t numRows;       /**< \brief Number of rows. */
  size_t numColumns;    /**< \brief Number of columns. */
  size_t numNonzeros;   /**< \brief Number of and memory allocated for nonzeros. */
  size_t * rowSlice;    /**< \brief Array mapping each row to the index of its first entry. */
  size_t* entryColumns; /**< \brief Array mapping each entry to its column.*/
  int* entryValues;     /**< \brief Array mapping each entry to its value. */
} CMR_INTMAT;

/**
 * \brief Row-wise representation of sparse char matrix.
 *
 * The nonzeros of the matrix are stored in the arrays \ref entryColumns and \ref \entryValues, each of length
 * \ref numNonzeros.
 * The nonzeros of row \c r are stored at positions \f$ p \in \{a, a+1, \dotsc, b-2, b-1 \}\f$, where
 * \f$ a := \f$ \ref rowSlice[r] and \f$ b := \f$ \ref rowSlice[r+1].
 * The column is \ref entryColumns[p] while the actual matrix entry is \ref entryValues[p].
 * The nonzeros of each row must be sorted by column in ascending order.
 * Moreover, no duplicates are allowed and all stored entries must indeed be nonzero.
 */

typedef struct
{
  size_t numRows;       /**< \brief Number of rows. */
  size_t numColumns;    /**< \brief Number of columns. */
  size_t numNonzeros;   /**< \brief Number of and memory allocated for nonzeros. */
  size_t* rowSlice;     /**< \brief Array mapping each row to the index of its first entry. */
  size_t* entryColumns; /**< \brief Array mapping each entry to its column.*/
  char* entryValues;    /**< \brief Array mapping each entry to its value. */
} CMR_CHRMAT;

/**
 * \brief Creates a double matrix of size \p numRows \f$ \times \f$ \p numColumns with \p numNonzeros
 *        nonzeros. The actual arrays are allocated but not initialized.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatCreate(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT** presult, /**< Pointer for storing the created matrix. */
  int numRows,          /**< Number of rows. */
  int numColumns,       /**< Number of columns. */
  int numNonzeros       /**< Number of nonzeros. */
);

/**
 * \brief Creates an int matrix of size \p numRows \f$ \times \f$ \p numColumns with \p numNonzeros
 *        nonzeros. The actual arrays are allocated but not initialized.
 */

CMR_EXPORT
CMR_ERROR CMRintmatCreate(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT** presult, /**< Pointer for storing the created matrix. */
  int numRows,          /**< Number of rows. */
  int numColumns,       /**< Number of columns. */
  int numNonzeros       /**< Number of nonzeros. */
);

/**
 * \brief Creates a char matrix of size \p numRows \f$ \times \f$ \p numColumns with \p numNonzeros
 *        nonzeros. The actual arrays are allocated but not initialized.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatCreate(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT** presult, /**< Pointer for storing the created matrix. */
  int numRows,          /**< Number of rows. */
  int numColumns,       /**< Number of columns. */
  int numNonzeros       /**< Number of nonzeros. */
);

/**
 * \brief Frees the memory of a double matrix.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatFree(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT** pmatrix  /**< Pointer to matrix. */
);

/**
 * \brief Frees the memory of an int matrix.
 */

CMR_EXPORT
CMR_ERROR CMRintmatFree(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT** pmatrix  /**< Pointer to matrix. */
);

/**
 * \brief Frees the memory of a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatFree(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT** pmatrix  /**< Pointer to matrix. */
);

/**
 * \brief Changes the number of nonzeros and reallocates corresponding arrays.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatChangeNumNonzeros(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,   /**< Given matrix. */
  size_t newNumNonzeros /**< New number of nonzeros. */ 
);

/**
 * \brief Changes the number of nonzeros and reallocates corresponding arrays.
 */

CMR_EXPORT
CMR_ERROR CMRintmatChangeNumNonzeros(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT* matrix,   /**< A matrix. */
  size_t newNumNonzeros /**< New number of nonzeros. */ 
);

/**
 * \brief Changes the number of nonzeros and reallocates corresponding arrays.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatChangeNumNonzeros(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< A matrix. */
  size_t newNumNonzeros /**< New number of nonzeros. */ 
);

/**
 * \brief Sorts the nonzeros of a double matrix by column in ascending order.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatSortNonzeros(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_DBLMAT* matrix  /**< A matrix. */
);

/**
 * \brief Sorts the nonzeros of an int matrix by column in ascending order.
 */

CMR_EXPORT
CMR_ERROR CMRintmatSortNonzeros(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_INTMAT* matrix  /**< A matrix. */
);

/**
 * \brief Sorts the nonzeros of a char matrix by column in ascending order.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatSortNonzeros(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix  /**< A matrix. */
);

/**
 * \brief Copies a double matrix to a newly allocated one.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatCopy(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,   /**< A matrix. */
  CMR_DBLMAT** presult  /**< Pointer for storing a copy of \p matrix. */
);

/**
 * \brief Copies an int matrix to a newly allocated one.
 */

CMR_EXPORT
CMR_ERROR CMRintmatCopy(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT* matrix,   /**< A matrix. */
  CMR_INTMAT** presult  /**< Pointer for storing a copy of \p matrix. */
);

/**
 * \brief Copies a char matrix to a newly allocated one.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatCopy(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< A matrix. */
  CMR_CHRMAT** presult  /**< Pointer for storing a copy of \p matrix. */
);

/**
 * \brief Creates the transpose of a double matrix.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatTranspose(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,   /**< A matrix. */
  CMR_DBLMAT** presult  /**< Pointer for storing the transpose of \p matrix. */
);

/**
 * \brief Creates the transpose of an int matrix.
 */

CMR_EXPORT
CMR_ERROR CMRintmatTranspose(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT* matrix,   /**< Given matrix. */
  CMR_INTMAT** presult  /**< Pointer for storing the transpose of \p matrix. */
);

/**
 * \brief Creates the transpose of a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatTranspose(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< Given matrix. */
  CMR_CHRMAT** presult  /**< Pointer for storing the transpose of \p matrix. */
);

/**
 * \brief Prints a double matrix in sparse format.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatPrintSparse(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_DBLMAT* matrix, /**< A matrix. */
  FILE* stream        /**< File stream to print to. */
);

/**
 * \brief Prints an int matrix in sparse format.
 */

CMR_EXPORT
CMR_ERROR CMRintmatPrintSparse(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_INTMAT* matrix, /**< A matrix. */
  FILE* stream        /**< File stream to print to. */
);

/**
 * \brief Prints a char matrix in sparse format.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatPrintSparse(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix, /**< A matrix. */
  FILE* stream        /**< File stream to print to. */
);

/**
 * \brief Prints a double matrix in dense format.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatPrintDense(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_DBLMAT* matrix, /**< A matrix. */
  FILE* stream,       /**< File stream to print to. */
  char zeroChar,      /**< Character to print for a zero. */
  bool header         /**< Whether to print row and column indices. */
);

/**
 * \brief Prints an int matrix in dense format.
 */

CMR_EXPORT
CMR_ERROR CMRintmatPrintDense(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_INTMAT* matrix, /**< A matrix. */
  FILE* stream,       /**< File stream to print to. */
  char zeroChar,      /**< Character to print for a zero. */
  bool header         /**< Whether to print row and column indices. */
);

/**
 * \brief Prints a char matrix in dense format.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatPrintDense(
  CMR* cmr,           /**< \ref CMR environment. */
  CMR_CHRMAT* matrix, /**< A matrix. */
  FILE* stream,       /**< File stream to print to. */
  char zeroChar,      /**< Character to print for a zero. */
  bool header         /**< Whether to print row and column indices. */
);

/**
 * \brief Reads a double matrix from a file \p stream in sparse format.
 * 
 * Returns \ref CMR_ERROR_INPUT in case of errors. In this case, *\p presult will be \c NULL.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatCreateFromSparseStream(
  CMR* cmr,             /**< \ref CMR environment. */
  FILE* stream,         /**< File stream to read from. */
  CMR_DBLMAT** presult  /**< Pointer for storing the matrix. */
);

/**
 * \brief Reads an int matrix from a file \p stream in sparse format.
 * 
 * Returns \ref CMR_ERROR_INPUT in case of errors. In this case, *\p presult will be \c NULL.
 */

CMR_EXPORT
CMR_ERROR CMRintmatCreateFromSparseStream(
  CMR* cmr,             /**< \ref CMR environment. */
  FILE* stream,         /**< File stream to read from. */
  CMR_INTMAT** presult  /**< Pointer for storing the matrix. */
);

/**
 * \brief Reads a char matrix from a file \p stream in sparse format.
 * 
 * Returns \ref CMR_ERROR_INPUT in case of errors. In this case, *\p presult will be \c NULL.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatCreateFromSparseStream(
  CMR* cmr,             /**< \ref CMR environment. */
  FILE* stream,         /**< File stream to read from. */
  CMR_CHRMAT** presult  /**< Pointer for storing the matrix. */
);

/**
 * \brief Reads a double matrix from a file \p stream in dense format.
 *
 * Returns \ref CMR_ERROR_INPUT in case of errors. In this case, *\p presult will be \c NULL.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatCreateFromDenseStream(
  CMR* cmr,             /**< \ref CMR environment. */
  FILE* stream,         /**< File stream to read from. */
  CMR_DBLMAT** presult  /**< Pointer for storing the matrix. */
);

/**
 * \brief Reads an int matrix from a file \p stream in dense format.
 *
 * Returns \ref CMR_ERROR_INPUT in case of errors. In this case, *\p presult will be \c NULL.
 */

CMR_EXPORT
CMR_ERROR CMRintmatCreateFromDenseStream(
  CMR* cmr,             /**< \ref CMR environment. */
  FILE* stream,         /**< File stream to read from. */
  CMR_INTMAT** presult  /**< Pointer for storing the matrix. */
);

/**
 * \brief Reads a char matrix from a file \p stream in dense format.
 *
 * Returns \ref CMR_ERROR_INPUT in case of errors. In this case, *\p presult will be \c NULL.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatCreateFromDenseStream(
  CMR* cmr,             /**< \ref CMR environment. */
  FILE* stream,         /**< File stream to read from. */
  CMR_CHRMAT** presult  /**< Pointer for storing the matrix. */
);

/**
 * \brief Checks whether two double matrices are equal.
 */

CMR_EXPORT
bool CMRdblmatCheckEqual(
  CMR_DBLMAT* matrix1,  /**< First matrix. */
  CMR_DBLMAT* matrix2   /**< Second matrix. */
);

/**
 * \brief Checks whether two int matrices are equal.
 */

CMR_EXPORT
bool CMRintmatCheckEqual(
  CMR_INTMAT* matrix1, /**< First matrix. */
  CMR_INTMAT* matrix2  /**< Second matrix. */
);

/**
 * \brief Checks whether two char matrices are equal.
 */

CMR_EXPORT
bool CMRchrmatCheckEqual(
  CMR_CHRMAT* matrix1,  /**< First matrix. */
  CMR_CHRMAT* matrix2   /**< Second matrix. */
);

/**
 * \brief Checks whether two double matrices are transposes of each other.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatCheckTranspose(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT* matrix1,  /**< First matrix. */
  CMR_DBLMAT* matrix2,  /**< Second matrix. */
  bool* pareTranspose   /**< Pointer for storing whether \p matrix1 and \p matrix2 are tranposes of each other. */
);

/**
 * \brief Checks whether two int matrices are transposes of each other.
 */

CMR_EXPORT
CMR_ERROR CMRintmatCheckTranspose(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT* matrix1, /**< First matrix */
  CMR_INTMAT* matrix2, /**< Second matrix */
  bool* pareTranspose   /**< Pointer for storing whether \p matrix1 and \p matrix2 are tranposes of each other. */
);

/**
 * \brief Checks whether two char matrices are transposes of each other.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatCheckTranspose(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix1, /**< First matrix */
  CMR_CHRMAT* matrix2, /**< Second matrix */
  bool* pareTranspose   /**< Pointer for storing whether \p matrix1 and \p matrix2 are tranposes of each other. */
);

/**
 * \brief Checks a double matrix for consistency.
 *
 * Checks whether the entries of a row are sorted by column in ascending order.
 * Checks for duplicate entries.
 * Checks for zero entries.
 *
 * \returns Explanation of inconsistency, or \c NULL if consistent.
 */

CMR_EXPORT
char* CMRdblmatConsistency(
  CMR_DBLMAT* matrix /**< A matrix. */
);

/**
 * \brief Checks an int matrix for consistency.
 *
 * Checks whether the entries of a row are sorted by column in ascending order.
 * Checks for duplicate entries.
 * Checks for zero entries.
 *
 * \returns Explanation of inconsistency, or \c NULL if consistent.
 */

CMR_EXPORT
char* CMRintmatConsistency(
  CMR_INTMAT* matrix /**< A matrix. */
);

/**
 * \brief Checks a char matrix for consistency.
 *
 * Checks whether the entries of a row are sorted by column in ascending order.
 * Checks for duplicate entries.
 * Checks for zero entries.
 *
 * \returns Explanation of inconsistency, or \c NULL if consistent.
 */

CMR_EXPORT
char* CMRchrmatConsistency(
  CMR_CHRMAT* matrix /**< A matrix. */
);

/**
 * \brief Creates a submatrix of a double matrix as an explicit matrix.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatFilterSubmat(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,     /**< A matrix */
  CMR_SUBMAT* submatrix,  /**< A submatrix of \p matrix. */
  CMR_DBLMAT** presult    /**< Pointer for storing the resulting matrix. */
);

/**
 * \brief Creates a submatrix of an int matrix as an explicit matrix.
 */

CMR_EXPORT
CMR_ERROR CMRintmatFilterSubmat(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_INTMAT* matrix,     /**< A matrix */
  CMR_SUBMAT* submatrix,  /**< A submatrix of \p matrix. */
  CMR_INTMAT** presult    /**< Pointer for storing the resulting matrix. */
);

/**
 * \brief Creates a submatrix of a char matrix as an explicit matrix.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatFilterSubmat(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< A matrix */
  CMR_SUBMAT* submatrix,  /**< A submatrix of \p matrix. */
  CMR_CHRMAT** presult    /**< Pointer for storing the resulting double matrix. */
);

/**
 * \brief Checks if a double matrix has only entries in \f$ \{0,1\} \f$ with absolute error tolerance \p epsilon.
 */

CMR_EXPORT
bool CMRdblmatIsBinary(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,     /**< A matrix. */
  double epsilon,         /**< Absolute error tolerance. */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a non-binary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Checks if an int matrix has only entries in \f$ \{0,1\} \f$.
 */

CMR_EXPORT
bool CMRintmatIsBinary(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_INTMAT* matrix,     /**< A matrix. */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a non-binary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Checks if a char matrix has only entries in \f$ \{0,1\} \f$.
 */

CMR_EXPORT
bool CMRchrmatIsBinary(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< A matrix. */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a non-binary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Checks if a double matrix has only entries in \f$ \{-1,0,+1\} \f$ with absolute error tolerance \p epsilon.
 */

CMR_EXPORT
bool CMRdblmatIsTernary(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,     /**< A matrix. */
  double epsilon,         /**< Absolute error tolerance. */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a non-ternary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Checks if an int matrix has only entries in \f$ \{-1,0,+1\} \f$.
 */

CMR_EXPORT
bool CMRintmatIsTernary(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_INTMAT* matrix,     /**< A matrix. */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a non-ternary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Checks if a double matrix has only entries in \f$ \{-1,0,+1\} \f$.
 */

CMR_EXPORT
bool CMRchrmatIsTernary(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< A matrix. */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a non-ternary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Creates the (binary) support matrix of a double matrix as a char matrix with absolute error tolerance \p epsilon.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatSupport(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,   /**< Double matrix */
  double epsilon,       /**< Absolute error tolerance. */
  CMR_CHRMAT** presult  /**< Pointer for storing the support matrix of \p matrix. */
);

/**
 * \brief Creates the (binary) support matrix of an int matrix as a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRintmatSupport(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT* matrix,   /**< Double matrix */
  CMR_CHRMAT** presult  /**< Pointer for storing the support matrix of \p matrix. */
);

/**
 * \brief Creates the (binary) support matrix of a char matrix as a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatSupport(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< Double matrix */
  CMR_CHRMAT** presult  /**< Pointer for storing the support matrix of \p matrix. */
);

/**
 * \brief Creates the (ternary) signed support matrix of a double matrix as a char matrix with absolute error tolerance
 *        \p epsilon.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatSignedSupport(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,   /**< Double matrix */
  double epsilon,       /**< Absolute error tolerance. */
  CMR_CHRMAT** presult  /**< Pointer for storing the support matrix of \p matrix. */
);

/**
 * \brief Creates the (ternary) signed support matrix of an int matrix as a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRintmatSignedSupport(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT* matrix,   /**< Double matrix */
  CMR_CHRMAT** presult  /**< Pointer for storing the signed support matrix of \p matrix. */
);

/**
 * \brief Creates the (ternary) signed support matrix of a char matrix as a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatSignedSupport(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,   /**< Double matrix */
  CMR_CHRMAT** presult  /**< Pointer for storing the signed support matrix of \p matrix. */
);

/**
 * \brief Converts an int matrix to a char matrix.
 *
 * \returns \ref CMR_ERROR_INPUT in case of overflow.
 */

CMR_EXPORT
CMR_ERROR CMRintmatToChr(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT* matrix,   /**< Input matrix. */
  CMR_CHRMAT** presult  /**< Pointer for storing the output matrix. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_MATRIX_H */
