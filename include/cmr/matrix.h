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
  /**
   * \brief Number of rows
   */
  size_t numRows;
  /**
   * \brief Array with row indices
   */
  size_t* rows;
  /**
   * \brief Number of columns
   */
  size_t numColumns;
  /**
   * \brief Array with column indices
   */
  size_t* columns;
} CMR_SUBMAT;

/**
 * \brief Row-wise representation of sparse double matrix.
 * 
 * The columns and values of all nonzeros are stored in \ref entryColumns and \ref entryValues,
 * respectively.
 * Those of row \c r are stored from \ref rowStarts[r] until (but not including)
 * \ref rowStarts[r+1]. The last row is an exception, since \ref rowStarts[\ref numRows] need not
 * be defined.
 * For convenience, one may store this additional entry.
 * In particular \ref CMRdblmatCreate allocates sufficient space for it.
 * However, all public methods use \ref numRows to determine the last row's number of nonzeros via
 * \ref numNonzeros.
 */

typedef struct
{
  int numRows;          /**< \brief Number of rows. */
  int numColumns;       /**< \brief Number of columns. */
  int numNonzeros;      /**< \brief Number of and memory allocated for nonzeros. */
  int* rowStarts;       /**< \brief Array mapping each row to the index of its first entry. */
  int* entryColumns;    /**< \brief Array mapping each entry to its column.*/
  double* entryValues;  /**< \brief Array mapping each entry to its value. */
} CMR_DBLMAT;

/**
 * \brief Row-wise representation of sparse int matrix.
 * 
 * The columns and values of all nonzeros are stored in \ref entryColumns and \ref entryValues,
 * respectively.
 * Those of row \c r are stored from \ref rowStarts[r] until (but not including)
 * \ref rowStarts[r+1]. The last row is an exception, since \ref rowStarts[\ref numRows] need not
 * be defined.
 * For convenience, one may store this additional entry.
 * In particular \ref CMRintmatCreate allocates sufficient space for it.
 * However, all public methods use \ref numRows to determine the last row's number of nonzeros via
 * \ref numNonzeros.
 */

typedef struct
{
  int numRows;        /**< \brief Number of rows. */
  int numColumns;     /**< \brief Number of columns. */
  int numNonzeros;    /**< \brief Number of and memory allocated for nonzeros. */
  int* rowStarts;     /**< \brief Array mapping each row to the index of its first entry. */
  int* entryColumns;  /**< \brief Array mapping each entry to its column.*/
  int* entryValues;   /**< \brief Array mapping each entry to its value. */
} CMR_INTMAT;

/**
 * \brief Row-wise representation of sparse char matrix.
 * 
 * The columns and values of all nonzeros are stored in \ref entryColumns and \ref entryValues,
 * respectively.
 * Those of row \c r are stored from \ref rowStarts[r] until (but not including)
 * \ref rowStarts[r+1]. The last row is an exception, since \ref rowStarts[\ref numRows] need not
 * be defined.
 * For convenience, one may store this additional entry.
 * In particular \ref CMRchrmatCreate allocates sufficient space for it.
 * However, all public methods use \ref numRows to determine the last row's number of nonzeros via
 * \ref numNonzeros.
 */

typedef struct
{
  int numRows;        /**< \brief Number of rows. */
  int numColumns;     /**< \brief Number of columns. */
  int numNonzeros;    /**< \brief Number of and memory allocated for nonzeros. */
  int* rowStarts;     /**< \brief Array mapping each row to the index of its first entry. */
  int* entryColumns;  /**< \brief Array mapping each entry to its column.*/
  char* entryValues;  /**< \brief Array mapping each entry to its value. */
} CMR_CHRMAT;



















/**
 * \brief Creates a double matrix of size \p numRows times \p numColumns with \p numNonzeros
 *        nonzeros. The row starts and entries are allocated but not initialized.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatCreate(
  CMR* cmr,              /**< \ref CMR environment. */
  CMR_DBLMAT** matrix,  /**< Pointer for storing the created matrix. */
  int numRows,         /**< Number of rows. */
  int numColumns,      /**< Number of columns. */
  int numNonzeros      /**< Number of nonzeros. */
);





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
 * \brief Creates a 1x1 submatrix.
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
 * \brief Frees the memory of a double matrix.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatFree(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT** pmatrix  /**< Pointer to matrix. */
);

/**
 * \brief Changes the number of nonzeros and reallocates corresponding arrays.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatChangeNumNonzeros(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,  /**< Given matrix. */
  int newNumNonzeros  /**< New number of nonzeros. */ 
);

/**
 * \brief Copies a double matrix to a newly allocated one.
 * 
 * Allocates *\p result and copies \p matrix there.
 */
CMR_EXPORT
CMR_ERROR CMRdblmatCopy(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,  /**< Given matrix. */
  CMR_DBLMAT** result  /**< Pointer to store a copy of \p matrix. */
);

/**
 * \brief Creates the transpose of a double matrix.
 */
CMR_EXPORT
CMR_ERROR CMRdblmatTranspose(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,  /**< Given matrix. */
  CMR_DBLMAT** result  /**< Pointer to store the transpose of \p matrix. */
);

/**
 * \brief Prints a double matrix.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatPrintSparse(
  FILE* stream,       /**< File stream to print to. */
  CMR_DBLMAT* matrix   /**< Double matrix. */
);


/**
 * \brief Prints a double matrix.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatPrintDense(
  FILE* stream,       /**< File stream to print to. */
  CMR_DBLMAT* matrix,  /**< Double matrix. */
  char zeroChar,      /**< Character to print for a zero. */
  bool header         /**< Whether to print row and column indices. */
);

/**
 * \brief Reads a sparse double matrix from a file \p stream.
 * 
 * Zero entries are ignored, and multiple occurences of (row,column) pairs are considered as errors.
 * Returns \ref CMR_ERROR_INPUT in case of errors. In this case, *\p pmatrix will be \c NULL.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatCreateFromSparseStream(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_DBLMAT** pmatrix,  /**< Pointer for storing the matrix. */
  FILE* stream          /**< File stream to read from. */
);

/**
 * \brief Reads a densely stored double matrix from a file \p stream.
 */

CMR_EXPORT
CMR_ERROR CMRdblmatCreateFromDenseStream(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_DBLMAT** pmatrix,  /**< Pointer for storing the matrix. */
  FILE* stream          /**< File stream to read from. */
);

/**
 * \brief Checks whether two double matrices are equal.
 */

CMR_EXPORT
bool CMRdblmatCheckEqual(
  CMR_DBLMAT* matrix1,  /**< First matrix */
  CMR_DBLMAT* matrix2   /**< Second matrix */
);

/**
 * \brief Checks whether two double matrices are transposes of each other.
 */

CMR_EXPORT
bool CMRdblmatCheckTranspose(
  CMR_DBLMAT* matrix1,  /**< First matrix */
  CMR_DBLMAT* matrix2   /**< Second matrix */
);

/**
 * \brief Checks whether double matrix has each row sorted by minor.
 */

CMR_EXPORT
bool CMRdblmatCheckSorted(
  CMR_DBLMAT* matrix /**< Double matrix */
);

/**
 * \brief Creates a submatrix of a double matrix as an explicit matrix.
 */
CMR_EXPORT
CMR_ERROR CMRdblmatFilterSubmat(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,    /**< Given matrix */
  CMR_SUBMAT* submatrix, /**< Specified submatrix */
  CMR_DBLMAT** result    /**< Pointer for storing the resulting double matrix. */
);

/**
 * \brief Checks if double matrix has only entries in {0, 1} with absolute error tolerance \p epsilon.
 */

CMR_EXPORT
bool CMRisBinaryDbl(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,      /**< Double matrix */
  double epsilon,         /**< Absolute error tolerance */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a non-binary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Checks if double matrix has only entries in {-1, 0, +1} with absolute error tolerance \p epsilon.
 */

CMR_EXPORT
bool CMRisTernaryDbl(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,      /**< Double matrix */
  double epsilon,         /**< Absolute error tolerance */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a non-ternary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Creates the support matrix of a double \p matrix as a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRsupportDbl(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,    /**< Double matrix */
  double epsilon,       /**< Absolute error tolerance */
  CMR_CHRMAT** psupport  /**< Pointer for storing the support matrix of \p matrix. */
);

/**
 * \brief Creates the signed support matrix of a double \p matrix as a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRsignedSupportDbl(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_DBLMAT* matrix,    /**< Double matrix */
  double epsilon,       /**< Absolute error tolerance */
  CMR_CHRMAT** psupport  /**< Pointer for storing the support matrix of \p matrix. */
);











/**
 * \brief Creates a double matrix of size \p numRows times \p numColumns with \p numNonzeros
 *        nonzeros. The row starts and entries are allocated but not initialized.
 */

CMR_EXPORT
CMR_ERROR CMRintmatCreate(
  CMR* cmr,              /**< \ref CMR environment. */
  CMR_INTMAT** matrix,  /**< Pointer for storing the created matrix. */
  int numRows,         /**< Number of rows. */
  int numColumns,      /**< Number of columns. */
  int numNonzeros      /**< Number of nonzeros. */
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
 * \brief Changes the number of nonzeros and reallocates corresponding arrays.
 */

CMR_EXPORT
CMR_ERROR CMRintmatChangeNumNonzeros(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT* matrix,  /**< Given matrix. */
  int newNumNonzeros  /**< New number of nonzeros. */ 
);

/**
 * \brief Copies an int matrix to a newly allocated one.
 * 
 * Allocates *\p result and copies \p matrix there.
 */
CMR_EXPORT
CMR_ERROR CMRintmatCopy(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT* matrix,  /**< Given matrix. */
  CMR_INTMAT** result  /**< Pointer to store a copy of \p matrix. */
);

/**
 * \brief Creates the transpose of an int matrix.
 */
CMR_EXPORT
CMR_ERROR CMRintmatTranspose(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_INTMAT* matrix,  /**< Given matrix. */
  CMR_INTMAT** result  /**< Pointer to store the transpose of \p matrix. */
);

/**
 * \brief Prints an int matrix.
 */

CMR_EXPORT
CMR_ERROR CMRintmatPrintSparse(
  FILE* stream,       /**< File stream to print to. */
  CMR_INTMAT* matrix   /**< Int matrix. */
);

/**
 * \brief Prints an int matrix.
 */

CMR_EXPORT
CMR_ERROR CMRintmatPrintDense(
  FILE* stream,       /**< File stream to print to. */
  CMR_INTMAT* matrix,  /**< Int matrix. */
  char zeroChar,      /**< Character to print for a zero. */
  bool header         /**< Whether to print row and column indices. */
);

/**
 * \brief Reads a sparse int matrix from a file \p stream.
 *  
 * Zero entries are ignored, and multiple occurences of (row,column) pairs are considered as errors.
 * Returns \ref CMR_ERROR_INPUT in case of errors. In this case, *\p pmatrix will be \c NULL.
 */

CMR_EXPORT
CMR_ERROR CMRintmatCreateFromSparseStream(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_INTMAT** pmatrix,  /**< Pointer for storing the matrix. */
  FILE* stream          /**< File stream to read from. */
);

/**
 * \brief Reads a densely stored int matrix from a file \p stream.
 */

CMR_EXPORT
CMR_ERROR CMRintmatCreateFromDenseStream(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_INTMAT** pmatrix,  /**< Pointer for storing the matrix. */
  FILE* stream          /**< File stream to read from. */
);

/**
 * \brief Checks whether two int matrices are equal.
 */

CMR_EXPORT
bool CMRintmatCheckEqual(
  CMR_INTMAT* matrix1, /**< First matrix */
  CMR_INTMAT* matrix2  /**< Second matrix */
);

/**
 * \brief Checks whether two int matrices are transposes of each other.
 */

CMR_EXPORT
bool CMRintmatCheckTranspose(
  CMR_INTMAT* matrix1, /**< First matrix */
  CMR_INTMAT* matrix2  /**< Second matrix */
);

/**
 * \brief Checks whether int matrix has each row sorted by minor.
 */

CMR_EXPORT
bool CMRintmatCheckSorted(
  CMR_INTMAT* matrix /**< Int matrix */
);

/**
 * \brief Creates a submatrix of an int matrix as an explicit matrix.
 */
CMR_EXPORT
CMR_ERROR CMRintmatFilterSubmat(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_INTMAT* matrix,    /**< Given matrix */
  CMR_SUBMAT* submatrix, /**< Specified submatrix */
  CMR_INTMAT** result    /**< Pointer for storing the resulting int matrix. */
);

/**
 * \brief Checks if int matrix has only entries in {0, 1}.
 */

CMR_EXPORT
bool CMRisBinaryInt(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_INTMAT* matrix,      /**< Int matrix */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a non-binary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Checks if int matrix has only entries in {-1, 0, +1}.
 */

CMR_EXPORT
bool CMRisTernaryInt(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_INTMAT* matrix,      /**< Int matrix */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a non-ternary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Creates the support matrix of an int \p matrix as a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRsupportInt(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_INTMAT* matrix,    /**< Int matrix */
  CMR_CHRMAT** psupport  /**< Pointer for storing the support matrix of \p matrix. */
);

/**
 * \brief Creates the signed support matrix of an int \p matrix as a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRsignedSupportInt(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_INTMAT* matrix,    /**< Int matrix */
  CMR_CHRMAT** psupport  /**< Pointer for storing the support matrix of \p matrix. */
);

















/**
 * \brief Creates a char matrix of size \p numRows times \p numColumns with \p numNonzeros
 *        nonzeros. The row starts and entries are allocated but not initialized.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatCreate(
  CMR* cmr,              /**< \ref CMR environment. */
  CMR_CHRMAT** matrix,  /**< Pointer for storing the created matrix. */
  int numRows,         /**< Number of rows. */
  int numColumns,      /**< Number of columns. */
  int numNonzeros      /**< Number of nonzeros. */
);

/**
 * \brief Frees the memory of an int matrix.
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
CMR_ERROR CMRchrmatChangeNumNonzeros(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,  /**< Given matrix. */
  int newNumNonzeros  /**< New number of nonzeros. */ 
);

/**
 * \brief Copies an int matrix to a newly allocated one.
 * 
 * Allocates *\p result and copies \p matrix there.
 */
CMR_EXPORT
CMR_ERROR CMRchrmatCopy(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,  /**< Given matrix. */
  CMR_CHRMAT** result  /**< Pointer to store a copy of \p matrix. */
);

/**
 * \brief Creates the transpose of an int matrix.
 */
CMR_EXPORT
CMR_ERROR CMRchrmatTranspose(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,  /**< Given matrix. */
  CMR_CHRMAT** result  /**< Pointer to store the transpose of \p matrix. */
);



/**
 * \brief Prints a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatPrintSparse(
  FILE* stream,       /**< File stream to print to. */
  CMR_CHRMAT* matrix   /**< Char matrix. */
);



/**
 * \brief Prints a char matrix.
 */

CMR_ERROR CMRchrmatPrintDense(
  CMR* cmr,           /**< \ref CMR environment. */
  FILE* stream,       /**< File stream to print to. */
  CMR_CHRMAT* matrix, /**< Char matrix. */
  char zeroChar,      /**< Character to print for a zero. */
  bool header         /**< Whether to print row and column indices. */
);



/**
 * \brief Reads a sparse char matrix from a file \p stream.
 *  
 * Zero entries are ignored, and multiple occurences of (row,column) pairs are considered as errors.
 * Returns \ref CMR_ERROR_INPUT in case of errors. In this case, *\p pmatrix will be \c NULL.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatCreateFromSparseStream(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT** pmatrix,  /**< Pointer for storing the matrix. */
  FILE* stream          /**< File stream to read from. */
);


/**
 * \brief Reads a densely stored char matrix from a file \p stream.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatCreateFromDenseStream(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT** pmatrix,  /**< Pointer for storing the matrix. */
  FILE* stream          /**< File stream to read from. */
);

/**
 * \brief Checks whether two char matrices are equal.
 */

CMR_EXPORT
bool CMRchrmatCheckEqual(
  CMR_CHRMAT* matrix1,  /**< First matrix */
  CMR_CHRMAT* matrix2   /**< Second matrix */
);

/**
 * \brief Checks whether two char matrices are transposes of each other.
 */

CMR_EXPORT
bool CMRchrmatCheckTranspose(
  CMR_CHRMAT* matrix1, /**< First matrix */
  CMR_CHRMAT* matrix2  /**< Second matrix */
);

/**
 * \brief Checks whether char matrix has each row sorted by minor.
 */

CMR_EXPORT
bool CMRchrmatCheckSorted(
  CMR_CHRMAT* matrix /**< Char matrix */
);

/**
 * \brief Creates a submatrix of a char matrix as an explicit matrix.
 */

CMR_EXPORT
CMR_ERROR CMRchrmatFilterSubmat(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Given matrix */
  CMR_SUBMAT* submatrix,  /**< Specified submatrix */
  CMR_CHRMAT** presult    /**< Pointer for storing the resulting char matrix. */
);

/**
 * \brief Checks if matrix has only entries in {0, 1}.
 */

CMR_EXPORT
bool CMRisBinaryChr(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,      /**< Char matrix */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a non-binary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Checks if char matrix has only entries in {-1, 0, +1}.
 */

CMR_EXPORT
bool CMRisTernaryChr(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,      /**< Char matrix */
  CMR_SUBMAT** psubmatrix  /**< Pointer for storing a non-ternary entry as a submatrix (may be \c NULL). */
);

/**
 * \brief Creates the support matrix of a char \p matrix as a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRsupportChr(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,    /**< Char matrix */
  CMR_CHRMAT** psupport  /**< Pointer for storing the support matrix of \p matrix. */
);

/**
 * \brief Creates the signed support matrix of a char \p matrix as a char matrix.
 */

CMR_EXPORT
CMR_ERROR CMRsignedSupportChr(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,    /**< Char matrix */
  CMR_CHRMAT** psupport  /**< Pointer for storing the support matrix of \p matrix. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_MATRIX_H */
