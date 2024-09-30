#ifndef CMR_SEPARATION_H
#define CMR_SEPARATION_H

#include <cmr/env.h>
#include <cmr/matrix.h>
#include <cmr/element.h>

#include <assert.h>
#include <stdint.h>

/**
 * \file separation.h
 *
 * \author Matthias Walter
 *
 * \brief Data structures for k-separations and k-sums.
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef enum
{
  CMR_SEPA_FIRST = 0,
    /**< This row/column belongs to the first child. */
  CMR_SEPA_SECOND = 1,
    /**< This row/column belongs to the second child. */
  CMR_SEPA_FLAG_RANK1 = 2,
    /**< This bit flag indicates that the row/column also belongs to the other child. */
  CMR_SEPA_FLAG_RANK2 = 4,
    /**< This bit flag indicates that the row/column also belongs to the other child and is independent of the one
     *   from CMR_SEPA_FLAG_RANK1. */

  CMR_SEPA_MASK_CHILD = 1,
    /**< Bit mask for the ownership. */
  CMR_SEPA_MASK_EXTRA = CMR_SEPA_FLAG_RANK1 | CMR_SEPA_FLAG_RANK2
    /**< Bit mask extra belongings. */
} CMR_SEPA_FLAGS;

typedef enum
{
  CMR_SEPA_TYPE_TWO = 2,
    /**< 2-separation whose bottom-left part has rank 1. */
  CMR_SEPA_TYPE_THREE_DISTRIBUTED_RANKS = 3,
    /**< 3-separation with distributed ranks. */
  CMR_SEPA_TYPE_THREE_CONCENTRATED_RANK = 4
    /**< 3-separation whose botom-left part has rank 2. */
} CMR_SEPA_TYPE;

typedef struct
{
  size_t numRows;     /**< \brief Number of rows of the matrix. */
  size_t numColumns;  /**< \brief Number of columns of the matrix. */
  int* rowsFlags;     /**< \brief Array with each row's flags; logical OR of CMR_SEPA_TYPE. */
  int* columnsFlags;  /**< \brief Array with each column's flags; logical OR of CMR_SEPA_TYPE. */
  CMR_SEPA_TYPE type; /**< \brief Type of separation. */
} CMR_SEPA;

/**
 * \brief Creates a 2- or 3-separation.
 *
 * Only the memory is allocated. The rowsFlags and columnsFlags arrays must be filled properly.
 */

CMR_EXPORT
CMR_ERROR CMRsepaCreate(
  CMR* cmr,           /**< \ref CMR environment. */
  size_t numRows,     /**< Number of rows. */
  size_t numColumns,  /**< Number of columns. */
  CMR_SEPA** psepa    /**< Pointer for storing the created separation. */
);

/**
 * \brief Frees a separation.
 */

CMR_EXPORT
CMR_ERROR CMRsepaFree(
  CMR* cmr,         /**< \ref CMR environment. */
  CMR_SEPA** psepa  /**< Pointer to separation. */
);

/**
 * \brief Computes the sizes of the top-left and bottom-right parts.
 */

CMR_EXPORT
CMR_ERROR CMRsepaComputeSizes(
  CMR_SEPA* sepa,                 /**< Separation. */
  size_t* pnumRowsTopLeft,        /**< Pointer for storing the number of rows of the top-left part. */
  size_t* pnumColumnsTopLeft,     /**< Pointer for storing the number of columns of the top-left part. */
  size_t* pnumRowsBottomRight,    /**< Pointer for storing the number of rows of the bottom-right part. */
  size_t* pnumColumnsBottomRight  /**< Pointer for storing the number of columns of the bottom-right part. */
);

/**
 * \brief Scans the support of \p matrix to compute all representative rows/columns for \p sepa and sets the type.
 *
 * Assumes that the sum of the ranks of the off-diagonal blocks is at most 2. Potentially swaps parts to ensure that
 * the rank of the bottom-left submatrix is at least that of the top-right submatrix.
 * Sets the rowsFlags and columnsFlags attributes accordingly.
 */

CMR_EXPORT
CMR_ERROR CMRsepaFindBinaryRepresentatives(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEPA* sepa,         /**< Separation. */
  CMR_CHRMAT* matrix,     /**< Matrix. */
  CMR_CHRMAT* transpose,  /**< Transpose of \p matrix. */
  bool* pswapped,         /**< Pointer for storing whether parts were swapped (may be \c NULL). */
  CMR_SUBMAT** pviolator  /**< Pointer for storing a violator submatrix if the ternary rank differs (may be \c NULL). */
);

/**
 * \brief Scans the support of \p submatrix of \p matrix to compute all representative rows/columns for \p sepa and
 *        sets the type.
 *
 * Assumes that the sum of the ranks of the off-diagonal blocks is at most 2. Potentially swaps parts to ensure that
 * the rank of the bottom-left submatrix is at least that of the top-right submatrix.
 * Sets the rowsFlags and columnsFlags attributes accordingly.
 */

CMR_EXPORT
CMR_ERROR CMRsepaFindBinaryRepresentativesSubmatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEPA* sepa,         /**< Separation. */
  CMR_CHRMAT* matrix,     /**< Matrix. */
  CMR_CHRMAT* transpose,  /**< Transpose of \p matrix. */
  CMR_SUBMAT* submatrix,  /**< Submatrix of \p matrix. */
  bool* pswapped,         /**< Pointer for storing whether parts were swapped (may be \c NULL). */
  CMR_SUBMAT** pviolator  /**< Pointer for storing a violator submatrix if the ternary rank differs (may be \c NULL). */
);

/**
 * \brief Returns representative rows/columns of the low-rank submatrices.
 */

CMR_EXPORT
CMR_ERROR CMRsepaGetRepresentatives(
  CMR_SEPA* sepa,           /**< Separation. */
  size_t reprRows[2][3],    /**< Array mapping child to arrays of (at most 3) different representative rows. */
  size_t reprColumns[2][3]  /**< Array mapping child to arrays of (at most 3) different representative columns. */
);

/**
 * \brief Creates mappings from rows/columns to those of \p part; also maps up to 3 representative rows/columns.
 */

CMR_EXPORT
CMR_ERROR CMRsepaGetProjection(
  CMR_SEPA* sepa,         /**< Separation. */
  size_t part,            /**< Part to project. */
  size_t* rowsToPart,     /**< Array for storing the mapping from rows to those of \p part. Must be large enough. */
  size_t* columnsToPart,  /**< Array for storing the mapping from columns to those of \p part. Must be large enough. */
  size_t* pnumPartRows,   /**< Pointer for storing the number of rows of \p part (excluding representatives). */
  size_t* pnumPartColumns /**< Pointer for storing the number of columns of \p part (excluding representatives). */
);

/**
 * \brief Checks for a given matrix whether the binary k-separation is also a ternary one.
 *
 * Checks, for a ternary input matrix \f$ M \f$ and a k-separation (\f$ k \in \{2,3\} \f$) of the (binary) support
 * matrix of \f$ M \f$, whether it is also a k-separation of \f$ M \f$ itself. The result is stored in \p *pisTernary.
 *
 * If the check fails, a violating 2-by-2 submatrix is returned.
 */

CMR_EXPORT
CMR_ERROR CMRsepaCheckTernary(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEPA* sepa,         /**< Separation. */
  CMR_CHRMAT* matrix,     /**< Matrix. */
  bool* pisTernary,       /**< Pointer for storing whether the check passed. */
  CMR_SUBMAT** pviolator  /**< Pointer for storing a violator submatrix (may be \c NULL). */
);

/**
 * \brief Checks for a submatrix of a given matrix whether the binary k-separation is also a ternary one.
 *
 * Checks, for a ternary input matrix \f$ M \f$ and a k-separation (\f$ k \in \{2,3\} \f$) of the (binary) support
 * matrix of \f$ M \f$, whether it is also a k-separation of \f$ M \f$ itself. The result is stored in \p *pisTernary.
 *
 * If the check fails, a certifying submatrix is returned in \p *pviolator. Its row/column indices refer to
 * \p submatrix.
 */

CMR_EXPORT
CMR_ERROR CMRsepaCheckTernarySubmatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEPA* sepa,         /**< Separation. */
  CMR_CHRMAT* matrix,     /**< Matrix. */
  CMR_SUBMAT* submatrix,  /**< Submatrix to consider. */
  bool* pisTernary,       /**< Pointer for storing whether the check passed. */
  CMR_SUBMAT** pviolator  /**< Pointer for storing a violator submatrix (may be \c NULL). */
);

/**
 * \brief Composes the 1-sum of the several \p matrices.
 *
 * Let \f$ A_1, A_2, \dotsc, A_k \f$ denote the matrices given by the array \p matrices.
 * Their 1-sum is the matrix
 * \f[
 *   B := \begin{bmatrix}
 *     A_1        & \mathbb{O} & \dots      & \mathbb{O} \\
 *     \mathbb{O} & A_2        &            & \vdots  \\
 *     \vdots     &            & \ddots     & \mathbb{O} \\
 *     \mathbb{O} & \dots      & \mathbb{O} & A_k
 *   \end{bmatrix}.
 * \f]
 * The resulting matrix \f$ B \f$ is created and stored in \p *presult.
 *
 * \see \ref CMRdecomposeBlocks for a decomposition of a given matrix \f$ B \f$ into \f$ A_i \f$.
 */

CMR_EXPORT
CMR_ERROR CMRoneSumCompose(
  CMR* cmr,               /**< \ref CMR environment. */
  size_t numMatrices,     /**< Number \f$ k \f$ of matrices in the sum. */
  CMR_CHRMAT** matrices,  /**< First matrix. */
  CMR_CHRMAT** presult    /**< Pointer for storing the result. */
);

/**
 * \brief Composes the 2-sum of the two matrices \p first and \p second with connecting elements \p firstMarker and
 *        \p secondMarker.
 *
 * If \p firstMarker indexes a row of \p first with row vector \f$ c^{\textsf{T}} \f$ then \p secondMarker must index a
 * column of \p second with column vector \f$ d \f$. In this case, let \f$ M_1 \f$ and \f$ M_2 \f$ denote the matrices
 * given by \p first and \p second, reordered such that \f$ M_1 = \begin{bmatrix} A \\ c^{\textsf{T}} \end{bmatrix} \f$
 * and \f$ M_2 = \begin{bmatrix} d & D \end{bmatrix} \f$ holds. Then the 2-sum is the matrix
 * \f[
 *   M := \begin{bmatrix}
 *     A & \mathbb{O} \\
 *     d c^{\textsf{T}} & D
 *   \end{bmatrix}.
 * \f]
 * Otherwise, \p firstMarker must index a column of \p first with column vector \f$ a \f$ and \p secondMarker indexes a
 * row of \p second with row vector \f$ b^{\textsf{T}} \f$. Let \f$ M_1 \f$ and \f$ M_2 \f$ denote the matrices given by
 * \p first and \p second, reordered such that \f$ M_1 = \begin{bmatrix} A & a \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} b^{\textsf{T}} \\ D \end{bmatrix} \f$ holds. Then the 2-sum is the matrix
 * \f[
 *   M := \begin{bmatrix}
 *     A & a b^{\textsf{T}} \\
 *     \mathbb{O} & D
 *   \end{bmatrix}.
 * \f]
 * The calculations are done modulo \p characteristic, where the value \f$ 3 \f$ yields numbers from
 * \f$ \{-1,0,+1\} \f$.
 *
 * The resulting matrix \f$ M \f$ is created and stored in \p *presult.
 */

CMR_EXPORT
CMR_ERROR CMRtwoSumCompose(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_CHRMAT* first,        /**< First matrix. */
  CMR_CHRMAT* second,       /**< Second matrix. */
  CMR_ELEMENT firstMarker,  /**< Marker element of first matrix. */
  CMR_ELEMENT secondMarker, /**< Marker element of second matrix. */
  int8_t characteristic,    /**< Field characteristic. */
  CMR_CHRMAT** presult      /**< Pointer for storing the result. */
);

/**
 * \brief Decomposes \p matrix as a 2-sum according to the 2-separation \p sepa and computing the first component.
 *
 * The input \p matrix \f$ M \f$ must have a 2-separation that is given by \p sepa, i.e., it can be reordered to look
 * like \f$ M = \begin{bmatrix} A & B \\ C & D \end{bmatrix} \f$, where \f$ \rank(B) + \rank(C) = 1 \f$.
 * If \f$ \rank(B) = \mathbb{O} \f$ then the two components of the 2-sum are matrices
 * \f$ M_1 = \begin{bmatrix} A \\ c^{\textsf{T}} \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} d & D \end{bmatrix} \f$ such that \f$ C = d c^{\textsf{T}} \f$ holds.
 * Otherwise, the two parts of the 2-sum are matrices
 * \f$ M_1 = \begin{bmatrix} A & a \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} b^{\textsf{T}} \\ D \end{bmatrix} \f$ such that \f$ B = a b^{\textsf{T}} \f$ holds.
 *
 * This function computes \f$ M_1 \f$, while \f$ M_2 \f$ can be computed by \ref CMRtwoSumDecomposeSecond.
 */

CMR_EXPORT
CMR_ERROR CMRtwoSumDecomposeFirst(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,         /**< Input matrix \f$ M \f$. */
  CMR_SEPA* sepa,             /**< 2-separation to decompose at. */
  CMR_CHRMAT** pfirst,        /**< Pointer for storing the first matrix \f$ M_1 \f$. */
  size_t* firstRowsOrigin,    /**< Array for storing the mapping from rows of \f$ M_1 \f$ to rows of \f$ M \f$ or to
                               **  \c SIZE_MAX; may be \c NULL. */
  size_t* firstColumnsOrigin, /**< Array for storing the mapping from columns of \f$ M_1 \f$ to columns of \f$ M \f$ or
                               **  to \c SIZE_MAX; may be \c NULL. */
  size_t* rowsToFirst,        /**< Array for storing the mapping from rows of \f$ M \f$ to rows of \f$ M_1 \f$ or to
                               **  \c SIZE_MAX; may be \c NULL. */
  size_t* columnsToFirst,     /**< Array for storing the mapping from columns of \f$ M \f$ to columns of \f$ M_1 \f$ or
                               **  to \c SIZE_MAX; may be \c NULL. */
  CMR_ELEMENT* pfirstMarker   /**< Pointer for storing the row index of \f$ c^{\textsf{T}} \f$ in \f$ M_1 \f$ or the
                               **  column index of \f$ a \f$ in \f$ M_1 \f$; may be \c NULL. */
);

/**
 * \brief Decomposes \p matrix as a 2-sum according to the 2-separation \p sepa and computing the second component.
 *
 * The input \p matrix \f$ M \f$ must have a 2-separation that is given by \p sepa, i.e., it can be reordered to look
 * like \f$ M = \begin{bmatrix} A & B \\ C & D \end{bmatrix} \f$, where \f$ \rank(B) + \rank(C) = 1 \f$.
 * If \f$ \rank(B) = \mathbb{O} \f$ then the two parts of the 2-sum are matrices
 * \f$ M_1 = \begin{bmatrix} A \\ c^{\textsf{T}} \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} d & D \end{bmatrix} \f$ such that \f$ C = d c^{\textsf{T}} \f$ holds.
 * Otherwise, the two parts of the 2-sum are matrices
 * \f$ M_1 = \begin{bmatrix} A & a \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} b^{\textsf{T}} \\ D \end{bmatrix} \f$ such that \f$ B = a b^{\textsf{T}} \f$ holds.
 *
 * This function computes \f$ M_2 \f$, while \f$ M_1 \f$ can be computed by \ref CMRtwoSumDecomposeFirst.
 */

CMR_EXPORT
CMR_ERROR CMRtwoSumDecomposeSecond(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,           /**< Input matrix \f$ M \f$. */
  CMR_SEPA* sepa,               /**< 2-separation to decompose at. */
  CMR_CHRMAT** psecond,         /**< Pointer for storing the second matrix \f$ M_2 \f$. */
  size_t* secondRowsOrigin,     /**< Array for storing the mapping from rows of \f$ M_2 \f$ to rows of \f$ M \f$ or to
                                 **  \c SIZE_MAX; may be \c NULL. */
  size_t* secondColumnsOrigin,  /**< Array for storing the mapping from columns of \f$ M_2 \f$ to columns of \f$ M \f$ or
                                 **  to \c SIZE_MAX; may be \c NULL. */
  size_t* rowsToSecond,         /**< Array for storing the mapping from rows of \f$ M \f$ to rows of \f$ M_2 \f$ or to
                                 **  \c SIZE_MAX; may be \c NULL. */
  size_t* columnsToSecond,      /**< Array for storing the mapping from columns of \f$ M \f$ to columns of \f$ M_2 \f$ or
                                 **  to \c SIZE_MAX; may be \c NULL. */
  CMR_ELEMENT* psecondMarker    /**< Pointer for storing the row index of \f$ b^{\textsf{T}} \f$ in \f$ M_2 \f$ or the
                                 **  column index of \f$ d \f$ in \f$ M_2 \f$; may be \c NULL. */
);

/**
 * \brief Constructs the 3-sum of the two matrices \p first and \p second via \p firstMarker1, \p firstMarker2,
 *        \p secondMarker1 and \p secondMarker2.
 *
 * Let \f$ A \f$ and \f$ B \f$ denote the matrices given by \p first and \p second, let \f$ A' \f$ be the matrix
 * \f$ A \f$ without the rows or columns indexed by \p firstMarker1 and \p firstMarker2, and let \f$ B' \f$ be the
 * matrix \f$ B \f$ without the rows or columns indexed by \p secondMarker1 and \p secondMarker2.
 * If \p firstMarker1 and \p firstMarker2 both index row vectors \f$ a_1^{\textsf{T}}, a_2^{\textsf{T}} \f$ of \f$ A \f$
 * then \p secondMarker1 and \p secondMarker2 must index column vectors \f$ b_1, b_2 \f$ of
 * \f$ B \f$. In this case the 3-sum is the matrix
 * \f[
 *   C := \begin{bmatrix}
 *     A' & \mathbb{O} \\
 *     b_1 a_1^{\textsf{T}} + b_2 a_2^{\textsf{T}} & B'
 *   \end{bmatrix}.
 * \f]
 * Otherwise, if \p firstMarker1 and \p firstMarker2 both index column vectors \f$ a_1, a_2 \f$ of \f$ A \f$ then
 * \p secondMarker1 and \p secondMarker2 must index row vectors \f$ b_1^{\textsf{T}}, b_2^{\textsf{T}} \f$ of \f$ B \f$.
 * In this case the 3-sum is the matrix
 * \f[
 *   C := \begin{bmatrix}
 *     A'         & a_1 b_1^{\textsf{T}} + a_2 b_2^{\textsf{T}}  \\
 *     \mathbb{O} & B'
 *   \end{bmatrix}.
 * \f]
 * Otherwise, if \p firstMarker1 indexes a row vector \f$ a_1^{\textsf{T}} \f$ and \p firstMarker2 indexes a column
 * vector \f$ a_2 \f$ of \f$ A \f$ then \p secondMarker1 must index a column vector \f$ b_1 \f$ of \f$ B \f$ and
 * \p secondMarker2 must index a row vector \f$ b_2^{\textsf{T}} \f$ of \f$ B \f$.
 * In this case the 3-sum is the matrix
 * \f[
 *   C := \begin{bmatrix}
 *     A'                   & a_2 b_2^{\textsf{T}}  \\
 *     b_1 a_1^{\textsf{T}} & B'
 *   \end{bmatrix}.
 * \f]
 * The remaining case is identical to the previous one, except that \p firstMarker1 and \p firstMarker2 as well as
 * \p secondMarker1 and \p secondMarker2 change roles.
 * The calculations are done modulo \p characteristic, where the value \f$ 3 \f$ yields numbers from \f$ \{-1,0,+1\} \f$.
 *
 * The resulting matrix \f$ C \f$ is created and stored in \p *presult.
 */

CMR_EXPORT
CMR_ERROR CMRthreeSum(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* first,          /**< First matrix. */
  CMR_CHRMAT* second,         /**< Second matrix. */
  CMR_ELEMENT firstMarker1,   /**< First marker element of first matrix. */
  CMR_ELEMENT secondMarker1,  /**< Second marker element of first matrix. */
  CMR_ELEMENT firstMarker2,   /**< First marker element of second matrix. */
  CMR_ELEMENT secondMarker2,  /**< Second marker element of second matrix. */
  int8_t characteristic,      /**< Field characteristic. */
  CMR_CHRMAT** presult        /**< Pointer for storing the result. */
);


#ifdef __cplusplus
}
#endif

#endif /* CMR_SEPARATION_H */
