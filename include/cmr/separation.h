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
 * Depending on the input, one of two variants of the 2-sum of two matrices \f$ M_1 \f$, given by \p first and
 * \f$ M_2 \f$, given by \p second are computed. For the first variant, consider the matrices
 * \f$ M_1 = \begin{bmatrix} A \\ c^{\textsf{T}} \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} d & D \end{bmatrix} \f$. Then the 2-sum is the matrix
 * \f[
 *   M := \begin{bmatrix}
 *     A & \mathbb{O} \\
 *     d c^{\textsf{T}} & D
 *   \end{bmatrix}.
 * \f]
 * To obtain this 2-sum, the arrays \p firstSpecialRows and \p secondSpecialColumns must be arrays of length 1,
 * consisting of the row index of \f$ c^{\textsf{T}} \f$ of \f$ M_1 \f$ and of the column index of \f$ d \f$ of
 * \f$ M_2 \f$.
 * For the other variant, consider the matrices \f$ M_1 = \begin{bmatrix} A & a \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} b^{\textsf{T}} \\ D \end{bmatrix} \f$. Then the 2-sum is the matrix
 * \f[
 *   M := \begin{bmatrix}
 *     A & a b^{\textsf{T}} \\
 *     \mathbb{O} & D
 *   \end{bmatrix}.
 * \f]
 * To obtain this 2-sum, the arrays \p firstSpecialColumns and \p secondSpecialRows must be arrays of length 1,
 * consisting of the column index of \f$ a \f$ of \f$ M_1 \f$ and of the row index of \f$ b^{\textsf{T}} \f$ of
 * \f$ M_2 \f$.
 *
 * The calculations are done modulo \p characteristic, where the value \f$ 3 \f$ yields numbers from
 * \f$ \{-1,0,+1\} \f$.
 *
 * The resulting matrix \f$ M \f$ is created and stored in \p *presult.
 */

CMR_EXPORT
CMR_ERROR CMRtwoSumCompose(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* first,            /**< First matrix \f$ M_1 \f$. */
  CMR_CHRMAT* second,           /**< Second matrix \f$ M_2 \f$. */
  size_t* firstSpecialRows,     /**< Array of length 1 with row index of \f$ c^{\textsf{T}} \f$ or \c NULL. */
  size_t* firstSpecialColumns,  /**< Array of length 1 with column index of \f$ a \f$ or \c NULL. */
  size_t* secondSpecialRows,    /**< Array of length 1 with row index of \f$ b^{\textsf{T}} \f$ or \c NULL. */
  size_t* secondSpecialColumns, /**< Array of length 1 with column index of \f$ d \f$ or \c NULL. */
  int8_t characteristic,        /**< Field characteristic. */
  CMR_CHRMAT** presult          /**< Pointer for storing the result. */
);

/**
 * \brief Decomposes \p matrix as a 2-sum according to the 2-separation \p sepa, computing the first component.
 *
 * The input \p matrix \f$ M \f$ must have a 2-separation that is given by \p sepa, i.e., it can be reordered to look
 * like \f$ M = \begin{bmatrix} A & B \\ C & D \end{bmatrix} \f$, where \f$ \text{rank}(B) + \text{rank}(C) = 1 \f$.
 * If \f$ \text{rank}(B) = \mathbb{O} \f$ then the two components of the 2-sum are matrices
 * \f$ M_1 = \begin{bmatrix} A \\ c^{\textsf{T}} \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} d & D \end{bmatrix} \f$ such that \f$ C = d c^{\textsf{T}} \f$ holds and such that
 * \f$ c^{\textsf{T}} \f$ is an actual row of \f$ M \f$.
 * Otherwise, the two components of the 2-sum are matrices
 * \f$ M_1 = \begin{bmatrix} A & a \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} b^{\textsf{T}} \\ D \end{bmatrix} \f$ such that \f$ B = a b^{\textsf{T}} \f$ holds and such
 * that \f$ a \f$ is an actual column of \f$ M \f$.
 *
 * This function computes \f$ M_1 \f$, while \f$ M_2 \f$ can be computed by \ref CMRtwoSumDecomposeSecond.
 *
 * \note For the first variant, if \p firstSpecialRows is not \c NULL then \p *firstSpecialRows[0] will be the last row.
 *
 * \note For the second variant, if \p firstSpecialColumns is not \c NULL then \p *firstSpecialColumns[0] will be the
 *       last column.
 */

CMR_EXPORT
CMR_ERROR CMRtwoSumDecomposeFirst(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,         /**< Input matrix \f$ M \f$. */
  CMR_SEPA* sepa,             /**< 2-separation to decompose at. */
  CMR_CHRMAT** pfirst,        /**< Pointer for storing the first matrix \f$ M_1 \f$. */
  size_t* firstRowsOrigin,    /**< Array for storing the mapping from rows of \f$ M_1 \f$ to rows of \f$ M \f$;
                               **  also set for the extra row if applicable; may be \c NULL. */
  size_t* firstColumnsOrigin, /**< Array for storing the mapping from columns of \f$ M_1 \f$ to columns of \f$ M \f$;
                               **  also set for the extra column if applicable; may be \c NULL. */
  size_t* rowsToFirst,        /**< Array for storing the mapping from rows of \f$ M \f$ to rows of \f$ M_1 \f$ or to
                               **  \c SIZE_MAX; may be \c NULL. */
  size_t* columnsToFirst,     /**< Array for storing the mapping from columns of \f$ M \f$ to columns of \f$ M_1 \f$ or
                               **  to \c SIZE_MAX; may be \c NULL. */
  size_t* firstSpecialRows,   /**< Either \c NULL or an array for storing the row index of \f$ c^{\textsf{T}} \f$ in
                               **  \f$ M_1 \f$. */
  size_t* firstSpecialColumns /**< Either \c NULL or an array for storing the column index of \f$ a \f$ in
                               **  \f$ M_1 \f$. */
);

/**
 * \brief Decomposes \p matrix as a 2-sum according to the 2-separation \p sepa, computing the second component.
 *
 * The input \p matrix \f$ M \f$ must have a 2-separation that is given by \p sepa, i.e., it can be reordered to look
 * like \f$ M = \begin{bmatrix} A & B \\ C & D \end{bmatrix} \f$, where \f$ \text{rank}(B) + \text{rank}(C) = 1 \f$.
 * If \f$ \text{rank}(B) = \mathbb{O} \f$ then the two components of the 2-sum are matrices
 * \f$ M_1 = \begin{bmatrix} A \\ c^{\textsf{T}} \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} d & D \end{bmatrix} \f$ such that \f$ C = d c^{\textsf{T}} \f$ holds and such that
 * \f$ c^{\textsf{T}} \f$ is an actual row of \f$ M \f$.
 * Otherwise, the two components of the 2-sum are matrices
 * \f$ M_1 = \begin{bmatrix} A & a \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} b^{\textsf{T}} \\ D \end{bmatrix} \f$ such that \f$ B = a b^{\textsf{T}} \f$ holds and such
 * that \f$ a \f$ is an actual column of \f$ M \f$.
 *
 * This function computes \f$ M_2 \f$, while \f$ M_1 \f$ can be computed by \ref CMRtwoSumDecomposeFirst.
 *
 * \note For the first variant, if \p secondSpecialColumns is not \c NULL then \p *secondSpecialColumns[0] will be the
 *       first column.
 *
 * \note For the second variant, if \p secondSpecialRows is not \c NULL then \p *secondSpecialRows[0] will be the first
 *       row.
 */

CMR_EXPORT
CMR_ERROR CMRtwoSumDecomposeSecond(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,           /**< Input matrix \f$ M \f$. */
  CMR_SEPA* sepa,               /**< 2-separation to decompose at. */
  CMR_CHRMAT** psecond,         /**< Pointer for storing the second matrix \f$ M_2 \f$. */
  size_t* secondRowsOrigin,     /**< Array for storing the mapping from rows of \f$ M_2 \f$ to rows of \f$ M \f$;
                                 **  also set for the extra row if applicable; may be \c NULL. */
  size_t* secondColumnsOrigin,  /**< Array for storing the mapping from columns of \f$ M_2 \f$ to columns of \f$ M \f$;
                                 **  also set for the extra column if applicable; may be \c NULL. */
  size_t* rowsToSecond,         /**< Array for storing the mapping from rows of \f$ M \f$ to rows of \f$ M_2 \f$ or to
                                 **  \c SIZE_MAX; may be \c NULL. */
  size_t* columnsToSecond,      /**< Array for storing the mapping from columns of \f$ M \f$ to columns of \f$ M_2 \f$ or
                                 **  to \c SIZE_MAX; may be \c NULL. */
  size_t* secondSpecialRows,    /**< Either \c NULL or an array for storing the row index of \f$ b^{\textsf{T}} \f$ in
                                 **  \f$ M_2 \f$. */
  size_t* secondSpecialColumns  /**< Either \c NULL or an array for storing the column index of \f$ d \f$ in
                                 **  \f$ M_2 \f$. */
);

/**
 * \brief Constructs the Seymour 3-sum of the two matrices \p first and \p second.
 *
 * Let \f$ M_1 \f$ and \f$ M_2 \f$ denote the matrices given by \p first and \p second, let \f$ A \f$ be the matrix
 * \f$ M_1 \f$ without the row indexed by \p firstSpecialRows[0] and the columns indexed by \p firstSpecialColumns[0]
 * and \p firstSpecialColumns[1]. After reordering these to be last, \f$ M_1 \f$ must be of the form
 * \f$
 *   M_1 = \begin{bmatrix}
 *     A & a & a \\
 *     c^{\textsf{T}} & 0 & \varepsilon
 *   \end{bmatrix},
 * \f$
 * where \f$ \varepsilon \in \{-1,+1 \} \f$ (otherwise, \c CMR_ERROR_STRUCTURE is returned).
 * Similarly, let \f$ D \f$ be the matrix \f$ M_2 \f$ without the row indexed by \p secondSpecialRows[0] and the columns
 * indexed by \p secondSpecialColumns[0] and \p secondSpecialColumns[1]. After reordering these to be first, \f$ M_2 \f$
 * must be of the form
 * \f$
 *   M_2 = \begin{bmatrix}
 *     \varepsilon & 0 & b^{\textsf{T}} \\
 *     d & d & D
 *   \end{bmatrix}
 * \f$
 * with the same \f$ \varepsilon \f$ (otherwise, \c CMR_ERROR_STRUCTURE is returned).
 * The 3-sum of \f$ M_1 \f$ and \f$ M_2 \f$ (at these special rows/columns) is the matrix
 * \f[
 *   M = \begin{bmatrix}
 *     A & a b^{\textsf{T}} \\
 *     d c^{\textsf{T}} & D
 *   \end{bmatrix}.
 * \f]
 * The calculations are done modulo \p characteristic, where the value \f$ 3 \f$ yields numbers from
 * \f$ \{-1,0,+1\} \f$.
 *
 * The resulting matrix \f$ M \f$ is created and stored in \p *presult.
 */

CMR_EXPORT
CMR_ERROR CMRthreeSumSeymourCompose(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* first,            /**< First matrix. */
  CMR_CHRMAT* second,           /**< Second matrix. */
  size_t* firstSpecialRows,     /**< Array of length 1 with row index of
                                 **  \f$ \begin{bmatrix} c^{\textsf{T}} & 0 & \varepsilon \end{bmatrix} \f$
                                 ** in \f$ M_1 \f$. */
  size_t* firstSpecialColumns,  /**< Array of length 2 with column indices of
                                 **  \f$ \begin{bmatrix} a \\ 0 \end{bmatrix} \f$ and
                                 **  \f$ \begin{bmatrix} a \\ \varepsilon \end{bmatrix} \f$
                                 ** in \f$ M_1 \f$. */
  size_t* secondSpecialRows,    /**< Array of length 1 with row index of \f$ b^{\textsf{T}} \f$
                                 **  \f$ \begin{bmatrix} \varepsilon & 0 & b^{\textsf{T}} \end{bmatrix} \f$
                                 ** in \f$ M_2 \f$. */
  size_t* secondSpecialColumns, /**< Array of length 2 with column indices
                                 **  \f$ \begin{bmatrix} \varepsilon \\ d \end{bmatrix} \f$ and
                                 **  \f$ \begin{bmatrix} 0 \\ d \end{bmatrix} \f$
                                 **  in \f$ M_2 \f$. */
  int8_t characteristic,      /**< Field characteristic. */
  CMR_CHRMAT** presult        /**< Pointer for storing the result. */
);

/**
 * \brief Decomposes \p matrix as a Seymour 3-sum according to the 3-separation \p sepa, computing \f$ \varepsilon \f$.
 *
 * The input \p matrix \f$ M \f$ must have a 3-separation that is given by \p sepa, i.e., it can be reordered to look
 * like \f$ M = \begin{bmatrix} A & B \\ C & D \end{bmatrix} \f$, where \f$ \text{rank}(B) = \text{rank}(C) = 1 \f$.
 * The two components of the 3-sum are matrices
 * \f$ M_1 = \begin{bmatrix} A & a & a \\ c^{\textsf{T}} & 0 & \varepsilon \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} \varepsilon & 0 & b^{\textsf{T}} \\ d & d & D \end{bmatrix} \f$ such that
 * \f$ B = a b^{\textsf{T}} \f$ and \f$ C = d c^{\textsf{T}} \f$ hold and such that
 * \f$ a \f$ and \f$ c^{\textsf{T}} \f$ are an actual column and row of \f$ M \f$.
 * Consequently, \f$ b^{\textsf{T}} \f$ and \f$ d \f$ are (possibly negated) rows and columns of \f$ M \f$.
 *
 * The value of \f$ \varepsilon \in \{-1,+1\} \f$ must be so that there exists a singular submatrix of \f$ M_1 \f$ with
 * exactly two nonzeros per row and per column that covers the top-left \f$ \varepsilon \f$-entry.
 *
 * This function only computes \f$ \varepsilon \f$; the matrices \f$ M_1 \f$ and \f$ M_2 \f$ can be computed by
 * \ref CMRthreeSumSeymourDecomposeFirst and \ref CMRthreeSumSeymourDecomposeSecond, respectively.
 */

CMR_EXPORT
CMR_ERROR CMRthreeSumSeymourDecomposeEpsilon(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Input matrix \f$ M \f$. */
  CMR_CHRMAT* transpose,  /**< Transpose matrix \f$ M^{\textsf{T}} \f$. */
  CMR_SEPA* sepa,         /**< 3-separation to decompose at. */
  char* pepsilon          /**< Pointer for storing a correct value of \f$ \varepsilon \f$. */
);

/**
 * \brief Decomposes \p matrix as a Seymour 3-sum according to the 3-separation \p sepa, computing the first component.
 *
 * The input \p matrix \f$ M \f$ must have a 3-separation that is given by \p sepa, i.e., it can be reordered to look
 * like \f$ M = \begin{bmatrix} A & B \\ C & D \end{bmatrix} \f$, where \f$ \text{rank}(B) = \text{rank}(C) = 1 \f$.
 * The two components of the 3-sum are matrices
 * \f$ M_1 = \begin{bmatrix} A & a & a \\ c^{\textsf{T}} & 0 & \varepsilon \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} \varepsilon & 0 & b^{\textsf{T}} \\ d & d & D \end{bmatrix} \f$ such that
 * \f$ B = a b^{\textsf{T}} \f$ and \f$ C = d c^{\textsf{T}} \f$ hold and such that
 * \f$ a \f$ and \f$ c^{\textsf{T}} \f$ are an actual column and row of \f$ M \f$.
 * Consequently, \f$ b^{\textsf{T}} \f$ and \f$ d \f$ are (possibly negated) rows and columns of \f$ M \f$.
 *
 * The value of \f$ \varepsilon \in \{-1,+1\} \f$ must be given by \p epsilon and should be computed by
 * \ref CMRthreeSumSeymourDecomposeEpsilon.
 *
 * This function computes \f$ M_1 \f$, while \f$ M_2 \f$ can be computed by \ref CMRthreeSumSeymourDecomposeSecond.
 *
 * If \p firstSpecialRows is not \c NULL, then \p firstSpecialRows[0] will refer to the last row of \f$ M_1 \f$.
 * If \p firstSpecialColumns is not \c NULL, then \p firstSpecialColumns[0] will refer to the second-to last column and
 * \p firstSpecialColumns[1] will refer to the last column of \f$ M_2 \f$.
 */

CMR_EXPORT
CMR_ERROR CMRthreeSumSeymourDecomposeFirst(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,         /**< Input matrix \f$ M \f$. */
  CMR_SEPA* sepa,             /**< 3-separation to decompose at. */
  char epsilon,               /**< Value of \f$ \varepsilon \f$. */
  CMR_CHRMAT** pfirst,        /**< Pointer for storing the first matrix \f$ M_1 \f$. */
  size_t* firstRowsOrigin,    /**< Array for storing the mapping from rows of \f$ M_1 \f$ to rows of \f$ M \f$;
                               **  also set for the extra row if applicable, even if negated; may be \c NULL. */
  size_t* firstColumnsOrigin, /**< Array for storing the mapping from columns of \f$ M_1 \f$ to columns of \f$ M \f$;
                               **  also set for the extra column if applicable, even if negated; may be \c NULL. */
  size_t* rowsToFirst,        /**< Array for storing the mapping from rows of \f$ M \f$ to rows of \f$ M_1 \f$ or to
                               **  \c SIZE_MAX; may be \c NULL. */
  size_t* columnsToFirst,     /**< Array for storing the mapping from columns of \f$ M \f$ to columns of \f$ M_1 \f$
                               **  or to \c SIZE_MAX; may be \c NULL. */
  size_t* firstSpecialRows,   /**< Array of length 1 for storing the row index of
                               **  \f$ \begin{bmatrix} c^{\textsf{T}} & 0 & \varepsilon \end{bmatrix} \f$
                               ** in \f$ M_1 \f$; may be \c NULL. */
  size_t* firstSpecialColumns /**< Array of length 2 for storing the column indices of
                               **  \f$ \begin{bmatrix} a \\ 0 \end{bmatrix} \f$ and
                               **  \f$ \begin{bmatrix} a \\ \varepsilon \end{bmatrix} \f$
                               ** in \f$ M_1 \f$; may be \c NULL. */
);

/**
 * \brief Decomposes \p matrix as a Seymour 3-sum according to the 3-separation \p sepa, computing the second component.
 *
 * The input \p matrix \f$ M \f$ must have a 3-separation that is given by \p sepa, i.e., it can be reordered to look
 * like \f$ M = \begin{bmatrix} A & B \\ C & D \end{bmatrix} \f$, where \f$ \text{rank}(B) = \text{rank}(C) = 1 \f$.
 * The two components of the 3-sum are matrices
 * \f$ M_1 = \begin{bmatrix} A & a & a \\ c^{\textsf{T}} & 0 & \varepsilon \end{bmatrix} \f$ and
 * \f$ M_2 = \begin{bmatrix} \varepsilon & 0 & b^{\textsf{T}} \\ d & d & D \end{bmatrix} \f$ such that
 * \f$ B = a b^{\textsf{T}} \f$ and \f$ C = d c^{\textsf{T}} \f$ hold and such that
 * \f$ a \f$ and \f$ c^{\textsf{T}} \f$ are an actual column and row of \f$ M \f$.
 * Consequently, \f$ b^{\textsf{T}} \f$ and \f$ d \f$ are (possibly negated) rows and columns of \f$ M \f$.
 *
 * The value of \f$ \varepsilon \in \{-1,+1\} \f$ must be given by \p epsilon and should be computed by
 * \ref CMRthreeSumSeymourDecomposeEpsilon.
 *
 * This function computes \f$ M_2 \f$, while \f$ M_1 \f$ can be computed by \ref CMRthreeSumSeymourDecomposeFirst.
 *
 * \note If \p secondSpecialRows is not \c NULL, then \p secondSpecialRows[0] will refer to the first row of
 *       \f$ M_2 \f$.
 * \note If \p secondSpecialColumns is not \c NULL, then \p secondSpecialColumns[0] will refer to the first column and
 *       \p secondSpecialColumns[1] will refer to the second column of \f$ M_2 \f$.
 */

CMR_EXPORT
CMR_ERROR CMRthreeSumSeymourDecomposeSecond(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,           /**< Input matrix \f$ M \f$. */
  CMR_SEPA* sepa,               /**< 3-separation to decompose at. */
  char epsilon,                 /**< Value of \f$ \varepsilon \f$. */
  CMR_CHRMAT** psecond,         /**< Pointer for storing the second matrix \f$ M_2 \f$. */
  size_t* secondRowsOrigin,     /**< Array for storing the mapping from rows of \f$ M_2 \f$ to rows of \f$ M \f$;
                                 **  also set for the extra row if applicable, even if negated; may be \c NULL. */
  size_t* secondColumnsOrigin,  /**< Array for storing the mapping from columns of \f$ M_2 \f$ to columns of \f$ M \f$;
                                 **  also set for the extra column if applicable, even if negated; may be \c NULL. */
  size_t* rowsToSecond,         /**< Array for storing the mapping from rows of \f$ M \f$ to rows of \f$ M_2 \f$ or to
                                 **  \c SIZE_MAX; may be \c NULL. */
  size_t* columnsToSecond,      /**< Array for storing the mapping from columns of \f$ M \f$ to columns of \f$ M_2 \f$
                                 **  or to \c SIZE_MAX; may be \c NULL. */
  size_t* secondSpecialRows,    /**< Array of length 1 for storing the row index of
                                 **  \f$ \begin{bmatrix} \varepsilon & 0 & b^{\textsf{T}} \end{bmatrix} \f$
                                 ** in \f$ M_2 \f$; may be \c NULL. */
  size_t* secondSpecialColumns  /**< Array of length 2 for storing the column indices of
                                 **  \f$ \begin{bmatrix} \varepsilon & d \end{bmatrix} \f$ and
                                 **  \f$ \begin{bmatrix} 0 & d \end{bmatrix} \f$
                                 ** in \f$ M_2 \f$; may be \c NULL. */
);

/**
 * \brief Constructs the Truemper 3-sum of the two matrices \p first and \p second at connecting rows
 *        \p firstSpecialRows and \p secondSpecialRows and columns \p firstSpecialColumns and \p secondSpecialColumns.
 *
 * Let \f$ M_1 \f$ and \f$ M_2 \f$ denote the matrices given by \p first and \p second, let \f$ A \f$ be the matrix
 * \f$ M_1 \f$ without the rows \p firstSpecialRows[0] and \p firstSpecialRows[1] and column \p firstSpecialColumns[2].
 * After permuting these to be last, \f$ M_1 \f$ must be of the form
 * \f[
 *   M_1 = \begin{bmatrix}
 *     A & \mathbb{O} \\
 *     C_{i,\star} & \alpha \\
 *     C_{j,\star} & \beta
 *   \end{bmatrix},
 * \f]
 * where \f$ \alpha,\beta \in \{-1,+1 \} \f$ (otherwise, \c CMR_ERROR_STRUCTURE is returned). Let \f$ D \f$ be the
 * matrix \f$ M_2 \f$ without the row \p secondSpecialRows[0] and columns \p secondSpecialColumns[0] and
 * \p secondSpecialColumns[1]. After reordering these to be first, \f$ M_2 \f$ must be of the form
 * \f[
 *   M_2 = \begin{bmatrix}
 *     \gamma & \delta & \mathbb{O}^{\textsf{T}} \\
 *     C_{\star,k} & C_{\star,\ell} & D
 *   \end{bmatrix},
 * \f]
 * where \f$ \gamma,\delta \in \{ -1,+1 \} \f$ (otherwise, \c CMR_ERROR_STRUCTURE is returned) and such that the matrix
 * \f[
 *   N = \begin{bmatrix}
 *     \gamma & \delta & 0 \\
 *     C_{i,k} & C_{i,\ell} & \alpha \\
 *     C_{j,k} & C_{j,\ell} & \beta
 *   \end{bmatrix}
 * \f]
 * is totally unimodular (otherwise, \c CMR_ERROR_STRUCTURE is returned). The columns \p firstSpecialColumns[0] and
 * \p firstSpecialColumns[1] indicate the columns of \f$ M_1 \f$ that shall correspond to \f$ C_{\star,k}\f$ and
 * \f$ C_{\star,\ell} \f$, respectively. Similarly, the rows \p secondSpecialRows[1] and \p secondSpecialRows[2]
 * indicate the rows of \f$ M_2 \f$ that shall correspond to \f$ C_{i,\star}\f$ and \f$ C_{j,\star} \f$, respectively.
 *
 * \note The 2-by-2 submatrix of \f$ M_1 \f$ indexed by rows \p firstSpecialRows[0] and \p firstSpecialRows[1] and
 *       columns \p firstSpecialColumns[0] and \p firstSpecialColumns[1] must be identical to the submatrix of
 *       \f$ M_2 \f$ indexed by rows \p secondSpecialRows[1] and \p secondSpecialRows[2] and columns
 *       \p secondSpecialColumns[0] and \p secondSpecialColumns[1], which is the matrix \f$ C_{\{i,j\},\{k,\ell\}} \f$.
 *
 * The 3-sum of \f$ M_1 \f$ and \f$ M_2 \f$ (at these rows/columns) is the matrix
 * \f[
 *   M = \begin{bmatrix}
 *     A & \mathbb{O} \\
 *     C & D
 *   \end{bmatrix},
 * \f]
 * where \f$ C \f$ is the unique rank-2 matrix having linearly independent rows \f$ C_{i,\star} \f$ and
 * \f$ C_{j,\star} \f$ and linearly independent columns \f$ C_{\star,k} \f$ and \f$ C_{\star,\ell} \f$.
 * The calculations are done modulo \p characteristic, where the value \f$ 3 \f$ yields numbers from
 * \f$ \{-1,0,+1\} \f$.
 *
 * The resulting matrix \f$ M \f$ is created and stored in \p *presult.
 */

CMR_EXPORT
CMR_ERROR CMRthreeSumTruemperCompose(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* first,            /**< First matrix. */
  CMR_CHRMAT* second,           /**< Second matrix. */
  size_t* firstSpecialRows,     /**< Array of length 2 with the last two rows of the connecting submatrix in \p first. */
  size_t* firstSpecialColumns,  /**< Array of length 3 with all columns of the connecting submatrix in \p first. */
  size_t* secondSpecialRows,    /**< Array of length 3 with the all rows of the connecting submatrix in \p second. */
  size_t* secondSpecialColumns, /**< Array of length 2 with the first two columns of connecting submatrix in \p second. */
  int8_t characteristic,        /**< Field characteristic. */
  CMR_CHRMAT** presult          /**< Pointer for storing the result. */
);


/**
 * \brief Decomposes \p matrix as a Truemper 3-sum according to the 3-separation \p sepa, computing the connecing
 *        matrix.
 *
 * The input \p matrix \f$ M \f$ must have a 3-separation that is given by \p sepa, i.e., it can be reordered to look
 * like \f$ M = \begin{bmatrix} A & \mathbb{O} \\ C & D \end{bmatrix} \f$, where \f$ \text{rank}(C) = 2 \f$.
 * The two components of the 3-sum are matrices
 * \f[
 *   M_1 = \begin{bmatrix}
 *     A & \mathbb{O} \\
 *     C_{i,\star} & 1 \\
 *     C_{j,\star} & \beta
 *   \end{bmatrix}
 * \f]
 * and
 * \f[
 *   M_2 = \begin{bmatrix}
 *     \gamma & 1 & \mathbb{O}^{\textsf{T}} \\
 *     C_{\star,k} & C_{\star,\ell} & D
 *   \end{bmatrix},
 * \f]
 * where \f$ \beta,\gamma \in \{-1,+1 \} \f$, \f$\text{rank}(C_{\{i,j\},\{k,\ell\}}) = 2\f$ and such that
 * \f[
 *   N := \begin{bmatrix}
 *     \gamma & 1 & 0 \\
 *     C_{i,k} & C_{i,\ell} & 1 \\
 *     C_{j,k} & C_{j,\ell} & \beta
 *   \end{bmatrix}
 * \f]
 * is totally unimodular.
 *
 * The value of \f$ \beta \in \{-1,+1\} \f$ must be so that there exists a singular submatrix of \f$ M_1 \f$ with
 * exactly two nonzeros per row and per column that covers the bottom-right \f$ \beta \f$-entry.
 *
 * This function only computes the indices \f$ i,j,k,\ell \f$ as well as values for \f$ \beta \f$ and \f$ \gamma \f$;
 * the matrices \f$ M_1 \f$ and \f$ M_2 \f$ can be computed by \ref CMRthreeSumTruemperDecomposeFirst and
 * \ref CMRthreeSumTruemperDecomposeSecond, respectively.
 */

CMR_EXPORT
CMR_ERROR CMRthreeSumTruemperDecomposeConnecting(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,     /**< Input matrix \f$ M \f$. */
  CMR_CHRMAT* transpose,  /**< Transpose matrix \f$ M^{\textsf{T}} \f$. */
  CMR_SEPA* sepa,         /**< 3-separation to decompose at. */
  size_t* specialRows,    /**< Array of length 2 for storing the rows \f$ i \f$ and \f$ j \f$ as rows of \f$ M \f$. */
  size_t* specialColumns, /**< Array of length 2 for storing the columns \f$ k \f$ and \f$ \ell \f$ as rows of
                           **  \f$ M \f$. */
  char* pgamma,           /**< Pointer for storing a correct value of \f$ \gamma \f$; may be \c NULL. */
  char* pbeta             /**< Pointer for storing a correct value of \f$ \beta \f$; may be \c NULL. */
);


/**
 * \brief Decomposes \p matrix as a Truemper 3-sum according to the 3-separation \p sepa, computing the first component.
 *
 * The input \p matrix \f$ M \f$ must have a 3-separation that is given by \p sepa, i.e., it can be reordered to look
 * like \f$ M = \begin{bmatrix} A & \mathbb{O} \\ C & D \end{bmatrix} \f$, where \f$ \text{rank}(C) = 2 \f$.
 * The two components of the 3-sum are matrices
 * \f[
 *   M_1 = \begin{bmatrix}
 *     A & \mathbb{O} \\
 *     C_{i,\star} & 1 \\
 *     C_{j,\star} & \beta
 *   \end{bmatrix}
 * \f]
 * and
 * \f[
 *   M_2 = \begin{bmatrix}
 *     \gamma & 1 & \mathbb{O}^{\textsf{T}} \\
 *     C_{\star,k} & C_{\star,\ell} & D
 *   \end{bmatrix},
 * \f]
 * where \f$ \beta,\gamma \in \{-1,+1 \} \f$, \f$\text{rank}(C_{\{i,j\},\{k,\ell\}}) = 2\f$ and such that
 * \f[
 *   N := \begin{bmatrix}
 *     \gamma & 1 & 0 \\
 *     C_{i,k} & C_{i,\ell} & 1 \\
 *     C_{j,k} & C_{j,\ell} & \beta
 *   \end{bmatrix}
 * \f]
 * is totally unimodular.
 *
 * The value of \f$ \beta \in \{-1,+1\} \f$, given via \p beta. The row indices \f$ i,j \f$, given via \p specialRows,
 * and column indices \f$ k,\ell \f$, given via \p specialColumns, must index a rank-2 submatrix of \f$ C \f$. They
 * should be computed by \ref CMRthreeSumTruemperDecomposeConnecting.
 *
 * This function computes \f$ M_1 \f$, while \f$ M_2 \f$ can be computed by \ref CMRthreeSumTruemperDecomposeSecond.
 *
 * If \p firstSpecialRows is not \c NULL, then \p firstSpecialRows[0] and \p firstSpecialRows[1] will refer to the
 * last two rows of \f$ M_1 \f$.
 * If \p firstSpecialColumns is not \c NULL, then \p firstSpecialColumns[0] and \p firstSpecialColumns[1] will refer to
 * the two columns containing \f$ C_{\star,k}\f$ and \f$ C_{\star,\ell} \f$ as columns of \f$ M_1 \f$, and
 * \p firstSpecialColumns[2] will refer to the last (artificial) column.
 */

CMR_EXPORT
CMR_ERROR CMRthreeSumTruemperDecomposeFirst(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,         /**< Input matrix \f$ M \f$. */
  CMR_SEPA* sepa,             /**< 3-separation to decompose at. */
  size_t* specialRows,        /**< Array of length 2 with the rows \f$ i \f$ and \f$ j \f$ as rows of \f$ M \f$. */
  size_t* specialColumns,     /**< Array of length 2 with the columns \f$ k \f$ and \f$ \ell \f$ as rows of \f$ M \f$. */
  char beta,                  /**< Value of \f$ \beta \f$. */
  CMR_CHRMAT** pfirst,        /**< Pointer for storing the first matrix \f$ M_1 \f$. */
  size_t* firstRowsOrigin,    /**< Array for storing the mapping from rows of \f$ M_1 \f$ to rows of \f$ M \f$;
                               **  may be \c NULL. */
  size_t* firstColumnsOrigin, /**< Array for storing the mapping from columns of \f$ M_1 \f$ to columns of \f$ M \f$;
                               **  set to \c SIZE_MAX for the last column; may be \c NULL. */
  size_t* rowsToFirst,        /**< Array for storing the mapping from rows of \f$ M \f$ to rows of \f$ M_1 \f$ or to
                               **  \c SIZE_MAX; may be \c NULL. */
  size_t* columnsToFirst,     /**< Array for storing the mapping from columns of \f$ M \f$ to columns of \f$ M_1 \f$
                               **  or to \c SIZE_MAX; may be \c NULL. */
  size_t* firstSpecialRows,   /**< Array of length 2 for storing the row indices of \f$ C_{i,\star} \f$ and of
                               **  \f$ C_{j,\star} \f$ in \f$ M_1 \f$; may be \c NULL. */
  size_t* firstSpecialColumns /**< Array of length 3 for storing the column indices of \f$ C_{\star,k} \f$, of
                               **  \f$ C_{\star,\ell} \f$ and of the artificial column in \f$ M_1 \f$;
                               **  may be \c NULL. */
);

/**
 * \brief Decomposes \p matrix as a Truemper 3-sum according to the 3-separation \p sepa, computing the second
 *        component.
 *
 * The input \p matrix \f$ M \f$ must have a 3-separation that is given by \p sepa, i.e., it can be reordered to look
 * like \f$ M = \begin{bmatrix} A & \mathbb{O} \\ C & D \end{bmatrix} \f$, where \f$ \text{rank}(C) = 2 \f$.
 * The two components of the 3-sum are matrices
 * \f[
 *   M_1 = \begin{bmatrix}
 *     A & \mathbb{O} \\
 *     C_{i,\star} & 1 \\
 *     C_{j,\star} & \beta
 *   \end{bmatrix}
 * \f]
 * and
 * \f[
 *   M_2 = \begin{bmatrix}
 *     \gamma & 1 & \mathbb{O}^{\textsf{T}} \\
 *     C_{\star,k} & C_{\star,\ell} & D
 *   \end{bmatrix},
 * \f]
 * where \f$ \beta,\gamma \in \{-1,+1 \} \f$, \f$\text{rank}(C_{\{i,j\},\{k,\ell\}}) = 2\f$ and such that
 * \f[
 *   N := \begin{bmatrix}
 *     \gamma & 1 & 0 \\
 *     C_{i,k} & C_{i,\ell} & 1 \\
 *     C_{j,k} & C_{j,\ell} & \beta
 *   \end{bmatrix}
 * \f]
 * is totally unimodular.
 *
 * The value of \f$ \gamma \in \{-1,+1\} \f$, given via \p gamma. The row indices \f$ i,j \f$, given via \p specialRows,
 * and column indices \f$ k,\ell \f$, given via \p specialColumns, must index a rank-2 submatrix of \f$ C \f$. They
 * should be computed by \ref CMRthreeSumTruemperDecomposeConnecting.
 *
 * This function computes \f$ M_2 \f$, while \f$ M_1 \f$ can be computed by \ref CMRthreeSumTruemperDecomposeFirst.
 *
 * If \p secondSpecialRows is not \c NULL, then \p secondSpecialRows[0] will refer to the first (artificial) row and
 * \p secondSpecialRows[1] and \p secondSpecialRows[2] will refer to the two rows containing \f$ C_{i,\star} \f$ and
 * \f$ C_{j,\star} \f$ as rows of \f$ M_2 \f$.
 * If \p secondSpecialColumns is not \c NULL, then \p secondSpecialColumns[0] and \p secondSpecialColumns[1] will refer
 * to the first two columns of \f$ M_2 \f$.
 */

CMR_EXPORT
CMR_ERROR CMRthreeSumTruemperDecomposeSecond(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,           /**< Input matrix \f$ M \f$. */
  CMR_SEPA* sepa,               /**< 3-separation to decompose at. */
  size_t* specialRows,          /**< Array of length 2 with the rows \f$ i \f$ and \f$ j \f$ as rows of \f$ M \f$. */
  size_t* specialColumns,       /**< Array of length 2 with the columns \f$ k \f$ and \f$ \ell \f$ as rows of
                                 **  \f$ M \f$. */
  char gamma,                   /**< Value of \f$ \gamma \f$. */
  CMR_CHRMAT** psecond,         /**< Pointer for storing the second matrix \f$ M_2 \f$. */
  size_t* secondRowsOrigin,     /**< Array for storing the mapping from rows of \f$ M_2 \f$ to rows of \f$ M \f$;
                                 **  set to \c SIZE_MAX for the first row; may be \c NULL. */
  size_t* secondColumnsOrigin,  /**< Array for storing the mapping from columns of \f$ M_2 \f$ to columns of \f$ M \f$;
                                 **  may be \c NULL. */
  size_t* rowsToSecond,         /**< Array for storing the mapping from rows of \f$ M \f$ to rows of \f$ M_2 \f$ or to
                                 **  \c SIZE_MAX; may be \c NULL. */
  size_t* columnsToSecond,      /**< Array for storing the mapping from columns of \f$ M \f$ to columns of \f$ M_2 \f$
                                 **  or to \c SIZE_MAX; may be \c NULL. */
  size_t* secondSpecialRows,    /**< Array of length 3 for storing the row indices of the artificial row, of
                                 ** \f$ C_{i,\star} \f$ and of \f$ C_{j,\star} \f$ in \f$ M_2 \f$; may be \c NULL. */
  size_t* secondSpecialColumns  /**< Array of length 2 for storing the column indices of \f$ C_{\star,k} \f$ and
                                 **  of \f$ C_{\star,\ell} \f$ in \f$ M_2 \f$; may be \c NULL. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_SEPARATION_H */
