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
 * \brief Constructs the Seymour 3-sum of the two matrices \p first and \p second via \p firstMarker1, \p firstMarker2,
 *        \p firstMarker3, \p secondMarker1, \p secondMarker2 and \p secondMarker3.
 *
 * Let \f$ M_1 \f$ and \f$ M_2 \f$ denote the matrices given by \p first and \p second, let \f$ A \f$ be the matrix
 * \f$ M_1 \f$ without the rows and columns indexed by \p firstMarker1, \p firstMarker2 and \p firstMarker3, and let
 * \f$ D \f$ be the matrix \f$ M_2 \f$ without the rows and columns indexed by \p secondMarker1, \p secondMarker2 and
 * \p secondMarker3. Exactly one of the first markers must index a row and the other two (distinct) columns of \p first.
 * After reordering these to be last, \f$ M_1 \f$ must be of the form
 * \f$
 *   M_1 = \begin{bmatrix}
 *     A & a & a \\
 *     c^{\textsf{T}} & 0 & \varepsilon
 *   \end{bmatrix},
 * \f$
 * where \f$ \varepsilon \in \{-1,+1 \} \f$.
 * Similarly, exactly one of the second markers must index a row and the other two (distinct) columns of \p second.
 * After reordering these to be first, \f$ M_2 \f$ must be of the form
 * \f$
 *   M_2 = \begin{bmatrix}
 *     \varepsilon & 0 & b^{\textsf{T}} \\
 *     d & d & D
 *   \end{bmatrix}
 * \f$
 * with the same \f$ \varepsilon \f$.
 * The 3-sum of \f$ M_1 \f$ and \f$ M_2 \f$ (at the markers) is the matrix
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
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* first,          /**< First matrix. */
  CMR_CHRMAT* second,         /**< Second matrix. */
  CMR_ELEMENT firstMarker1,   /**< 1st marker element of first matrix. */
  CMR_ELEMENT firstMarker2,   /**< 2nd marker element of first matrix. */
  CMR_ELEMENT firstMarker3,   /**< 3rd marker element of first matrix. */
  CMR_ELEMENT secondMarker1,  /**< 1st marker element of second matrix. */
  CMR_ELEMENT secondMarker2,  /**< 2nd marker element of second matrix. */
  CMR_ELEMENT secondMarker3,  /**< 3rd marker element of second matrix. */
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
 * Consequently, \f$ b^{\textsf{T}} \f$ and \f$ d \f$ are possibly negated rows and columns of \f$ M \f$.
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
 * Consequently, \f$ b^{\textsf{T}} \f$ and \f$ d \f$ are possibly negated rows and columns of \f$ M \f$.
 *
 * The value of \f$ \varepsilon \in \{-1,+1\} \f$ must be given by \p epsilon and should be computed by
 * \ref CMRthreeSumSeymourDecomposeEpsilon.
 *
 * This function computes \f$ M_1 \f$, while \f$ M_2 \f$ can be computed by \ref CMRthreeSumSeymourDecomposeSecond.
 *
 * \note If \p pfirstMarker1 is not \c NULL, then \p *pfirstMarker1 will always refer to the last row.
 * \note If \p pfirstMarker2 is not \c NULL, then \p *pfirstMarker2 will always refer to the second-to-last column.
 * \note If \p pfirstMarker3 is not \c NULL, then \p *pfirstMarker3 will always refer to the last column
 */

CMR_EXPORT
CMR_ERROR CMRthreeSumSeymourDecomposeFirst(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,           /**< Input matrix \f$ M \f$. */
  CMR_SEPA* sepa,               /**< 3-separation to decompose at. */
  char epsilon,                 /**< Value of \f$ \varepsilon \f$. */
  CMR_CHRMAT** pfirst,          /**< Pointer for storing the first matrix \f$ M_1 \f$. */
  size_t* firstRowsOrigin,      /**< Array for storing the mapping from rows of \f$ M_1 \f$ to rows of \f$ M \f$;
                                 **  also set for the extra row if applicable, even if negated; may be \c NULL. */
  size_t* firstColumnsOrigin,   /**< Array for storing the mapping from columns of \f$ M_1 \f$ to columns of \f$ M \f$;
                                 **  also set for the extra column if applicable, even if negated; may be \c NULL. */
  size_t* rowsToFirst,          /**< Array for storing the mapping from rows of \f$ M \f$ to rows of \f$ M_1 \f$ or to
                                 **  \c SIZE_MAX; may be \c NULL. */
  size_t* columnsToFirst,       /**< Array for storing the mapping from columns of \f$ M \f$ to columns of \f$ M_1 \f$
                                 **  or to \c SIZE_MAX; may be \c NULL. */
  CMR_ELEMENT* pfirstMarker1,   /**< Pointer for storing the row index of \f$ c^{\textsf{T}} \f$ in \f$ M_1 \f$;
                                 **  may be \c NULL. */
  CMR_ELEMENT* pfirstMarker2,   /**< Pointer for storing the column index of the 1st copy of \f$ a \f$ in \f$ M_1 \f$;
                                 **  may be \c NULL. */
  CMR_ELEMENT* pfirstMarker3    /**< Pointer for storing the column index of the 2nd copy of \f$ a \f$ in \f$ M_1 \f$;
                                 **  may be \c NULL. */
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
 * Consequently, \f$ b^{\textsf{T}} \f$ and \f$ d \f$ are possibly negated rows and columns of \f$ M \f$.
 *
 * The value of \f$ \varepsilon \in \{-1,+1\} \f$ must be given by \p epsilon and should be computed by
 * \ref CMRthreeSumSeymourDecomposeEpsilon.
 *
 * This function computes \f$ M_2 \f$, while \f$ M_1 \f$ can be computed by \ref CMRthreeSumSeymourDecomposeFirst.
 *
 * \note If \p psecondMarker1 is not \c NULL, then \p *psecondMarker1 will always refer to the 1st row.
 * \note If \p psecondMarker2 is not \c NULL, then \p *psecondMarker2 will always refer to the 1st column.
 * \note If \p psecondMarker3 is not \c NULL, then \p *psecondMarker3 will always refer to the 2nd column
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
  CMR_ELEMENT* psecondMarker1,  /**< Pointer for storing the row index of \f$ b^{\textsf{T}} \f$ in \f$ M_2 \f$;
                                 **  may be \c NULL. */
  CMR_ELEMENT* psecondMarker2,  /**< Pointer for storing the column index of the 1st copy of \f$ d \f$ in \f$ M_2 \f$;
                                 **  may be \c NULL. */
  CMR_ELEMENT* psecondMarker3   /**< Pointer for storing the column index of the 2nd copy of \f$ d \f$ in \f$ M_2 \f$;
                                 **  may be \c NULL. */
);

/**
 * \brief Constructs the Truemper 3-sum of the two matrices \p first and \p second via \p firstRows, \p firstColumns,
 *        \p secondRows and \p secondColumns.
 *
 * Let \f$ M_1 \f$ and \f$ M_2 \f$ denote the matrices given by \p first and \p second, let \f$ A \f$ be the matrix
 * \f$ M_1 \f$ without the rows \p firstRows[0] and \p firstRows[1] and column \p firstColumns[2], and let \f$ D \f$ be
 * the matrix \f$ M_2 \f$ without the row \p secondRows[0] and columns \p secondColumns[0] and \p secondColumns[1].
 * After permuting rows \p firstRows[0] and \p firstRows[1] to be the last two rows and permuting column
 * \p firstColumns[2] to be the last, \f$ M_1 \f$ must be of the form
 * \f[
 *   M_1 = \begin{bmatrix}
 *     A & \mathbb{O} \\
 *     C_{i,\star} & \alpha \\
 *     C_{j,\star} & \beta
 *   \end{bmatrix},
 * \f]
 * where \f$ \alpha,\beta \in \{-1,+1 \} \f$.
 * Similarly, after reordering the rows / columns of \p second such that \p secondRows[0] is first and columns
 * \p secondColumns[0] and \p secondColumns[1] are first, \f$ M_2 \f$ must be of the form
 * \f[
 *   M_2 = \begin{bmatrix}
 *     \gamma & \delta & \mathbb{O}^{\textsf{T}} \\
 *     C_{\star,k} & C_{\star,\ell} & D
 *   \end{bmatrix},
 * \f]
 * where \f$ \gamma,\delta \in \{ -1,+1 \} \f$ and such that the matrix
 * \f[
 *   N = \begin{bmatrix}
 *     \gamma & \delta & 0 \\
 *     C_{i,k} & C_{i,\ell} & \alpha \\
 *     C_{j,k} & C_{j,\ell} & \beta
 *   \end{bmatrix}
 * \f]
 * is totally unimodular. The columns \p firstColumns[0] and \p firstColumns[1] indicate the columns of \f$ M_1 \f$
 * that shall correspond to \f$ C_{\star,k}\f$ and \f$ C_{\star,\ell} \f$, respectively.
 * Similarly, the rows \p secondRows[1] and \p secondRows[2] indicate the rows of \f$ M_2 \f$
 * that shall correspond to \f$ C_{i,\star}\f$ and \f$ C_{j,\star} \f$, respectively.
 *
 * \note The 2-by-2 submatrix of \f$ M_1 \f$ indexed by rows \p firstRows[0] and \p firstRows[1] and columns
 *       \p firstColumns[0] and \p firstColumns[1] must be identical to the submatrix of \f$ M_2 \f$ indexed by rows
 *       \p secondRows[1] and \p secondRows[2] and columns \p secondColumns[0] and \p secondColumns[1], which is the
 *       matrix \f$ C_{\{i,j\},\{k,\ell\}} \f$.
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
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* first,          /**< First matrix. */
  CMR_CHRMAT* second,         /**< Second matrix. */
  size_t firstRows[2],        /**< Array of last two rows of connecting submatrix in \p first. */
  size_t firstColumns[3],     /**< Array of all columns of connecting submatrix in \p first. */
  size_t secondRows[3],       /**< Array of all rows of connecting submatrix in \p second. */
  size_t secondColumns[2],    /**< Array of first two columns of connecting submatrix in \p second. */
  int8_t characteristic,      /**< Field characteristic. */
  CMR_CHRMAT** presult        /**< Pointer for storing the result. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_SEPARATION_H */
