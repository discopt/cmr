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

typedef struct
{
  unsigned char* rowsToPart;        /**< \brief Indicates to which block each row belongs. Values above 1 are ignored. */
  unsigned char* columnsToPart;     /**< \brief Indicates to which block each column belongs. Values above 1 are ignored. */
  size_t numRows[2];                /**< \brief Indicates the number of rows of each part. */
  size_t numColumns[2];             /**< \brief Indicates the number of columns of each part. */
  size_t* rows[2];                  /**< \brief Array of sorted rows for each part. */
  size_t* columns[2];               /**< \brief Array of sorted columns for each part. */
  size_t extraRows[2][2];           /**< \brief For each part, array of extra rows; may be \c SIZE_MAX. */
  size_t extraColumns[2][2];        /**< \brief For each part, array of extra columns; may be \c SIZE_MAX. */
  unsigned char* indicatorMemory;   /**< \brief Memory for \ref rowsToPart and \ref columnsToPart. */
  size_t* elementMemory;            /**< \brief Memory for \ref rows and \ref columns. */
} CMR_SEPA;

/**
 * \brief Creates a separation.
 *
 * Only the memory is allocated. The usualy way to initialize it is to fill the arrays \ref rowsToPart and
 * \ref columnsToPart and then call \ref CMRsepaInitialize.
 */

CMR_EXPORT
CMR_ERROR CMRsepaCreate(
  CMR* cmr,           /**< \ref CMR environment. */
  size_t numRows,     /**< Number of rows. */
  size_t numColumns,  /**< Number of columns. */
  CMR_SEPA** psepa    /**< Pointer for storing the created separation. */
);

/**
 * \brief Initializes a separation.
 *
 * Assumes that \p separation was created via \ref CMRsepaCreate and that all entries of \ref rowsToPart and
 * \ref columnsToPart are set to either 0 or 1.
 * The connecting elements are given by the caller.
 */

CMR_EXPORT
CMR_ERROR CMRsepaInitialize(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_SEPA* separation,       /**< Already created separation. */
  size_t firstExtraRow0,      /**< First extra row for part 0 or \c SIZE_MAX if bottom-left rank is 0. */
  size_t firstExtraColumn1,   /**< First extra column for part 1 or \c SIZE_MAX if bottom-left rank is 0. */
  size_t firstExtraRow1,      /**< First extra row for part 1 or \c SIZE_MAX if top-right rank is 0. */
  size_t firstExtraColumn0,   /**< First extra column for part 0 or \c SIZE_MAX if top-right rank is 0. */
  size_t secondExtraRow0,     /**< Second extra row for part 0 or \c SIZE_MAX if bottom-left rank is at most 1. */
  size_t secondExtraColumn1,  /**< Second extra column for part 1 or \c SIZE_MAX if bottom-left rank is at most 1. */
  size_t secondExtraRow1,     /**< Second extra row for part 1 or \c SIZE_MAX if top-right rank is at most 1. */
  size_t secondExtraColumn0   /**< Second extra column for part 0 or \c SIZE_MAX if top-right rank is at most 1. */
);

/**
 * \brief Initializes a separation.
 *
 * Assumes that \p separation was created via \ref CMRsepaCreate and that all entries of \ref rowsToPart and
 * \ref columnsToPart are set to either 0 or 1.
 * The connecting elements will be searched by inspecting \p matrix.
 */

CMR_EXPORT
CMR_ERROR CMRsepaInitializeMatrix(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEPA* separation,   /**< Already created separation. */
  CMR_CHRMAT* matrix,     /**< Matrix this separation belongs to. */
  unsigned char totalRank /**< Total rank of separation. */
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
 * \brief Returns rank of bottom-left submatrix.
 */

static inline
unsigned char CMRsepaRankBottomLeft(
  CMR_SEPA* sepa  /**< Separation. */
)
{
  assert(sepa);

  return sepa->extraRows[0][0] == SIZE_MAX ? 0
    : (sepa->extraRows[0][1] == SIZE_MAX ? 1 : 2);
}

/**
 * \brief Returns rank of top-right submatrix.
 */

static inline
unsigned char CMRsepaRankTopRight(
  CMR_SEPA* sepa  /**< Separation. */
)
{
  return sepa->extraRows[1][0] == SIZE_MAX ? 0
    : (sepa->extraRows[1][1] == SIZE_MAX ? 1 : 2);
}

/**
 * \brief Returns rank sum of bottom-left and top-right submatrices.
 */

static inline
unsigned char CMRsepaRank(
  CMR_SEPA* sepa  /**< Separation. */
)
{
  return CMRsepaRankBottomLeft(sepa) + CMRsepaRankTopRight(sepa);
}

/**
 * \brief Checks for a given matrix whether the binary k-separation is also a ternary one.
 *
 * Checks, for a ternary input matrix \f$ M \f$ and a k-separation (\f$ k \in \{1,2,3\} \f$) of the (binary) support
 * matrix of \f$ M \f$, whether it is also a k-separation of \f$ M \f$ itself. The result is stored in \p *pisTernary.
 *
 * If the check fails, a certifying submatrix is returned.
 */

CMR_EXPORT
CMR_ERROR CMRsepaCheckTernary(
  CMR* cmr,               /**< \ref CMR environment. */
  CMR_SEPA* sepa,         /**< Separation. */
  CMR_CHRMAT* matrix,     /**< Matrix. */
  CMR_SUBMAT* submatrix,  /**< Submatrix of \p matrix to consider (may be \c NULL). */
  bool* pisTernary,       /**< Pointer for storing whether the check passed. */
  CMR_SUBMAT** psubmatrix /**< Pointer for storing a violator submatrix (may be \c NULL). */
);

/**
 * \brief Constructs the 1-sum of the two matrices \p first and \p second.
 *
 * Let \f$ A \f$ and \f$ B \f$ denote the matrices given by \p first and \p second.
 * Their 2-sum is the matrix
 * \f[
 *   C := \begin{pmatrix}
 *     A & \mathbb{O} \\
 *     \mathbb{O} & B
 *   \end{pmatrix}.
 * \f]
 * The resulting matrix \f$ C \f$ is created and stored in \p *presult.
 */

CMR_EXPORT
CMR_ERROR CMRoneSum(
  CMR* cmr,             /**< \ref CMR environment. */
  CMR_CHRMAT* first,    /**< First matrix. */
  CMR_CHRMAT* second,   /**< Second matrix. */
  CMR_CHRMAT** presult  /**< Pointer for storing the result. */
);

/**
 * \brief Constructs the 2-sum of the two matrices \p first and \p second via \p firstMarker and \p secondMarker.
 *
 * Let \f$ A \f$ and \f$ B \f$ denote the matrices given by \p first and \p second and let \f$ A' \f$ and \f$ B' \f$ be
 * these matrices without the row or column indexed by \p firstMarker and \p secondMarker, respectively.
 * If \p firstMarker indexes a row vector \f$ a^{\textsf{T}} \f$ of \f$ A \f$ then \p secondMarker must index a column
 * vector \f$ b \f$ of \f$ B \f$. In this case the 2-sum is the matrix
 * \f[
 *   C := \begin{pmatrix}
 *     A' & \mathbb{O} \\
 *     b a^{\textsf{T}} & B'
 *   \end{pmatrix}.
 * \f]
 * Otherwise, \p firstMarker must index a column vector \f$ a \f$ of \f$ A \f$ and \p secondMarker must index a row
 * vector \f$ b^{\textsf{T}} \f$ of \f$ B \f$, and the 2-sum is the matrix
 * \f[
 *   C := \begin{pmatrix}
 *     A' & a b^{\textsf{T}} \\
 *     \mathbb{O} & B'
 *   \end{pmatrix}.
 * \f]
 * The calculations are done modulo \p characteristic, where the value \f$ 3 \f$ yields numbers from \f$ \{-1,0,+1\} \f$.
 *
 * The resulting matrix \f$ C \f$ is created and stored in \p *presult.
 */

CMR_EXPORT
CMR_ERROR CMRtwoSum(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* first,          /**< First matrix. */
  CMR_CHRMAT* second,         /**< Second matrix. */
  CMR_ELEMENT firstMarker,    /**< Marker element of first matrix. */
  CMR_ELEMENT secondMarker,   /**< Marker element of second matrix. */
  int8_t characteristic,      /**< Field characteristic. */
  CMR_CHRMAT** presult        /**< Pointer for storing the result. */
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
 *   C := \begin{pmatrix}
 *     A' & \mathbb{O} \\
 *     b_1 a_1^{\textsf{T}} + b_2 a_2^{\textsf{T}} & B'
 *   \end{pmatrix}.
 * \f]
 * Otherwise, if \p firstMarker1 and \p firstMarker2 both index column vectors \f$ a_1, a_2 \f$ of \f$ A \f$ then
 * \p secondMarker1 and \p secondMarker2 must index row vectors \f$ b_1^{\textsf{T}}, b_2^{\textsf{T}} \f$ of \f$ B \f$.
 * In this case the 3-sum is the matrix
 * \f[
 *   C := \begin{pmatrix}
 *     A'         & a_1 b_1^{\textsf{T}} + a_2 b_2^{\textsf{T}}  \\
 *     \mathbb{O} & B'
 *   \end{pmatrix}.
 * \f]
 * Otherwise, if \p firstMarker1 indexes a row vector \f$ a_1^{\textsf{T}} \f$ and \p firstMarker2 indexes a column
 * vector \f$ a_2 \f$ of \f$ A \f$ then \p secondMarker1 must index a column vector \f$ b_1 \f$ of \f$ B \f$ and
 * \p secondMarker2 must index a row vector \f$ b_2^{\textsf{T}} \f$ of \f$ B \f$.
 * In this case the 3-sum is the matrix
 * \f[
 *   C := \begin{pmatrix}
 *     A'                   & a_2 b_2^{\textsf{T}}  \\
 *     b_1 a_1^{\textsf{T}} & B'
 *   \end{pmatrix}.
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
