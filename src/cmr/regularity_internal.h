#ifndef CMR_REGULAR_INTERNAL_H
#define CMR_REGULAR_INTERNAL_H

#include <cmr/regular.h>

#include "matroid_internal.h"

typedef struct DecompositionTask
{
  CMR_MATROID_DEC* dec;                     /**< \brief Decomposition node that shall be processed. */

  struct DecompositionTask* next;           /**< \brief Next task in queue. */

  CMR_REGULAR_PARAMS* params;               /**< Parameters for the computation. */
  CMR_REGULAR_STATS* stats;                 /**< \brief Statistics for the computation (may be \c NULL). */
  clock_t startClock;                       /**< \brief Clock for the start time. */
  double timeLimit;                         /**< \brief Time limit to impose. */
} DecompositionTask;

CMR_ERROR CMRregularityTaskCreateRoot(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec,           /**< Decomposition node. */
  DecompositionTask** ptask,      /**< Pointer for storing the new task. */
  CMR_REGULAR_PARAMS* params,     /**< Parameters for the computation. */
  CMR_REGULAR_STATS* stats,       /**< Statistics for the computation (may be \c NULL). */
  clock_t startClock,             /**< Clock for the start time. */
  double timeLimit                /**< Time limit to impose. */
);

CMR_ERROR CMRregularityTaskFree(
  CMR* cmr,                       /**< \ref CMR environment. */
  DecompositionTask** ptask       /**< Pointer to task. */
);



/**
 * \brief Enumerates 3-separations for a 3-connected matrix.
 */

CMR_ERROR CMRregularSearchThreeSeparation(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec,           /**< Decomposition node. */
  CMR_CHRMAT* transpose,          /**< Transpose of nested-minors matrix of \p dec. */
  bool ternary,                   /**< Whether to consider the signs of the matrix. */
  size_t firstNonCoGraphicMinor,  /**< Index of first nested minor that is neither graphic nor cographic. */
  CMR_SUBMAT** psubmatrix,        /**< Pointer for storing a violator matrix. */
  CMR_REGULAR_STATS* stats,       /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                /**< Time limit to impose. */
);

/**
 * \brief Tests sequence of nested 3-connected minors for graphicness.
 */

CMR_ERROR CMRregularSequenceGraphic(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,           /**< Matrix. */
  CMR_CHRMAT* transpose,        /**< Transpose. */
  size_t lengthSequence,        /**< Length of the sequence of nested minors. */
  size_t* sequenceNumRows,      /**< Array with number of rows of each minor. */
  size_t* sequenceNumColumns,   /**< Array with number of columns of each minor. */
  size_t* plastGraphicMinor,    /**< Pointer for storing the last graphic minor. */
  CMR_GRAPH** pgraph,           /**< Pointer for storing the graph. */
  CMR_ELEMENT** pedgeElements,  /**< Pointer for storing the mapping of edges to elements. */
  CMR_REGULAR_STATS* stats,     /**< Statistics for the computation (may be \c NULL). */
  double timeLimit              /**< Time limit to impose. */
);

/**
 * \brief Extends an incomplete sequence of nested 3-connected minors for the matrix of a decomposition node.
 *
 * In case the matrix is not 3-connected, a 2-separation is applied to \p dec and the function terminates, filling
 * the relevant variables of \p dec.
 */

CMR_ERROR CMRregularExtendNestedMinorSequence(
  CMR* cmr,                 /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec,     /**< Decomposition node. */
  bool ternary,             /**< Whether to consider the signs of the matrix. */
  CMR_SUBMAT** psubmatrix,  /**< Pointer for storing a violator matrix. */
  CMR_REGULAR_STATS* stats, /**< Statistics for the computation (may be \c NULL). */
  double timeLimit          /**< Time limit to impose. */
);

/**
 * \brief Constructs a sequence of nested 3-connected minors for the matrix of a decomposition node.
 *
 * In case the matrix is not 3-connected, a 2-separation is applied to \p dec and the function terminates, filling
 * the relevant variables of \p dec.
 */

CMR_ERROR CMRregularConstructNestedMinorSequence(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec,       /**< Decomposition node. */
  bool ternary,               /**< Whether to consider the signs of the matrix. */
  CMR_SUBMAT* wheelSubmatrix, /**< Wheel submatrix to start with. */
  CMR_SUBMAT** psubmatrix,    /**< Pointer for storing a violator matrix. */
  CMR_REGULAR_STATS* stats,   /**< Statistics for the computation (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);

/**
 * \brief Tests \p matrix for representing \f$ R_{10} \f$.
 */

CMR_ERROR CMRregularityTestR10(
  CMR* cmr,                         /**< \ref CMR environment. */
  DecompositionTask* task,          /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionTask** punprocessed  /**< Pointer to head of list of unprocessed tasks. */
);

/**
 * \brief Initializes a sequence of nested minors with one minor, for a given \p wheelSubmatrix.
 */

CMR_ERROR CMRregularityInitNestedMinorSequence(
  CMR* cmr,                   /**< \ref CMR environment. */
  DecompositionTask* task,    /**< Task to be processed; already removed from the list of unprocessed tasks. */
  CMR_SUBMAT* wheelSubmatrix  /**< Wheel submatrix of \p task->dec. */
);

/**
 * \brief Splits off series-parallel elements from the matrix of the decomposition node.
 *
 * In case the matrix is \ref series-parallel, then the node is declared to be planar.
 *
 * In case the matrix does not admit series-parallel reductions, then the node remains remain unchanged.
 * Otherwise, a child node is created whose matrix is the series-parallel-reduced one.
 * For that child node, a 2-separation may be found, and corresponding children will be created.
 */

CMR_ERROR CMRregularityDecomposeSeriesParallel(
  CMR* cmr,                         /**< \ref CMR environment. */
  DecompositionTask* task,          /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionTask** punprocessed  /**< Pointer to head of list of unprocessed tasks. */
);

/**
 * \brief Tests \p matrix for graphicness/network and stores it in \p dec.
 */

CMR_ERROR CMRregularityTestGraphicness(
  CMR* cmr,                         /**< \ref CMR environment. */
  DecompositionTask* task,          /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionTask** punprocessed  /**< Pointer to head of list of unprocessed tasks. */
);

/**
 * \brief Tests \p matrix for cographicness/conetwork and stores it in \p dec.
 */

CMR_ERROR CMRregularityTestCographicness(
  CMR* cmr,                         /**< \ref CMR environment. */
  DecompositionTask* task,          /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionTask** punprocessed  /**< Pointer to head of list of unprocessed tasks. */
);

/**
 * \brief Performs a 1-sum decomposition of \p matrix and stores it in \p dec.
 *
 * If \p matrix is 1-connected, then \p dec remains unchanged. Otherwise, \p dec will become a
 * \ref CMR_MATROID_DEC_TYPE_ONE_SUM node with children that are initialized to the 1-connected components. In this
 * case, the \c matrix and \c transpose members of the child nodes are set.
 */

CMR_ERROR CMRregularitySearchOneSum(
  CMR* cmr,                         /**< \ref CMR environment. */
  DecompositionTask* task,          /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionTask** punprocessed  /**< Pointer to head of list of unprocessed tasks. */
);

/**
 * \brief Tests ternary or binary linear matroid for regularity.
 *
 * If \p pdec is not \c NULL, \c *pdec will be a (partial) decomposition tree.
 *
 * If \p pminor is not \c NULL and \p matrix is not regular, then an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor or a
 * submatrix with non-ternary determinant is searched. This causes additional computational effort!
 */

CMR_ERROR CMRregularityTest(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,         /**< Input matrix. */
  bool ternary,               /**< Whether the matrix shall be considered ternary. */
  bool *pisRegular,           /**< Pointer for storing whether \p matrix is regular. */
  CMR_MATROID_DEC** pdec,     /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_MINOR** pminor,         /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
  CMR_REGULAR_PARAMS* params, /**< Parameters for the computation. */
  CMR_REGULAR_STATS* stats,   /**< Statistics for the computation (may be \c NULL). */
  double timeLimit            /**< Time limit to impose. */
);

#endif /* CMR_REGULAR_INTERNAL_H */
