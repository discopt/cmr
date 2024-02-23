#ifndef CMR_REGULAR_INTERNAL_H
#define CMR_REGULAR_INTERNAL_H

#include <cmr/regular.h>

#include "matroid_internal.h"

typedef struct DecompositionTask
{
  CMR_MATROID_DEC* dec;           /**< \brief Decomposition node that shall be processed. */
  struct DecompositionTask* next; /**< \brief Next task in queue. */
  CMR_REGULAR_PARAMS* params;     /**< \brief Parameters for the computation. */
  CMR_REGULAR_STATS* stats;       /**< \brief Statistics for the computation (may be \c NULL). */
  clock_t startClock;             /**< \brief Clock for the start time. */
  double timeLimit;               /**< \brief Time limit to impose. */
} DecompositionTask;

/**
 * \brief Creates a decomposition task for the root of the decomposition.
 */

CMR_ERROR CMRregularityTaskCreateRoot(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_MATROID_DEC* dec,           /**< Decomposition node. */
  DecompositionTask** ptask,      /**< Pointer for storing the new task. */
  CMR_REGULAR_PARAMS* params,     /**< Parameters for the computation. */
  CMR_REGULAR_STATS* stats,       /**< Statistics for the computation (may be \c NULL). */
  clock_t startClock,             /**< Clock for the start time. */
  double timeLimit                /**< Time limit to impose. */
);

/**
 * \brief Frees a decomposition task.
 */

CMR_ERROR CMRregularityTaskFree(
  CMR* cmr,                       /**< \ref CMR environment. */
  DecompositionTask** ptask       /**< Pointer to task. */
);

typedef struct DecompositionQueue
{
  DecompositionTask* head;  /**< \brief Next task to be processed. */
  bool foundIrregularity;   /**< \brief Whether irregularity was detected for some node. */
} DecompositionQueue;

/**
 * \brief Initializes a decomposition queue.
 */

CMR_ERROR CMRregularityQueueCreate(
  CMR* cmr,                   /**< \ref CMR environment. */
  DecompositionQueue** pqueue /**< Pointer for storing the queue. */
);

/**
 * \brief Frees the decomposition queue.
 */

CMR_ERROR CMRregularityQueueFree(
  CMR* cmr,                   /**< \ref CMR environment. */
  DecompositionQueue** pqueue /**< Pointer to queue. */
);

/**
 * \brief Returns whether a queue is empty.
 */

bool CMRregularityQueueEmpty(
  DecompositionQueue* queue /**< Queue. */
);

/**
 * \brief Removes a task from a decomposition queue.
 */

DecompositionTask* CMRregularityQueueRemove(
  DecompositionQueue* queue /**< Queue. */
);


/**
 * \brief Adds a task to a decomposition queue.
 */

void CMRregularityQueueAdd(
  DecompositionQueue* queue,  /**< Queue. */
  DecompositionTask* task     /**< Task. */
);

/**
 * \brief Applies a 3-sum decomposition.
 */

CMR_ERROR
CMRregularityDecomposeThreeSum(
  CMR* cmr,                   /**< \ref CMR environment. */
  DecompositionTask* task,    /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionQueue* queue,  /**< Queue of unprocessed nodes. */
  CMR_SEPA* separation        /**< 3-separation. */
);

/**
 * \brief Searches for 3-separations along the sequence of nested minors and decomposes as a 3-sum.
 */

CMR_ERROR
CMRregularityNestedMinorSequenceSearchThreeSeparation(
  CMR* cmr,                 /**< \ref CMR environment. */
  DecompositionTask* task,  /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionQueue* queue /**< Queue of unprocessed nodes. */
);

/**
 * \brief Tests each minor of the sequence of nested 3-connected minors for graphicness.
 *
 * Sets \c task->dec->CMRregularityNestedMinorSequenceGraphicness accordingly.
 */

CMR_ERROR
CMRregularityNestedMinorSequenceGraphicness(
  CMR* cmr,                 /**< \ref CMR environment. */
  DecompositionTask* task,  /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionQueue* queue /**< Queue of unprocessed nodes. */
);

/**
 * \brief Tests each minor of the sequence of nested 3-connected minors for graphicness.
 *
 * Sets \c task->dec->CMRregularityNestedMinorSequenceCographicness accordingly.
 */

CMR_ERROR
CMRregularityNestedMinorSequenceCographicness(
  CMR* cmr,                 /**< \ref CMR environment. */
  DecompositionTask* task,  /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionQueue* queue /**< Queue of unprocessed nodes. */
);

/**
 * \brief Extends an incomplete sequence of nested 3-connected minors for the matrix of a decomposition node.
 *
 * In case the matrix is not 3-connected, a 2-separation is applied to \p dec and the function terminates, filling
 * the relevant variables of \p dec.
 */

CMR_ERROR
CMRregularityExtendNestedMinorSequence(
  CMR* cmr,                 /**< \ref CMR environment. */
  DecompositionTask* task,  /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionQueue* queue /**< Queue of unprocessed nodes. */
);

/**
 * \brief Tests \p matrix for representing \f$ R_{10} \f$.
 */

CMR_ERROR CMRregularityTestR10(
  CMR* cmr,                 /**< \ref CMR environment. */
  DecompositionTask* task,  /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionQueue* queue /**< Queue of unprocessed nodes. */
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
  CMR* cmr,                 /**< \ref CMR environment. */
  DecompositionTask* task,  /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionQueue* queue /**< Queue of unprocessed nodes. */
);

/**
 * \brief Tests \p matrix for graphicness/network and stores it in \p dec.
 */

CMR_ERROR CMRregularityTestGraphicness(
  CMR* cmr,                 /**< \ref CMR environment. */
  DecompositionTask* task,  /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionQueue* queue /**< Queue of unprocessed nodes. */
);

/**
 * \brief Tests \p matrix for cographicness/conetwork and stores it in \p dec.
 */

CMR_ERROR CMRregularityTestCographicness(
  CMR* cmr,                 /**< \ref CMR environment. */
  DecompositionTask* task,  /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionQueue* queue /**< Queue of unprocessed nodes. */
);

/**
 * \brief Performs a 1-sum decomposition of \p matrix and stores it in \p dec.
 *
 * If \p matrix is 1-connected, then \p dec remains unchanged. Otherwise, \p dec will become a
 * \ref CMR_MATROID_DEC_TYPE_ONE_SUM node with children that are initialized to the 1-connected components. In this
 * case, the \c matrix and \c transpose members of the child nodes are set.
 */

CMR_ERROR CMRregularitySearchOneSum(
  CMR* cmr,                 /**< \ref CMR environment. */
  DecompositionTask* task,  /**< Task to be processed; already removed from the list of unprocessed tasks. */
  DecompositionQueue* queue /**< Queue of unprocessed nodes. */
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
