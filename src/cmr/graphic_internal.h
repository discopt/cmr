#ifndef CMR_GRAPHIC_INTERNAL_H
#define CMR_GRAPHIC_INTERNAL_H

#include <cmr/graphic.h>

/**
 * \brief Internal method for computing a graphic representation matrix of a graph or a network matrix of a digraph.
 */

CMR_ERROR CMRcomputeRepresentationMatrix(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_GRAPH* graph,               /**< Graph. */
  bool ternary,                   /**< Whether we need to compute correct signs. */
  CMR_CHRMAT** ptranspose,        /**< Pointer for storing the transpose of the matrix. */
  bool* edgesReversed,            /**< Indicates, for each edge {u,v}, whether we consider (u,v) (if \c false)
                                   **< or (v,u) (if \c true). */
  int numForestEdges,             /**< Length of \p forestEdges (0 if \c forestEdges is \c NULL). */
  CMR_GRAPH_EDGE* forestEdges,    /**< If not \c NULL, tries to use these edges for the basis. */
  int numCoforestEdges,           /**< Length of \p coforestEdges (0 if \c coforestEdges is \c NULL). */
  CMR_GRAPH_EDGE* coforestEdges,  /**< If not \c NULL, tries to order columns as specified. */
  bool* pisCorrectForest          /**< If not \c NULL, returns \c true if and only if \c forestEdges are spanning forest. */
);

#endif /* CMR_GRAPHIC_INTERNAL_H */
