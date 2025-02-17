#ifndef CMR_BIPARTITE_GRAPH_INTERNAL_H
#define CMR_BIPARTITE_GRAPH_INTERNAL_H

#include <cmr/element.h>
#include <cmr/matrix.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Finds a shortest path between different vertex groups in the bipartite graph of a submatrix of the given
 *        \p matrix.
 *
 * The bipartite graph of a matrix \f$ M \f$ has the rows and columns of as vertices and edges for all nonzeros of
 * \f$ M \f$. The rows and columns are assigned to groups, and some shortest path from any nonzero group to any larger
 * nonzero group will be returned. Vertices whose group is negative are ignored.
 *
 * A negative row/column group value means to disable the node. Positive values indicate different groups.
 */

CMR_ERROR CMRchrmatSubmatrixBipartitePath(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,               /**< Matrix. */
  CMR_CHRMAT* transpose,            /**< Transpose of \p matrix. */
  int* rowsGroup,                   /**< Array that specifies each row's group. */
  int* columnsGroup,                /**< Array that specifies each column's group. */
  bool* pconnected,                 /**< Pointer for storing whether such a path exists; may be \c NULL. */
  CMR_ELEMENT* ppathSource,         /**< Pointer for storing the source row/column; set to invalid if no path exists;
                                     **  may be \c NULL. */
  CMR_ELEMENT* ppathTarget,         /**< Pointer for storing the target row/column; set to invalid if no path exists;
                                     **  may be \c NULL. */
  CMR_ELEMENT* rowsPredecessor,     /**< Array for storing the predecessor of each row vertex; may be \c NULL. */
  CMR_ELEMENT* columnsPredecessor,  /**< Array for storing the predecessor of each column vertex; may be \c NULL. */
  int* psum                         /**< Pointer for storing the sum of the edge entries on the path; may be \c NULL. */
);

#ifdef __cplusplus
}
#endif

#endif /* CMR_BIPARTITE_GRAPH_INTERNAL_H */
