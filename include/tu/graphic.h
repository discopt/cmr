#ifndef TU_GRAPHIC_H
#define TU_GRAPHIC_H

#include <tu/env.h>
#include <tu/matrix.h>
#include <tu/graph.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \brief Computes the binary representation matrix M for a given graph.
 */

TU_EXPORT
TU_ERROR TUcomputeGraphBinaryRepresentationMatrix(
  TU* tu,                       /**< \ref TU environment. */
  TU_GRAPH* graph,              /**< Graph. */
  TU_CHRMAT** pmatrix,          /**< Pointer for storing M (may be be \c NULL). */
  TU_CHRMAT** ptranspose,       /**< Pointer for storing the transpose of M (may be be \c NULL). */
  int numBasisEdges,            /**< Length of \p basisEdges (0 if \c basisEdges is NULL). */
  TU_GRAPH_EDGE* basisEdges,    /**< If not \c NULL, tries to use these edges for the basis. */
  int numCobasisEdges,          /**< Length of \p cobasisEdges (0 if \c cobasisEdges is NULL). */
  TU_GRAPH_EDGE* cobasisEdges,  /**< If not \c NULL, tries to order columns as specified. */
  bool* pisCorrectBasis         /**< If not \c NULL, returns \c true if and only if \c basisEdges formed a basis. */
);

/**
 * \brief Computes the ternary representation matrix M for a graph.
 */

TU_EXPORT
TU_ERROR TUcomputeGraphTernaryRepresentationMatrix(
  TU* tu,                       /**< \ref TU environment. */
  TU_GRAPH* graph,              /**< Graph. */
  TU_CHRMAT** pmatrix,          /**< Pointer for storing M (may be be \c NULL). */
  TU_CHRMAT** ptranspose,       /**< Pointer for storing the transpose of M (may be be \c NULL). */
  bool* edgesReversed,          /**< Indicates, for each edge {u,v}, whether we consider (u,v) (if \c false) */
                                /**< or (v,u) (if \c true). */
  int numBasisEdges,            /**< Length of \p basisEdges (0 if \c basisEdges is \c NULL). */
  TU_GRAPH_EDGE* basisEdges,    /**< If not \c NULL, tries to use these edges for the basis. */
  int numCobasisEdges,          /**< Length of \p cobasisEdges (0 if \c cobasisEdges is \c NULL). */
  TU_GRAPH_EDGE* cobasisEdges,  /**< If not \c NULL, tries to order columns as specified. */
  bool* pisCorrectBasis         /**< If not \c NULL, returns \c true if and only if \c basisEdges formed a basis. */
);

/**
 * \brief Tests a binary matrix for graphicness.
 * 
 * The binary matrix is given by its \p transpose. It determines whether the matrix is the binary representation matrix
 * of some graph and sets \p *pisGraphic accordingly. The algorithm is based on one of Bixby and Wagner.
 *
 * If \p pgraph is not \c NULL and if the matrix is graphic, then *\p pgraph will point to a graph representing it, and
 * it is set to \c NULL otherwise. The caller must release the memory via \ref TUgraphFree.
 * If in addition to \p pgraph also \p pbasis (resp. \p pcobasis) is not \c NULL, then this will be an array of edges
 * forming the basis (resp. cobasis), i.e., indexing the rows (resp. columns) of the matrix.
 * The caller must release the memory via \ref TUfreeBlockArray.
 *
 * If \p psubmatrix is not \c NULL and the matrix is not graphic, then a minimal nongraphic submatrix will be searched,
 * which may cause extra computational effort. In this case, *\p submatrix will point to this submatrix. The caller must
 * release the memory via \ref TUsubmatFree. It is set to \c NULL otherwise.
 */

TU_EXPORT
TU_ERROR TUtestBinaryGraphic(
  TU* tu,                   /**< \ref TU environment. */
  TU_CHRMAT* transpose,     /**< Transpose of matrix to be tested. */
  bool* pisGraphic,         /**< Returns true if and only if the matrix is graphic. */
  TU_GRAPH** pgraph,        /**< Pointer for storing the graph (if graphic). */
  TU_GRAPH_EDGE** pbasis,   /**< Pointer for storing the basis (if graphic), i.e., a spanning forest.  */
  TU_GRAPH_EDGE** pcobasis, /**< Pointer for storing the cobasis (if graphic), i.e., complement of \p basis. */
  TU_SUBMAT** psubmatrix    /**< Pointer for storing a minimal nongraphic submatrix (if nongraphic). */
);

/**
 * \brief Tests a ternary matrix for graphicness.
 * 
 * The ternary matrix is given by its \p transpose. It determines whether the matrix is the ternary representation
 * matrix of some graph and sets \p *pisGraphic accordingly. The algorithm is based on one of Bixby and Wagner.
 *
 * If \p pgraph is not \c NULL and if the matrix is graphic, then *\p pgraph will point to a graph representing it, and
 * it is set to \c NULL otherwise. The caller must release the memory via \ref TUgraphFree.
 * If in addition to \p pgraph also \p pbasis (resp. \p pcobasis) is not \c NULL, then this will be an array of edges
 * forming the basis (resp. cobasis), i.e., indexing the rows (resp. columns) of the matrix.
 * Similarly, if \p pgraph is not \c NULL and if the matrix is graphic, then \p pedgesReversed will be an array that
 * indicates, for each edge, whether this shall be reversed.
 * The caller must release the memory of \p *pbasis, \p *pcobasis and \p *pedgesReversed via \ref TUfreeBlockArray.
 *
 * If \p psubmatrix is not \c NULL and the matrix is not graphic, then a minimal nongraphic submatrix will be searched,
 * which may cause extra computational effort. In this case, *\p submatrix will point to this submatrix. The caller must
 * release the memory via \ref TUsubmatFree. It is set to \c NULL otherwise.
 */

TU_EXPORT
TU_ERROR TUtestTernaryGraphic(
  TU* tu,                   /**< \ref TU environment. */
  TU_CHRMAT* transpose,     /**< Transpose of matrix to be tested. */
  bool* pisGraphic,         /**< Returns true if and only if the matrix is graphic. */
  TU_GRAPH** pgraph,        /**< Pointer for storing the graph (if graphic). */
  TU_GRAPH_EDGE** pbasis,   /**< Pointer for storing the basis (if graphic), i.e., a spanning forest.  */
  TU_GRAPH_EDGE** pcobasis, /**< Pointer for storing the cobasis (if graphic), i.e., complement of \p basis. */
  bool** pedgesReversed,    /**< Pointer for storing indicators which edges are reversed for the correct sign. */
  TU_SUBMAT** psubmatrix    /**< Pointer for storing a minimal nongraphic submatrix (if nongraphic). */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_GRAPHIC_H */
