#ifndef TU_REGULAR_INTERNAL_H
#define TU_REGULAR_INTERNAL_H

#include <tu/regular.h>

/**
 * \brief Creates a decomposition tree node.
 */
void TUcreateDec(
  TU* tu,       /**< TU environment. */
  TU_DEC** dec  /**< Pointer to decomposition tree. */
);

/**
 * \brief Creates a decomposition tree node from as minor.
 */
void TUcreateDecChild(
  TU* tu,                       /**< TU environment. */
  TU_DEC* dec,                  /**< Given decomposition node. */
  int numRows,                  /**< Number of rows of the child. */
  int* rows,                    /**< Rows of child node's matrix. */
  int numColumns,               /**< Number of columns of the child. */
  int* columns,                 /**< Columns of child node's matrix. */
  int numNonzeros,              /**< Number of nonzeros in parent rows/columns or \c 0. */
  int numExtraRows,             /**< Extra rows for the child matrix without a parent row. */
  int numExtraColumns,          /**< Extra columns for the child matrix without a parent row. */
  int numExtraNonzeros,         /**< Number of nonzeros in extra rows and columns. */
  bool constructDecomposition,  /**< Whether to construct a proper decomposition. */
  TU_DEC** result               /**< Pointer to decomposition child. */
);



/**
 * \brief Tests binary linear matroid for regularity. The given decomposition is initialized with
 *        a sequentially connected matrix.
 */

bool TUregularSequentiallyConnected(
  TU* tu,                 /**< TU environment. */
  TU_DEC* decomposition,  /**< Partial decomposition. */
  bool certify,           /**< Whether to compute a proper decomposition tree. */
  bool notGraphic,        /**< If \c true, the matroid is known to be not graphic. */
  bool notCographic       /**< If \c true, the matroid is known to be not cographic. */
);

#endif /* TU_REGULAR_INTERNAL_H */

