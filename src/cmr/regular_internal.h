#ifndef CMR_REGULAR_INTERNAL_H
#define CMR_REGULAR_INTERNAL_H

#include <cmr/regular.h>

/**
 * \brief Enumerates 3-separations for a 3-connected matrix.
 */

CMR_ERROR CMRregularSearchThreeSeparation(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_DEC* dec,                   /**< Decomposition node. */
  CMR_CHRMAT* transpose,          /**< Transpose of nested-minors matrix of \p dec. */
  bool ternary,                   /**< Whether to consider the signs of the matrix. */
  size_t firstNonCoGraphicMinor,  /**< Index of first nested minor that is neither graphic nor cographic. */
  CMR_SUBMAT** psubmatrix,        /**< Pointer for storing a violator matrix. */
  CMR_REGULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                /**< Time limit to impose. */
);

/**
 * \brief Tests whether given 3-connected matrix represents \f$ R_{10} \f$.
 */

CMR_ERROR CMRregularThreeConnectedIsR10(
  CMR* cmr,     /**< \ref CMR environment. */
  CMR_DEC* dec, /**< Decomposition. */
  bool* pisR10  /**< Pointer for storing whether \p matrix represents \f$ R_{10} \f$. */
);

/**
 * \brief Tests sequence of nested 3-connected minors for graphicness.
 */

CMR_ERROR CMRregularSequenceGraphic(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,             /**< Matrix. */
  CMR_CHRMAT* transpose,          /**< Transpose. */
  size_t lengthSequence,          /**< Length of the sequence of nested minors. */
  size_t* sequenceNumRows,        /**< Array with number of rows of each minor. */  
  size_t* sequenceNumColumns,     /**< Array with number of columns of each minor. */  
  size_t* plastGraphicMinor,      /**< Pointer for storing the last graphic minor. */
  CMR_GRAPH** pgraph,             /**< Pointer for storing the graph. */
  CMR_ELEMENT** pedgeElements,    /**< Pointer for storing the mapping of edges to elements. */
  CMR_REGULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                /**< Time limit to impose. */
);

/**
 * \brief Extends an incomplete sequence of nested 3-connected minors for the matrix of a decomposition node.
 *
 * In case the matrix is not 3-connected, a 2-separation is applied to \p dec and the function terminates, filling
 * the relevant variables of \p dec.
 */

CMR_ERROR CMRregularExtendNestedMinorSequence(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_DEC* dec,                   /**< Decomposition node. */
  bool ternary,                   /**< Whether to consider the signs of the matrix. */
  CMR_SUBMAT** psubmatrix,        /**< Pointer for storing a violator matrix. */
  CMR_REGULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                /**< Time limit to impose. */
);

/**
 * \brief Constructs a sequence of nested 3-connected minors for the matrix of a decomposition node.
 *
 * In case the matrix is not 3-connected, a 2-separation is applied to \p dec and the function terminates, filling
 * the relevant variables of \p dec.
 */

CMR_ERROR CMRregularConstructNestedMinorSequence(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_DEC* dec,                   /**< Decomposition node. */
  bool ternary,                   /**< Whether to consider the signs of the matrix. */
  CMR_SUBMAT* wheelSubmatrix,     /**< Wheel submatrix to start with. */
  CMR_SUBMAT** psubmatrix,        /**< Pointer for storing a violator matrix. */
  CMR_REGULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                /**< Time limit to impose. */
);

/**
 * \brief Tests matrix for graphicness.
 */

CMR_ERROR CMRregularTestGraphic(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_CHRMAT** pmatrix,           /**< Pointer to matrix. */
  CMR_CHRMAT** ptranspose,        /**< Pointer to transpose. */
  bool ternary,                   /**< Whether to also check the signs of the matrix. */
  bool* pisGraphic,               /**< Pointer for storing the result. */
  CMR_GRAPH** pgraph,             /**< Pointer for storing the graph if the matrix is graph (may be \c NULL). */
  CMR_GRAPH_EDGE** pforest,       /**< Pointer for storing the mapping of rows to forest edges. */
  CMR_GRAPH_EDGE** pcoforest,     /**< Pointer for storing the mapping of rows to forest edges. */
  bool** parcsReversed,           /**< Pointer for storing the array indicating which arcs are reversed. */ 
  CMR_SUBMAT** psubmatrix,        /**< Pointer for storing a minimal non-graphic submatrix (may be \c NULL). */
  CMR_REGULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                /**< Time limit to impose. */
);

/**
 * \brief Splits off series-parallel elements from the matrix of a decomposition node.
 *
 * In case the matrix is \ref series-parallel, then \p *pmatrix is set to \c NULL and \p *pdec is declared to be planar.
 *
 * In case the matrix does not admit series-parallel reductions, then \p *pmatrix and \p *pdec remain unchanged, except
 * The decomposition node \p *pdec and matrix \p *pmatrix are replaced by the SP-reduced ones, i.e., by the
 * corresponding (grand-) children of the given decomposition and the corresponding matrix.
 */

CMR_ERROR CMRregularDecomposeSeriesParallel(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_DEC** pdec,                 /**< Pointer to decomposition node. */
  bool ternary,                   /**< Whether to consider the signs of the matrix. */
  CMR_SUBMAT** psubmatrix,        /**< Pointer for storing a violator matrix. */
  CMR_REGULAR_PARAMETERS* params, /**< Parameters for the computation. */
  CMR_REGULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                /**< Time limit to impose. */
);

/**
 * \brief Performs a 1-sum decomposition of \p matrix and stores it in \p dec.
 *
 * If \p matrix is 1-connected, then \p dec remains unchanged. In particular, the \c matrix and \c transpose members remain
 * \c NULL. Otherwise, \p dec will become a \ref CMR_DEC_ONE_SUM node with children that are initialized to the
 * 1-connected components. In this case, the \c matrix and \c transpose members of the child nodes are set.
 */

CMR_ERROR CMRregularDecomposeOneSum(
  CMR* cmr,     /**< \ref CMR environment. */
  CMR_DEC* dec  /**< Initialized decomposition node with \ref matrix attribute. */
);

/**
 * \brief Tests ternary or binary linear matroid for regularity.
 *
 * If \p pdec is not \c NULL, \c *pdec will be a (partial) decomposition tree.
 *
 * If \p pminor is not \c NULL and \p matrix is not regular, then an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor or a
 * submatrix with non-ternary determinant is searched. This causes additional computational effort!
 */

CMR_ERROR CMRtestRegular(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_CHRMAT* matrix,             /**< Input matrix. */
  bool ternary,                   /**< Whether signs of \p matrix play a role. */
  bool *pisRegular,               /**< Pointer for storing whether \p matrix is regular. */
  CMR_DEC** pdec,                 /**< Pointer for storing the decomposition tree (may be \c NULL). */
  CMR_MINOR** pminor,             /**< Pointer for storing an \f$ F_7 \f$ or \f$ F_7^\star \f$ minor. */
  CMR_REGULAR_PARAMETERS* params, /**< Parameters for the computation. */
  CMR_REGULAR_STATISTICS* stats,  /**< Statistics for the computation (may be \c NULL). */
  double timeLimit                /**< Time limit to impose. */
);

#endif /* CMR_REGULAR_INTERNAL_H */
