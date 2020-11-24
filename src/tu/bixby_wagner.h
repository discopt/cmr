#ifndef TU_BIXBY_WAGNER_H
#define TU_BIXBY_WAGNER_H

#include <tu/graph.h>
#include <tu/matrix.h>

#include "one_sum.h"

bool testGraphicnessBixbyWagner(
  TU* tu,                       /*< TU environment. */
  TU_CHRMAT* matrix,            /*< 1-connected matrix to be test. */
  TU_CHRMAT* transpose,         /*< Transpose of \p matrix. */
  TU_LISTGRAPH* graph,          /*< If not \c NULL and graphic, a graph represented by the matrix. */
  TU_LISTGRAPH_EDGE* basis,     /*< If not \c NULL and graphic, a map from rows to basis edges. */
  TU_LISTGRAPH_EDGE* cobasis,   /*< If not \c NULL and graphic, a map from columns to cobasis edges. */
  TU_SUBMAT** psubmatrix        /*< If not \c NULL and not graphic, a minimal violating submatrix. */
);

#endif /* TU_BIXBY_WAGNER_H */
