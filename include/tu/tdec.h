#ifndef TU_TDEC_H
#define TU_TDEC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/graph.h>
#include <tu/matrix.h>

#include "one_sum.h"

typedef enum
{
  TDEC_MEMBER_TYPE_INVALID = 0,
  TDEC_MEMBER_TYPE_BOND = 1,
  TDEC_MEMBER_TYPE_POLYGON = 2,
  TDEC_MEMBER_TYPE_PRIME = 3
} TU_TDEC_MEMBER_TYPE;

typedef int TU_TDEC_EDGE;
typedef int TU_TDEC_NODE;
typedef int TU_TDEC_MEMBER;

struct _TU_TDEC;
typedef struct _TU_TDEC TU_TDEC;
typedef struct _TU_TDEC_NEWCOLUMN TU_TDEC_NEWCOLUMN;

TU_EXPORT
void TUtdecCreate(
  TU* tu,         /*< TU environment. */
  TU_TDEC** ptdec,   /*< Pointer to new t-decomposition. .*/
  int rootRow,    /*< Row of 1-entry of root member. */
  int memEdges,   /*< Initial memory for edges of the t-decomposition. */
  int memNodes,   /*< Initial memory for nodes of the t-decomposition. */
  int memMembers, /*< Initial memory for members of the t-decomposition. */
  int numRows,    /*< Number of rows. */
  int numColumns  /*< Number of columns. */
);

TU_EXPORT
void TUtdecFree(
  TU* tu,         /*< TU environment. */
  TU_TDEC** ptdec /*< Pointer to t-decomposition. .*/
);

TU_EXPORT
int TUtdecBasisSize(
  TU_TDEC* tdec   /*< t-decomposition. */
);

TU_EXPORT
int TUtdecCobasisSize(
  TU_TDEC* tdec   /*< t-decomposition. */
);

TU_EXPORT
int TUtdecNumEdges(
  TU_TDEC* tdec   /*< t-decomposition. */
);

/**
 * \brief Creates a graph represented by given t-decomposition.
 */

TU_EXPORT
void TUtdecToGraph(
  TU* tu,                 /*< TU environment. */
  TU_TDEC* tdec,          /*< t-decomposition. */
  TU_GRAPH* graph,        /*< Graph. */
  bool merge,             /*< Merge and remove corresponding parent and child markers. */
  TU_GRAPH_EDGE* basis,   /*< If not NULL, the edges of a spanning tree are stored here. */
  TU_GRAPH_EDGE* cobasis, /*< If not NULL, the non-basis edges are stored here. */
  int* edgeElements       /*< If not NULL, the elements for each edge are stored here. */
);

TU_EXPORT
void TUtdecnewcolumnCreate(
  TU* tu,                         /*< TU environment. */
  TU_TDEC_NEWCOLUMN** pnewcolumn  /*< new-column structure. */
);

TU_EXPORT
void TUtdecnewcolumnFree(
  TU* tu,                         /*< TU environment. */
  TU_TDEC_NEWCOLUMN** pnewcolumn  /*< new-column structure. */
);

TU_EXPORT
void TUtdecAddColumnCheck(
  TU* tu,                     /*< TU environment. */
  TU_TDEC* tdec,              /*< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcol,  /*< new-column structure. */
  int* entryRows,             /*< Array of rows with 1-entry in this column. */
  int numEntries              /*< Number of 1-entries in this column. */
);

TU_EXPORT
void TUtdecAddColumnApply(
  TU* tu,                     /*< TU environment. */
  TU_TDEC* tdec,              /*< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcol,  /*< new-column structure. */
  int column,                 /*< Index of new column to be added. */
  int* entryRows,             /*< Array of rows with 1-entry in this column. */
  int numEntries              /*< Number of 1-entries in this column. */
);

TU_EXPORT
bool testGraphicnessTDecomposition(
  TU* tu,                 /*< TU environment. */
  TU_CHRMAT* matrix,      /*< 1-connected matrix to be test. */
  TU_CHRMAT* transpose,   /*< Transpose of \p matrix. */
  TU_GRAPH* graph,        /*< If not \c NULL and graphic, a graph represented by the matrix. */
  TU_GRAPH_EDGE* basis,   /*< If not \c NULL and graphic, a map from rows to basis edges. */
  TU_GRAPH_EDGE* cobasis, /*< If not \c NULL and graphic, a map from columns to cobasis edges. */
  TU_SUBMAT** psubmatrix  /*< If not \c NULL and not graphic, a minimal violating submatrix. */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_TDEC_H */
