#ifndef TU_TDEC_H
#define TU_TDEC_H

#include <tu/graph.h>
#include <tu/matrix.h>

#include "one_sum.h"

typedef enum
{
  TDEC_MEMBER_TYPE_BOND = 1,
  TDEC_MEMBER_TYPE_POLYGON = 2,
  TDEC_MEMBER_TYPE_PRIME = 3
} TU_TDEC_MEMBER_TYPE;

typedef int TDEC_EDGE;
typedef int TDEC_NODE;
typedef int TDEC_MEMBER;

struct _TU_TDEC;
typedef struct _TU_TDEC TU_TDEC;
typedef struct _TU_TDEC_NEWCOLUMN TU_TDEC_NEWCOLUMN;

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

void TUtdecFree(
  TU* tu,         /*< TU environment. */
  TU_TDEC** ptdec /*< Pointer to t-decomposition. .*/
);

void TUtdecnewcolumnCreate(
  TU* tu,                         /*< TU environment. */
  TU_TDEC_NEWCOLUMN** pnewcolumn  /*< new-column structure. */
);

void TUtdecnewcolumnFree(
  TU* tu,                         /*< TU environment. */
  TU_TDEC_NEWCOLUMN** pnewcolumn  /*< new-column structure. */
);

void TUtdecAddColumnPrepare(
  TU* tu,                     /*< TU environment. */
  TU_TDEC* tdec,              /*< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcol,  /*< new-column structure. */
  int* entryRows,             /*< Array of rows with 1-entry in this column. */
  int numEntries              /*< Number of 1-entries in this column. */
);

void TUtdecAddColumnApply(
  TU* tu,                   /*< TU environment. */
  TU_TDEC* tdec,            /*< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcol /*< new-column structure. */
);

bool testGraphicnessTDecomposition(
  TU* tu,                 /*< TU environment. */
  TU_CHRMAT* matrix,      /*< 1-connected matrix to be test. */
  TU_CHRMAT* transpose,   /*< Transpose of \p matrix. */
  TU_GRAPH* graph,        /*< If not \c NULL and graphic, a graph represented by the matrix. */
  TU_GRAPH_EDGE* basis,   /*< If not \c NULL and graphic, a map from rows to basis edges. */
  TU_GRAPH_EDGE* cobasis, /*< If not \c NULL and graphic, a map from columns to cobasis edges. */
  TU_SUBMAT** psubmatrix  /*< If not \c NULL and not graphic, a minimal violating submatrix. */
);

#endif /* TU_TDEC_H */
