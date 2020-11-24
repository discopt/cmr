#include "tdec.h"
#include "env_internal.h"

#include <assert.h>

typedef struct
{
  TDEC_NODE nextNode; /*< Next representative of same node towards root. */
} TDEC_NODE_DATA;

typedef struct
{
  int name;                 /*< Name of this edge. */
  TDEC_MEMBER member;       /*< Member this edge belongs to. */
  TDEC_NODE head;           /*< Head node of this edge. */
  TDEC_NODE tail;           /*< Tail node of this edge. */
  TDEC_EDGE prev;           /*< Next edge of this member. Must be a directed cycle if member is a polygon. */
  TDEC_EDGE next;           /*< Previous edge of this member. Must be a directed cycle if member is a polygon. */
  TDEC_MEMBER childMember;  /*< Child member linked to this edge. */
} TDEC_EDGE_DATA;

typedef struct
{
  int name;                  /*< Unique name. */
  TDEC_MEMBER nextMember;    /*< Next representative of same member towards root. */
  TDEC_MEMBER parentMember;  /*< Parent member of this member. Only valid if root representative. */
  int numEdges;              /*< Number of edges. Only valid if root representative. */
  TU_TDEC_MEMBER_TYPE type;     /*< Type of member. Only valid if root representative. */
  TDEC_EDGE markerToParent;  /*< Parent marker edge. Only valid if root representative. */
  TDEC_EDGE markerOfParent;  /*< Child marker of parent to which this member is linked. Only valid if root representative. */
} TDEC_MEMBER_DATA;

typedef struct
{
  TDEC_EDGE edge;
} TDEC_ROW_DATA;

typedef struct
{
  TDEC_EDGE edge;
} TDEC_COLUMN_DATA;

struct _TU_TDEC
{
  int memMembers;                 /*< Allocated memory for members. */
  int numMembers;                 /*< Number of members. */
  TDEC_MEMBER_DATA* members;      /*< Array of members. */
  TDEC_MEMBER firstFreeMember;    /*< First member in free list or -1. */
  int rootRow;                    /*< Unique row element in member 0. */

  int memEdges;                   /*< Allocated memory for edges. */
  int numEdges;                   /*< Number of used edges. */
  TDEC_EDGE_DATA* edges;          /*< Array of edges. */
  TDEC_EDGE firstFreeEdge;        /*< First edge in free list or -1. */

  int memNodes;                   /*< Allocated memory for nodes. */
  int numNodes;                   /*< Number of nodes. */
  TDEC_NODE_DATA* nodes;          /*< Array of nodes. */
  TDEC_NODE firstFreeNode;        /*< First node in free list or -1. */

  int memRows;                    /*< Allocated memory for \c rowEdges. */
  int numRows;                    /*< Number of rows. */
  TDEC_ROW_DATA* rowEdges;        /*< Maps each row to its edge (or NULL). */

  int memColumns;                 /*< Allocated memory for \c columnEdges. */
  int numColumns;                 /*< Number of columns. */
  TDEC_COLUMN_DATA* columnEdges;  /*< Maps each column to its edge (or NULL). */
};

struct _TU_TDEC_NEWCOLUMN
{
  bool isGraphic;
};

void TUtdecCreate(TU* tu, TU_TDEC** ptdec, int rootRow, int memEdges, int memNodes,
  int memMembers, int numRows, int numColumns)
{
  assert(tu);
  assert(ptdec);
  assert(!*ptdec);

  TUallocBlock(tu, ptdec);
  TU_TDEC* tdec = *ptdec;
  if (memMembers < 1)
    memMembers = 1;
  tdec->memMembers = memMembers;
  tdec->numMembers = 0;
  tdec->members = NULL;
  TUallocBlockArray(tu, &tdec->members, tdec->memMembers);
  for (int m = 1; m < memMembers; ++m)
    tdec->members[m].nextMember = m+1;
  if (memMembers > 1)
  {
    tdec->members[memMembers-1].nextMember = -1;
    tdec->firstFreeMember = 1;
  }
  else
    tdec->firstFreeMember = -1;
  tdec->members[0].name = rootRow;
  tdec->members[0].nextMember = -1;
  tdec->members[0].parentMember = -1;
  tdec->members[0].numEdges = 1;
  tdec->members[0].type = TDEC_MEMBER_TYPE_BOND;
  tdec->members[0].markerToParent = -1;
  tdec->members[0].markerOfParent = -1;
  tdec->rootRow = rootRow;

  if (memNodes < 2)
    memNodes = 2;
  tdec->memNodes = memNodes;
  tdec->nodes = NULL;
  TUallocBlockArray(tu, &tdec->nodes, memNodes);
  tdec->nodes[0].nextNode = 0;
  tdec->nodes[1].nextNode = 1;
  tdec->numNodes = 2;
  if (tdec->memNodes > 2)
  {
    for (int v = 2; v < memNodes; ++v)
      tdec->nodes[v].nextNode = v+1;
    tdec->nodes[memNodes-1].nextNode = -1;
    tdec->firstFreeNode = 2;
  }
  else
    tdec->firstFreeNode = -1;

  if (memEdges < 1)
    memEdges = 1;
  tdec->memEdges = memEdges;
  tdec->edges = NULL;
  TUallocBlockArray(tu, &tdec->edges, memEdges);  
  tdec->edges[0].name = rootRow;
  tdec->edges[0].member = -1;
  tdec->edges[0].head = 0;
  tdec->edges[0].tail = 1;
  tdec->edges[0].childMember = -1;
  tdec->edges[0].prev = 0;
  tdec->edges[0].next = 0;
  if (memEdges > 1)
  {
    for (int e = 1; e < memEdges; ++e)
      tdec->edges[e].next = e+1;
    tdec->edges[memEdges-1].next = -1;
    tdec->firstFreeEdge = 1;
  }
  else
    tdec->firstFreeEdge = -1;

  tdec->numRows = numRows > rootRow ? numRows : rootRow + 1;
  TUallocBlockArray(tu, &tdec->rowEdges, tdec->numRows);
  for (int r = 0; r < tdec->numRows; ++r)
    tdec->rowEdges[r].edge = -1;
  tdec->rowEdges[rootRow].edge = 0;

  tdec->numColumns = numColumns > 0 ? numColumns : 1;
  TUallocBlockArray(tu, &tdec->columnEdges, tdec->numColumns);
  for (int c = 0; c < tdec->numColumns; ++c)
    tdec->columnEdges[c].edge = -1;
}

void TUtdecFree(TU* tu, TU_TDEC** ptdec)
{
  assert(ptdec);
  assert(*ptdec);

  TU_TDEC* tdec = *ptdec;
  TUfreeBlockArray(tu, &tdec->members);
  TUfreeBlockArray(tu, &tdec->edges);
  TUfreeBlockArray(tu, &tdec->nodes);
  TUfreeBlockArray(tu, &tdec->rowEdges);
  TUfreeBlockArray(tu, &tdec->columnEdges);
  TUfreeBlock(tu, ptdec);
}


void TUtdecnewcolumnCreate(TU* tu, TU_TDEC_NEWCOLUMN** pnewcolumn)
{
  assert(tu);
  TUallocBlock(tu, pnewcolumn);
}

void TUtdecnewcolumnFree(TU* tu, TU_TDEC_NEWCOLUMN** pnewcolumn)
{
  assert(tu);
  TUfreeBlock(tu, pnewcolumn);
}

void TUtdecAddColumnPrepare(TU* tu, TU_TDEC* tdec, TU_TDEC_NEWCOLUMN* newcol, int* entryRows,
  int numEntries)
{
  assert(tu);
  assert(tdec);
  assert(newcol);
  assert(entryRows);
  assert(numEntries >= 1);

  newcol->isGraphic = true;
}

void TUtdecAddColumnApply(TU* tu, TU_TDEC* tdec, TU_TDEC_NEWCOLUMN* newcol)
{
  assert(tu);
  assert(tdec);
  assert(newcol);
}

bool testGraphicnessBixbyWagner(TU* tu, TU_CHRMAT* matrix, TU_CHRMAT* transpose,
  TU_GRAPH* graph, TU_GRAPH_EDGE* basis, TU_GRAPH_EDGE* cobasis, TU_SUBMAT** psubmatrix)
{
  assert(tu);
  assert(matrix);
  assert(transpose);
  assert(!graph || (TUgraphNumNodes(graph) == 0 && TUgraphNumEdges(graph) == 0));
  assert(!psubmatrix || !*psubmatrix);
  assert(!basis || graph);
  assert(!cobasis || graph);

  bool isGraphic = true;
  printf("testGraphicnessBixbyWagner called for a 1-connected %dx%d matrix.\n",
    matrix->numRows, matrix->numColumns);
  TUchrmatPrintDense(stdout, (TU_CHRMAT*) matrix, ' ', true);
  
  if (matrix->numNonzeros == 0)
  {
    if (graph)
    {
      /* Construct a path with numRows edges and with numColumns loops at 0. */

      TU_GRAPH_NODE s = TUgraphAddNode(tu, graph);
      for (int c = 0; c < matrix->numColumns; ++c)
      {
        TU_GRAPH_EDGE e = TUgraphAddEdge(tu, graph, s, s);
        if (cobasis)
          *cobasis++ = e;
      }
      for (int r = 0; r < matrix->numRows; ++r)
      {
        TU_GRAPH_NODE t = TUgraphAddNode(tu, graph);
        TU_GRAPH_EDGE e = TUgraphAddEdge(tu, graph, s, t);
        if (basis)
          *basis++ = e;
        s = t;
      }

      printf("Constructed graph with %d nodes and %d edges.\n", TUgraphNumNodes(graph),
        TUgraphNumEdges(graph));
    }
    return true;
  }

  int rootRow = transpose->entryColumns[0];
  printf("  Root row is %d.\n", rootRow);
  TU_TDEC* tdec = NULL;
  TUtdecCreate(tu, &tdec, rootRow, 0, 0, 0, 0, 0); /* TODO: avoid reallocations. */

  /* Process each column. */
  TU_TDEC_NEWCOLUMN* newcol = NULL;
  TUtdecnewcolumnCreate(tu, &newcol);
  for (int column = 0; column < matrix->numColumns; ++column)
  {
    TUtdecAddColumnPrepare(tu, tdec, newcol,
      &transpose->entryColumns[transpose->rowStarts[column]],
      transpose->rowStarts[column+1] - transpose->rowStarts[column]);

    if (newcol->isGraphic)
    {
      TUtdecAddColumnApply(tu, tdec, newcol);
    }
    else
    {
      assert(!"Not implemented");
    }
  }
  TUtdecnewcolumnFree(tu, &newcol);

  TUtdecFree(tu, &tdec);

  return isGraphic;
}

