#include "bixby_wagner.h"

#include <assert.h>

typedef enum
{
  TDEC_MEMBER_TYPE_BOND = 1,
  TDEC_MEMBER_TYPE_POLYGON = 2,
  TDEC_MEMBER_TYPE_PRIME = 3
} TDEC_MEMBER_TYPE;

struct _TDEC_EDGE;

typedef struct _TDEC_MEMBER
{
  int name;                           /*< Unique name. */
  struct _TDEC_MEMBER* nextMember;    /*< Next representative of same member towards root. */
  struct _TDEC_MEMBER* parentMember;  /*< Parent member of this member. Only valid if root representative. */
  int numEdges;                       /*< Number of edges. Only valid if root representative. */
  TDEC_MEMBER_TYPE type;              /*< Type of member. Only valid if root representative. */
  struct _TDEC_EDGE* markerToParent;  /*< Parent marker edge. Only valid if root representative. */
  struct _TDEC_EDGE* markerOfParent;  /*< Child marker of parent to which this member is linked. Only valid if root representative. */
} TDEC_MEMBER;

typedef struct _TDEC_NODE
{
  int name;                     /*< Name of this node. */
  struct _TDEC_NODE * nextNode; /*< Next representative of same node towards root. */
} TDEC_NODE;

typedef struct _TDEC_EDGE
{
  int name;                 /*< Name of this edge. */
  TDEC_MEMBER* member;      /*< Member this edge belongs to. */
  TDEC_NODE* head;          /*< Head node of this edge. */
  TDEC_NODE* tail;          /*< Tail node of this edge. */
  struct _TDEC_EDGE* prev;  /*< Next node of this member. Must be a directed cycle if member is a polygon. */
  struct _TDEC_EDGE* next;  /*< Previous node of this member. Must be a directed cycle if member is a polygon. */
  TDEC_MEMBER* childMember; /*< Child member linked to this edge. */
  // Maybe: boolean flag for whether this edge belongs to P.
} TDEC_EDGE;

typedef struct
{
  int memMembers;           /*< Allocated space for proper members. */
  int numMembers;           /*< Number of proper members. */
  TDEC_MEMBER* members;     /*< Array of proper members. */
  int rootRow;              /*< Unique row element in the root member. */
  TDEC_MEMBER rootMember;  /*< Root member (corresponding to a unit column). */
  int memEdges;             /*< Allocated space for edges. */
  TDEC_EDGE* edges;         /*< Array of edges. */
  int memNodes;             /*< Allocated space for nodes. */
  TDEC_NODE* nodes;         /*< Array of nodes. */
} TDEC;

static void createTDecomposition(TU* tu, TDEC** ptdec, int rootRow, int memEdges, int memNodes)
{
  assert(tu);
  assert(ptdec);
  assert(!*ptdec);

  TUallocBlock(tu, ptdec);
  TDEC* tdec = *ptdec;
  tdec->memMembers = 32;
  tdec->numMembers = 0;
  tdec->members = NULL;
  TUallocBlockArray(tu, &tdec->members, tdec->memMembers);
  tdec->rootRow = rootRow;
  tdec->rootMember.name = rootRow;
  tdec->rootMember.nextMember = NULL;
  tdec->rootMember.parentMember = NULL;
  tdec->rootMember.numEdges = 1;
  tdec->rootMember.type = TDEC_MEMBER_TYPE_BOND;
  tdec->rootMember.markerToParent = NULL;
  tdec->rootMember.markerOfParent = NULL;
  tdec->memEdges = memEdges;
  tdec->edges = NULL;
  TUallocBlockArray(tu, &tdec->edges, memEdges);
  tdec->memNodes = memNodes;
  tdec->nodes = NULL;
  TUallocBlockArray(tu, &tdec->nodes, memNodes);
}

static void freeTDecomposition(TU* tu, TDEC** ptdec)
{
  assert(ptdec);
  assert(*ptdec);

  TDEC* tdec = *ptdec;
  TUfreeBlockArray(tu, &tdec->members);
  TUfreeBlockArray(tu, &tdec->edges);
  TUfreeBlockArray(tu, &tdec->nodes);
  TUfreeBlock(tu, ptdec);
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
  TDEC* tdec = NULL;
  createTDecomposition(tu, &tdec, rootRow, 2*matrix->numRows + matrix->numColumns,
    2*matrix->numRows);

  /* Process each column. */
  for (int column = 0; column < matrix->numColumns; ++column)
  {
    
  }

  freeTDecomposition(tu, &tdec);

  return isGraphic;
}

