#define TU_DEBUG_TDEC /* Uncomment to enable debugging of t-decompositions. */

#include <tu/tdec.h>
#include "env_internal.h"

#include <assert.h>
#include <limits.h>

typedef struct
{
  TU_TDEC_NODE nextNode; /**< \brief Next representative of same node towards root, or -1 if root. */
} TU_TDEC_NODE_DATA;

typedef struct
{
  int name;                   /**< \brief Name of this edge. */
  TU_TDEC_MEMBER member;      /**< \brief Member this edge belongs to or -1 if in free list. */
  TU_TDEC_NODE head;          /**< \brief Head node of this edge. */
  TU_TDEC_NODE tail;          /**< \brief Tail node of this edge. */
  TU_TDEC_EDGE prev;          /**< \brief Next edge of this member. Must be a directed cycle if member is a polygon. */
  TU_TDEC_EDGE next;          /**< \brief Previous edge of this member. Must be a directed cycle if member is a polygon. */
  TU_TDEC_MEMBER childMember; /**< \brief Child member linked to this edge, or -1. */
} TU_TDEC_EDGE_DATA;

typedef struct
{
  TU_TDEC_MEMBER_TYPE type;     /**< \brief Type of member. Only valid if root representative. */
  TU_TDEC_MEMBER nextMember;    /**< \brief Next representative of same member towards root, or -1 if root. */
  TU_TDEC_MEMBER parentMember;  /**< \brief Parent member of this member. Only valid if root representative. */
  int numEdges;                 /**< \brief Number of edges. Only valid if root representative. */
  TU_TDEC_EDGE markerToParent;  /**< \brief Parent marker edge. Only valid if root representative. */
  TU_TDEC_EDGE markerOfParent;  /**< \brief Child marker of parent to which this member is linked. Only valid if root representative. */
  TU_TDEC_EDGE firstEdge;       /**< \brief First edge in doubly-linked edge list of this member. */
} TU_TDEC_MEMBER_DATA;

typedef struct
{
  TU_TDEC_EDGE edge;  /**< \brief Edge or -1. */
} TU_TDEC_ROW_DATA;

typedef struct
{
  TU_TDEC_EDGE edge;  /**< \brief Edge or -1. */
} TU_TDEC_COLUMN_DATA;

struct _TU_TDEC
{
  int memMembers;                   /**< \brief Allocated memory for members. */
  int numMembers;                   /**< \brief Number of members. */
  TU_TDEC_MEMBER_DATA* members;     /**< \brief Array of members. */
  TU_TDEC_MEMBER firstFreeMember;   /**< \brief First member in free list or -1. */
  int rootRow;                      /**< \brief Unique row element in member 0. */

  int memEdges;                     /**< \brief Allocated memory for edges. */
  int numEdges;                     /**< \brief Number of used edges. */
  TU_TDEC_EDGE_DATA* edges;         /**< \brief Array of edges. */
  TU_TDEC_EDGE firstFreeEdge;       /**< \brief First edge in free list or -1. */

  int memNodes;                     /**< \brief Allocated memory for nodes. */
  int numNodes;                     /**< \brief Number of nodes. */
  TU_TDEC_NODE_DATA* nodes;         /**< \brief Array of nodes. */
  TU_TDEC_NODE firstFreeNode;       /**< \brief First node in free list or -1. */

  int memRows;                      /**< \brief Allocated memory for \c rowEdges. */
  int numRows;                      /**< \brief Number of rows. */
  TU_TDEC_ROW_DATA* rowEdges;       /**< \brief Maps each row to its edge. */

  int memColumns;                   /**< \brief Allocated memory for \c columnEdges. */
  int numColumns;                   /**< \brief Number of columns. */
  TU_TDEC_COLUMN_DATA* columnEdges; /**< \brief Maps each column to its edge. */

  int numMarkers;                   /**< \brief Number of marker edge pairs in t-decomposition. */
};

typedef enum
{
  TYPE_UNKNOWN = 0,
  TYPE_CLOSES_CYCLE = 1,
  TYPE_SHORTCUTS_PATH = 2,
  TYPE_EXTENDS_PATH = 3,
  TYPE_CONNECTS_PATHS = 4,
  TYPE_OTHER = 5
} Type;

/**
 * \brief Additional edge information specific to a path.
 */

typedef struct _ReducedEdge
{
  TU_TDEC_EDGE edge;          /**< \brief The edge in the t-decomposition. */
  struct _ReducedEdge* next;  /**< \brief Next edge of this reduced member, or NULL. */
} ReducedEdge;

/**
 * \brief Additional member information specfic to a given path.
 */

typedef struct _ReducedMember
{
  TU_TDEC_MEMBER member;                /**< \brief The member from the t-decomposition. */
  int depth;                            /**< \brief Depth of this member in the reduced t-decomposition. */
  Type type;                            /**< \brief Type of this member. */
  int numChildren;                      /**< \brief Number of children in the reduced t-decomposition. */
  struct _ReducedMember** children;     /**< \brief Children in the reduced t-decomposition. */
  ReducedEdge* firstReducedEdge;        /**< \brief First edge in linked list of edges of this reduced member. */
} ReducedMember;

struct _TU_TDEC_NEWCOLUMN
{
  bool remainsGraphic;                      /**< \brief Indicator whether adding this column maintains graphicness. */
  int memReducedMembers;                    /**< \brief Allocated memory for \c reducedMembers. */
  int numReducedMembers;                    /**< \brief Number of members in \c reducedMembers. */
  ReducedMember* reducedMembers;            /**< \brief Array of reduced members, sorted by increasing depth. */
  ReducedMember** membersToReducedMembers;  /**< \brief Array mapping members to members of the reduced t-decomposition. */

  ReducedEdge* reducedEdgeStorage;          /**< \brief Storage for edge lists of reduced members. */
  int memReducedEdgeStorage;                /**< \brief Allocated memory for \c reducedEdgeStorage. */
  int usedReducedEdgeStorage;               /**< \brief Number of stored edges in \c reducedEdgeStorage. */

  ReducedMember** childrenStorage;          /**< \brief Storage for members' arrays of children in reduced t-decomposition. */
  int usedChildrenStorage;                  /**< \brief Number of stored children in \c childrenStorage. */
  int memChildrenStorage;                   /**< \brief Allocated memory for \c childrenStorage. */

  TU_TDEC_NODE terminalNode1;               /**< \brief First terminal node of path. */
  TU_TDEC_NODE terminalNode2;               /**< \brief Second terminal node of path. */
  TU_TDEC_MEMBER terminalMember1;           /**< \brief First terminal member of path. */
  TU_TDEC_MEMBER terminalMember2;           /**< \brief Second terminal member of path. */
};

int compareMemberDepths(const void* a, const void* b)
{
  const ReducedMember* first = a;
  const ReducedMember* second = b;
  /* Negative depths are moved to the end. */
  if (first->depth <= 0)
    return +1;
  if (second->depth <= 0)
    return -1;
  return first->depth - second->depth;
}

static TU_TDEC_MEMBER findMember(TU_TDEC* tdec, TU_TDEC_MEMBER start)
{
  TU_TDEC_MEMBER current = start;
  TU_TDEC_MEMBER next;
  while ((next = tdec->members[current].nextMember) >= 0)
    current = next;
  TU_TDEC_MEMBER root = current;
  current = start;
  while ((next = tdec->members[current].nextMember) >= 0)
  {
    if (next != root)
      tdec->members[current].nextMember = root;
    current = next;
  }
  return root;
}

static inline
TU_TDEC_MEMBER findMemberParent(TU_TDEC* tdec, TU_TDEC_MEMBER member)
{
  TU_TDEC_MEMBER someParent = tdec->members[member].parentMember;
  if (someParent >= 0)
    return findMember(tdec, someParent);
  else
    return -1;
}

static inline
TU_TDEC_MEMBER findEdgeMember(TU_TDEC* tdec, TU_TDEC_EDGE edge)
{
  return findMember(tdec, tdec->edges[edge].member);
}

static TU_TDEC_NODE findNode(TU_TDEC* tdec, TU_TDEC_NODE start)
{
  TU_TDEC_NODE current = start;
  TU_TDEC_NODE next;
#if defined(TU_DEBUG_TDEC)
  printf("        findNode(%d)\n", start);
#endif /* TU_DEBUG_TDEC */
  while ((next = tdec->nodes[current].nextNode) >= 0)
    current = next;
  TU_TDEC_NODE root = current;
  current = start;
  while ((next = tdec->nodes[current].nextNode) >= 0)
  {
    if (next != root)
      tdec->nodes[current].nextNode = root;
    current = next;
  }
  return root;
}

static inline
TU_TDEC_NODE findEdgeHead(TU_TDEC* tdec, TU_TDEC_EDGE edge)
{
  return findNode(tdec, tdec->edges[edge].head);
}

static inline
TU_TDEC_NODE findEdgeTail(TU_TDEC* tdec, TU_TDEC_EDGE edge)
{
  return findNode(tdec, tdec->edges[edge].tail);
}

static TU_TDEC_NODE createNode(
  TU* tu,       /**< \ref TU environment . */
  TU_TDEC* tdec /**< t-decomposition. */
)
{
  assert(tu);
  assert(tdec);

  TU_TDEC_NODE node = tdec->firstFreeNode;
  if (node >= 0)
  {
#if defined(TU_DEBUG_TDEC)
    printf("        createNode returns free node %d.\n", node);
#endif /* TU_DEBUG_TDEC */
    tdec->firstFreeNode = tdec->nodes[node].nextNode;
  }
  else /* No member in free list, so we enlarge the array. */
  {
    int newSize = 2 * tdec->memNodes + 16;
    TUreallocBlockArray(tu, &tdec->nodes, newSize);
    for (int v = tdec->memNodes + 1; v < newSize; ++v)
      tdec->nodes[v].nextNode = v+1;
    tdec->nodes[newSize-1].nextNode = -1;
    tdec->firstFreeNode = tdec->memNodes + 1;
    node = tdec->memNodes;
    tdec->memNodes = newSize;
#if defined(TU_DEBUG_TDEC)
    printf("        createNode enlarges node array to %d and returns node %d.\n", newSize, node);
#endif /* TU_DEBUG_TDEC */
  }
  tdec->nodes[node].nextNode = -1;
  tdec->numNodes++;

  return node;
}

static void setRowEdge(
  TU* tu,           /**< \ref TU environment. */
  TU_TDEC* tdec,    /**< t-decomposition. */
  int row ,         /**< Row (index). */
  TU_TDEC_EDGE edge /**< Edge to be assigned to \p row. */
)
{
  assert(tu);
  assert(tdec);
  assert(row >= 0);
  assert(edge >= 0);

  if (row >= tdec->memRows)
  {
    int newSize = 2*tdec->memRows + 16;
    TUreallocBlockArray(tu, &tdec->rowEdges, newSize);
    for (int c = tdec->memRows; c < newSize; ++c)
      tdec->rowEdges[c].edge = -1;
    tdec->memRows = newSize;
  }

  assert(tdec->rowEdges[row].edge == -1);
  tdec->rowEdges[row].edge = edge;
  if (row >= tdec->numRows)
    tdec->numRows = row + 1;
}

static void setColumnEdge(
  TU* tu,           /**< \ref TU environment. */
  TU_TDEC* tdec,    /**< t-decomposition. */
  int column,       /**< Column (index). */
  TU_TDEC_EDGE edge /**< Edge to be assigned to \p column. */
)
{
  assert(tu);
  assert(tdec);
  assert(column >= 0);
  assert(edge >= 0);

  if (column >= tdec->memColumns)
  {
    int newSize = 2*tdec->memColumns + 16;
    TUreallocBlockArray(tu, &tdec->columnEdges, newSize);
    for (int c = tdec->memColumns; c < newSize; ++c)
      tdec->columnEdges[c].edge = -1;
    tdec->memColumns = newSize;
  }

  assert(tdec->columnEdges[column].edge == -1);
  tdec->columnEdges[column].edge = edge;
  if (column >= tdec->numColumns)
    tdec->numColumns = column + 1;
}

/**
 * \brief Adds \p edge to the edge list of \p member.
 */

static
TU_ERROR addEdgeToMember(
  TU* tu,               /**< \ref TU environment. */
  TU_TDEC* tdec,        /**< t-decomposition. */
  TU_TDEC_EDGE edge,    /**< Edge to be added. */ 
  TU_TDEC_MEMBER member /**< Member. */
)
{
  assert(tu);
  assert(tdec);
  assert(edge >= 0);
  assert(member >= 0);
  assert(member < tdec->numMembers);

  tdec->edges[edge].next = tdec->members[member].firstEdge;
  if (tdec->members[member].firstEdge >= 0)
  {
    assert(tdec->edges[tdec->members[member].firstEdge].prev == -1);
    tdec->edges[tdec->members[member].firstEdge].prev = edge;
  }
  tdec->edges[edge].prev = -1;
  tdec->members[member].firstEdge = edge;

  return TU_OKAY;
}

/**
 * \brief Creates a new edge.
 */

static
TU_ERROR createEdge(
  TU* tu,                 /**< \ref TU environment. */
  TU_TDEC* tdec,          /**< t-decomposition. */
  TU_TDEC_MEMBER member,  /**< Member this edge belongs to. */
  TU_TDEC_EDGE* pedge     /**< Pointer for storing the new edge. */
)
{
  assert(tu);
  assert(tdec);
  assert(pedge);

  TU_TDEC_EDGE edge = tdec->firstFreeEdge;
  if (edge >= 0)
  {
#if defined(TU_DEBUG_TDEC)
    printf("        createEdge returns free edge %d.\n", edge);
#endif /* TU_DEBUG_TDEC */
    tdec->firstFreeEdge = tdec->edges[edge].next;
  }
  else /* No edge in free list, so we enlarge the array. */
  {
    int newSize = 2 * tdec->memEdges + 16;
    TU_CALL( TUreallocBlockArray(tu, &tdec->edges, newSize) );
    for (int e = tdec->memEdges + 1; e < newSize; ++e)
    {
      tdec->edges[e].next = e+1;
      tdec->edges[e].member = -1;
    }
    tdec->edges[newSize-1].next = -1;
    tdec->firstFreeEdge = tdec->memEdges + 1;
    edge = tdec->memEdges;
    tdec->memEdges = newSize;
#if defined(TU_DEBUG_TDEC)
    printf("        createEdge enlarges edge array to %d and returns edge %d.\n", newSize, edge);
#endif /* TU_DEBUG_TDEC */
  }

  tdec->edges[edge].member = member;
  tdec->numEdges++;

  if (member >= 0)
  {
#if defined(TU_DEBUG_TDEC)
    printf("        createEdge inserts into member's edge list.\n");
#endif /* TU_DEBUG_TDEC */
    TU_CALL( addEdgeToMember(tu, tdec, edge, member) );
  }

  *pedge = edge;
  
  return TU_OKAY;
}

static
TU_ERROR createRowEdge(
  TU* tu,                 /**< \ref TU environment. */
  TU_TDEC* tdec,          /**< t-decomposition. */
  TU_TDEC_EDGE* pedge,    /**< Pointer for storing the new edge. */
  TU_TDEC_MEMBER member,  /**< Member this edge belongs to. */
  TU_TDEC_NODE head,      /**< Head node of this edge. */
  TU_TDEC_NODE tail,      /**< Tail node of this edge. */
  int row                 /**< Row (index) this edge corresponds to. */
)
{
  assert(tu);
  assert(tdec);
  assert(pedge);

  TU_CALL( createEdge(tu, tdec, member, pedge) );
  TU_TDEC_EDGE edge = *pedge;
  TU_TDEC_EDGE_DATA* data = &tdec->edges[edge];
  data->head = head;
  data->tail = tail;
  data->childMember = -1;
  data->name = row;
  setRowEdge(tu, tdec, row, edge);
#ifndef NDEBUG
  data->next = INT_MIN;
  data->prev = INT_MIN;
#endif /* !NDEBUG */

#if defined(TU_DEBUG_TDEC)
  printf("        Created row edge {%d,%d} of member %d for row %d.\n", head, tail, member, row);
#endif /* TU_DEBUG_TDEC */

  return TU_OKAY;
}

static
TU_ERROR createColumnEdge(
  TU* tu,                 /*< TU environment. */
  TU_TDEC* tdec,          /*< t-decomposition. */
  TU_TDEC_EDGE* pedge,    /**< Pointer for storing the new edge. */
  TU_TDEC_MEMBER member,  /*< Member this edge belongs to. */
  TU_TDEC_NODE head,      /*< Head node of this edge. */
  TU_TDEC_NODE tail,      /*< Tail node of this edge. */
  int column              /*< Column (index) this edge corresponds to. */
)
{
  assert(tu);
  assert(tdec);
  assert(pedge);

  TU_CALL( createEdge(tu, tdec, member, pedge) );
  TU_TDEC_EDGE edge = *pedge;
  TU_TDEC_EDGE_DATA* data = &tdec->edges[edge];
  data->head = head;
  data->tail = tail;
  data->childMember = -1;
  data->name = -1-column;
  setColumnEdge(tu, tdec, column, edge);
#ifndef NDEBUG
  data->next = INT_MIN;
  data->prev = INT_MIN;
#endif /* !NDEBUG */

#if defined(TU_DEBUG_TDEC)
  printf("        Created column edge {%d,%d} of member %d for column %d.\n", head, tail, member,
    column);
#endif /* TU_DEBUG_TDEC */

  return edge;
}

static
TU_ERROR createMarkerEdge(
  TU* tu,                 /**< \ref TU environment. */
  TU_TDEC* tdec,          /**< t-decomposition. */
  TU_TDEC_EDGE* pedge,    /**< Pointer for storing the new edge. */
  TU_TDEC_MEMBER member,  /**< Member this edge belongs to. */
  TU_TDEC_NODE head,      /**< Head node of this edge. */
  TU_TDEC_NODE tail,      /**< Tail node of this edge. */
  bool isParent           /**< Whether this is the parent marker edge. */
)
{
  assert(tu);
  assert(tdec);
  assert(pedge);

  TU_CALL( createEdge(tu, tdec, member, pedge) );
  TU_TDEC_EDGE edge = *pedge;
  TU_TDEC_EDGE_DATA* data = &tdec->edges[edge];
  data->head = head;
  data->tail = tail;
  data->childMember = -1;
  if (isParent)
    data->name = INT_MAX - tdec->numMarkers;
  else
    data->name = -(INT_MAX - tdec->numMarkers);
#ifndef NDEBUG
  data->next = INT_MIN;
  data->prev = INT_MIN;
#endif /* !NDEBUG */

#if defined(TU_DEBUG_TDEC)
  printf("        Created %s marker edge {%d,%d} of member %d.\n", isParent ? "parent" : "child",
    head, tail, member);
#endif /* TU_DEBUG_TDEC */

  return TU_OKAY;
}

static
TU_ERROR createMember(
  TU* tu,                   /**< \ref TU environment . */
  TU_TDEC* tdec,            /**< t-decomposition. */
  TU_TDEC_MEMBER_TYPE type, /**< Type of member. */
  TU_TDEC_MEMBER* pmember   /**< Created member. */
)
{
  assert(tu);
  assert(tdec);

  TU_TDEC_MEMBER member = tdec->firstFreeMember;
  if (member >= 0)
  {
#if defined(TU_DEBUG_TDEC)
    printf("        createMember returns free member %d.\n", member);
#endif /* TU_DEBUG_TDEC */
    tdec->firstFreeMember = tdec->members[member].nextMember;
  }
  else /* No member in free list, so we enlarge the array. */
  {
    int newSize = 2 * tdec->memMembers + 16;
    TU_CALL( TUreallocBlockArray(tu, &tdec->members, newSize) );
    for (int m = tdec->memMembers + 1; m < newSize; ++m)
      tdec->members[m].nextMember = m+1;
    tdec->members[newSize-1].nextMember = -1;
    tdec->firstFreeMember = tdec->memMembers + 1;
    member = tdec->memMembers;
    tdec->memMembers = newSize;
#if defined(TU_DEBUG_TDEC)
    printf("        createMember enlarges member array to %d and returns member %d.\n", newSize,
      member);
#endif /* TU_DEBUG_TDEC */
  }
  TU_TDEC_MEMBER_DATA* data = &tdec->members[member];
  data->markerOfParent = INT_MIN;
  data->markerToParent = INT_MIN;
  data->firstEdge = -1;
  data->nextMember = -1;
  data->numEdges = 0;
  data->parentMember = INT_MIN;
  data->type = type;

  tdec->numMembers++;

  *pmember = member;

  return TU_OKAY;
}

TU_ERROR TUtdecCreate(TU* tu, TU_TDEC** ptdec, int rootRow, int memEdges, int memNodes,
  int memMembers, int numRows, int numColumns)
{
  assert(tu);
  assert(ptdec);
  assert(!*ptdec);

  TU_CALL( TUallocBlock(tu, ptdec) );
  TU_TDEC* tdec = *ptdec;
  if (memMembers < 1)
    memMembers = 1;
  tdec->memMembers = memMembers;
  tdec->numMembers = 1;
  tdec->members = NULL;
  TU_CALL( TUallocBlockArray(tu, &tdec->members, tdec->memMembers) );
  for (int m = 1; m < memMembers; ++m)
    tdec->members[m].nextMember = m+1;
  if (memMembers > 1)
  {
    tdec->members[memMembers-1].nextMember = -1;
    tdec->firstFreeMember = 1;
  }
  else
    tdec->firstFreeMember = -1;
  tdec->members[0].nextMember = -1;
  tdec->members[0].parentMember = -1;
  tdec->members[0].numEdges = 2;
  tdec->members[0].type = TDEC_MEMBER_TYPE_BOND;
  tdec->members[0].markerToParent = -1;
  tdec->members[0].markerOfParent = -1;
  tdec->rootRow = rootRow;

  if (memNodes < 2)
    memNodes = 2;
  tdec->memNodes = memNodes;
  tdec->nodes = NULL;
  TU_CALL( TUallocBlockArray(tu, &tdec->nodes, memNodes) );
  tdec->nodes[0].nextNode = -1;
  tdec->nodes[1].nextNode = -1;
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

  if (memEdges < 2)
    memEdges = 2;
  tdec->memEdges = memEdges;
  tdec->edges = NULL;
  TU_CALL( TUallocBlockArray(tu, &tdec->edges, memEdges) );
  tdec->numEdges = 2;

  /* First edge is co-tree edge corresponding to artificial column. */
  tdec->edges[0].name = INT_MIN;
  tdec->numMarkers = 0;
  tdec->edges[0].member = 0;
  tdec->edges[0].head = 0;
  tdec->edges[0].tail = 1;
  tdec->edges[0].childMember = -1;
  tdec->edges[0].prev = 1;
  tdec->edges[0].next = 1;

  /* Second edge is tree edge corresponding to \c rootRow. */
  tdec->edges[1].name = rootRow;
  tdec->edges[1].member = 0;
  tdec->edges[1].head = 0;
  tdec->edges[1].tail = 1;
  tdec->edges[1].childMember = -1;
  tdec->edges[1].prev = 0;
  tdec->edges[1].next = 0;
  if (memEdges > 2)
  {
    for (int e = 2; e < memEdges; ++e)
    {
      tdec->edges[e].next = e+1;
      tdec->edges[e].member = -1;
    }
    tdec->edges[memEdges-1].next = -1;
    tdec->firstFreeEdge = 2;
  }
  else
    tdec->firstFreeEdge = -1;

  tdec->numRows = numRows > rootRow ? numRows : rootRow + 1;
  tdec->memRows = numRows;
  tdec->rowEdges = NULL;
  TU_CALL( TUallocBlockArray(tu, &tdec->rowEdges, tdec->numRows) );
  for (int r = 0; r < tdec->numRows; ++r)
    tdec->rowEdges[r].edge = -1;
  tdec->rowEdges[rootRow].edge = 1;

  tdec->numColumns = numColumns > 0 ? numColumns : 1;
  tdec->memColumns = tdec->numColumns;
  tdec->columnEdges = NULL;
  TU_CALL( TUallocBlockArray(tu, &tdec->columnEdges, tdec->numColumns) );
  for (int c = 0; c < tdec->numColumns; ++c)
    tdec->columnEdges[c].edge = -1;

  return TU_OKAY;
}

TU_ERROR TUtdecFree(TU* tu, TU_TDEC** ptdec)
{
  assert(ptdec);
  assert(*ptdec);

  TU_TDEC* tdec = *ptdec;
  TU_CALL( TUfreeBlockArray(tu, &tdec->members) );
  TU_CALL( TUfreeBlockArray(tu, &tdec->edges) );
  TU_CALL( TUfreeBlockArray(tu, &tdec->nodes) );
  TU_CALL( TUfreeBlockArray(tu, &tdec->rowEdges) );
  TU_CALL( TUfreeBlockArray(tu, &tdec->columnEdges) );
  TU_CALL( TUfreeBlock(tu, ptdec) );

  return TU_OKAY;
}

int TUtdecBasisSize(TU_TDEC* tdec)
{
  assert(tdec);

  return tdec->numRows;
}

int TUtdecCobasisSize(TU_TDEC* tdec)
{
  assert(tdec);

  return tdec->numColumns;
}

int TUtdecNumEdges(TU_TDEC* tdec)
{
  assert(tdec);

  return tdec->numEdges;
}

TU_ERROR TUtdecToGraph(TU* tu, TU_TDEC* tdec, TU_GRAPH* graph, bool merge, TU_GRAPH_EDGE* basis,
  TU_GRAPH_EDGE* cobasis, int* edgeElements)
{
  assert(tu);
  assert(tdec);
  assert(graph);

#if defined(TU_DEBUG_TDEC)
  printf("TUtdecToGraph for t-decomposition.\n");
#endif /* TU_DEBUG_TDEC */

  TU_CALL( TUgraphClear(tu, graph) );

  TU_GRAPH_EDGE* localEdgeElements = NULL;
  if (edgeElements)
    localEdgeElements = edgeElements;
  else if (basis || cobasis)
    TU_CALL( TUallocStackArray(tu, &localEdgeElements, tdec->numEdges) );
  TU_GRAPH_NODE* tdecNodesToGraphNodes = NULL;
  TU_CALL( TUallocStackArray(tu, &tdecNodesToGraphNodes, tdec->numNodes) );
  TU_GRAPH_EDGE* tdecEdgesToGraphEdges = NULL;
  TU_CALL( TUallocStackArray(tu, &tdecEdgesToGraphEdges, tdec->numEdges) );

  for (int v = 0; v < tdec->numNodes; ++v)
  {
    if (tdec->nodes[v].nextNode < 0)
    {
      TU_CALL( TUgraphAddNode(tu, graph, &tdecNodesToGraphNodes[v]) );
    }
    else
      tdecNodesToGraphNodes[v] = -1;
  }

  assert(tdec->edges[0].name == INT_MIN);
  for (int e = 1; e < tdec->numEdges; ++e)
  {
    if (tdec->edges[e].member >= 0)
    {
      TU_TDEC_NODE head = findEdgeHead(tdec, e);
      TU_TDEC_NODE tail = findEdgeTail(tdec, e);
      TU_TDEC_EDGE edge;
      TU_CALL( TUgraphAddEdge(tu, graph, tdecNodesToGraphNodes[head], tdecNodesToGraphNodes[tail],
        &edge) );
      tdecEdgesToGraphEdges[e] = edge;
      assert(edge < tdec->numEdges);
      if (localEdgeElements)
        localEdgeElements[edge] = tdec->edges[e].name;
    }
  }

  /* Merge respective parent and child edges. */

  if (merge)
  {
    for (int m = 1; m < tdec->numMembers; ++m)
    {
      if (tdec->members[m].type == TDEC_MEMBER_TYPE_INVALID)
        continue;

      TU_GRAPH_EDGE parent = tdecEdgesToGraphEdges[tdec->members[m].markerOfParent];
      TU_GRAPH_EDGE child = tdecEdgesToGraphEdges[tdec->members[m].markerToParent];
      TU_GRAPH_NODE parentU = TUgraphEdgeU(graph, parent);
      TU_GRAPH_NODE parentV = TUgraphEdgeV(graph, parent);
      TU_GRAPH_NODE childU = TUgraphEdgeU(graph, child);
      TU_GRAPH_NODE childV = TUgraphEdgeV(graph, child);

      TU_CALL( TUgraphMergeNodes(tu, graph, parentU, childU) );      
      TU_CALL( TUgraphDeleteNode(tu, graph, childU) );
      TU_CALL( TUgraphMergeNodes(tu, graph, parentV, childV) );
      TU_CALL( TUgraphDeleteNode(tu, graph, childV) );

      TU_CALL( TUgraphDeleteEdge(tu, graph, parent) );
      TU_CALL( TUgraphDeleteEdge(tu, graph, child) );
    }
  }

  /* Construct (co)basis. */

  if (basis || cobasis)
  {
    for (TU_GRAPH_ITER i = TUgraphEdgesFirst(graph); TUgraphEdgesValid(graph, i);
      i = TUgraphEdgesNext(graph, i))
    {
      TU_GRAPH_EDGE e = TUgraphEdgesEdge(graph, i);
      int element = localEdgeElements[e];
      if (element >= 0 && basis)
        basis[element] = e;
      else if (element < 0 && cobasis)
        cobasis[-1-element] = e;
    }
  }

  TU_CALL( TUfreeStackArray(tu, &tdecEdgesToGraphEdges) );
  TU_CALL( TUfreeStackArray(tu, &tdecNodesToGraphNodes) );
  if (localEdgeElements != edgeElements)
    TU_CALL( TUfreeStackArray(tu, &localEdgeElements) );

  return TU_OKAY;
}

TU_ERROR TUtdecnewcolumnCreate(TU* tu, TU_TDEC_NEWCOLUMN** pnewcolumn)
{
  assert(tu);

  TU_CALL( TUallocBlock(tu, pnewcolumn) );
  TU_TDEC_NEWCOLUMN* newcolumn = *pnewcolumn;
  newcolumn->remainsGraphic = true;
  newcolumn->memReducedMembers = 0;
  newcolumn->numReducedMembers = 0;
  newcolumn->reducedMembers = NULL;
  newcolumn->membersToReducedMembers = NULL;

  newcolumn->reducedEdgeStorage = NULL;
  newcolumn->memReducedEdgeStorage = 0;
  newcolumn->usedReducedEdgeStorage = 0;

  newcolumn->memChildrenStorage = 0;
  newcolumn->usedChildrenStorage = 0;
  newcolumn->childrenStorage = NULL;

  return TU_OKAY;
}

TU_ERROR TUtdecnewcolumnFree(TU* tu, TU_TDEC_NEWCOLUMN** pnewcolumn)
{
  assert(tu);
  assert(*pnewcolumn);

  TU_TDEC_NEWCOLUMN* newcolumn = *pnewcolumn;
  
  if (newcolumn->reducedMembers)
  {
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->reducedMembers) );
  }
  if (newcolumn->membersToReducedMembers)
  {
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->membersToReducedMembers) );
  }
  if (newcolumn->reducedEdgeStorage)
  {
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->reducedEdgeStorage) );
  }
  if (newcolumn->childrenStorage)
  {
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->childrenStorage) );
  }

  TU_CALL( TUfreeBlock(tu, pnewcolumn) );

  return TU_OKAY;
}

static
TU_ERROR initializeNewColumn(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn  /**< new column. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);

#ifndef NDEBUG
  newcolumn->terminalNode1 = INT_MIN;
  newcolumn->terminalNode2 = INT_MIN;
  newcolumn->terminalMember1 = INT_MIN;
  newcolumn->terminalMember2 = INT_MIN;
  newcolumn->numReducedMembers = 0;
#endif /* !NDEBUG */

  newcolumn->remainsGraphic = true;
  newcolumn->usedReducedEdgeStorage = 0;

  return TU_OKAY;
}

static
TU_ERROR findReducedDecomposition(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new column. */
  int* entryRows,               /**< Array of rows of new column's enries. */
  int numEntries                /**< Length of \p entryRows. */
)
{
  /* Enlarge members array. */
  if (newcolumn->memReducedMembers < tdec->numMembers)
  {
    newcolumn->memReducedMembers = tdec->memMembers;
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->reducedMembers, newcolumn->memReducedMembers) );
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->membersToReducedMembers,
      newcolumn->memReducedMembers) );
  }

  /* Identify all members on the path. For the induced sub-arborescence we also compute the
   * depths. After the computation, its root has depth pathRootDepth. */
#if defined(TU_DEBUG_TDEC)
  printf("    Finding reduced t-decomposition.\n");
#endif /* TU_DEBUG_TDEC */

  int* memberDepths = NULL;
  TU_CALL( TUallocStackArray(tu, &memberDepths, tdec->numMembers) );
  for (int m = 0; m < tdec->numMembers; ++m)
    memberDepths[m] = 0;
  TU_TDEC_MEMBER reducedRootMember = -1;
  int reducedRootDepth = 0;
  newcolumn->numReducedMembers = 0;
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    TU_TDEC_EDGE edge = (row < tdec->numRows) ? tdec->rowEdges[row].edge : -1;
    if (edge >= 0)
    {
#if defined(TU_DEBUG_TDEC)
      printf("      Edge %d exists.\n", edge);
#endif /* TU_DEBUG_TDEC */
      TU_TDEC_MEMBER member = findEdgeMember(tdec, edge);
      if (!reducedRootDepth)
      {
        /* The first member receives a #members, its parent #members-1, etc. */
        reducedRootMember = member;
        int depth = tdec->numMembers;
        reducedRootDepth = depth;
        while (member >= 0)
        {
#if defined(TU_DEBUG_TDEC)
          printf("        Member %d receives depth %d on initial path.\n", member, depth);
#endif /* TU_DEBUG_TDEC */
          assert(memberDepths[member] == 0);
          memberDepths[member] = depth;
          newcolumn->reducedMembers[newcolumn->numReducedMembers].member = member;
          newcolumn->reducedMembers[newcolumn->numReducedMembers].depth = depth;
          newcolumn->numReducedMembers++;
          --depth;
          member = findMemberParent(tdec, member);
        }
      }
      else
      {
        int count = 0;
        TU_TDEC_MEMBER m = member;
        while (memberDepths[m] == 0)
        {
          memberDepths[m] = 1;
          ++count;
          m = findMemberParent(tdec, m);
          assert(m >= 0);
        }
        for (m = member; count; --count)
        {
          memberDepths[m] += count;

#if defined(TU_DEBUG_TDEC)
          printf("        Member %d receives %d on subsequent path.\n", m, memberDepths[m]);
#endif /* TU_DEBUG_TDEC */
          newcolumn->reducedMembers[newcolumn->numReducedMembers].member = m;
          newcolumn->reducedMembers[newcolumn->numReducedMembers].depth = memberDepths[m];
          newcolumn->numReducedMembers++;
          m = findMemberParent(tdec, m);
        }
        if (memberDepths[m] < reducedRootDepth)
        {
          reducedRootMember = member;
          reducedRootDepth = memberDepths[m];
        }
      }
    }
    else
    {
#if defined(TU_DEBUG_TDEC)
      printf("      Edge %d is new.\n", edge);
#endif /* TU_DEBUG_TDEC */
    }
  }

#if defined(TU_DEBUG_TDEC)
  printf("      Root member is %d with temporary depth %d.\n", reducedRootMember, reducedRootDepth);
#endif /* TU_DEBUG_TDEC */

  for (int i = 0; i < newcolumn->numReducedMembers; ++i)
  {
#if defined(TU_DEBUG_TDEC)
    printf("        Shifting depth of %d: %d -> %d\n", newcolumn->reducedMembers[i].member,
      newcolumn->reducedMembers[i].depth,
      newcolumn->reducedMembers[i].depth + 1 - reducedRootDepth);
#endif /* TU_DEBUG_TDEC */
    newcolumn->reducedMembers[i].depth += 1 - reducedRootDepth;
  }
  
  qsort(newcolumn->reducedMembers, newcolumn->numReducedMembers, sizeof(ReducedMember),
    compareMemberDepths);

  for (int i = 0; i < newcolumn->numReducedMembers; ++i)
  {
    /* Remove members corresponding to non-positive depths since these are behind the root. */
    if (newcolumn->reducedMembers[i].depth <= 0)
      newcolumn->numReducedMembers = i;
#if defined(TU_DEBUG_TDEC)
    printf("        Member %d has depth %d.\n", newcolumn->reducedMembers[i].member,
      newcolumn->reducedMembers[i].depth);
#endif /* TU_DEBUG_TDEC */
  }
  assert(newcolumn->numReducedMembers > 0);
  assert(newcolumn->reducedMembers[0].member == reducedRootMember);
  
  /* We now create the mapping from members to reduced members. */
  for (int m = 0; m < tdec->numMembers; ++m)
    newcolumn->membersToReducedMembers[m] = NULL;
  for (int i = 0; i < newcolumn->numReducedMembers; ++i)
    newcolumn->membersToReducedMembers[newcolumn->reducedMembers[i].member] = &newcolumn->reducedMembers[i];

  TU_CALL( TUfreeStackArray(tu, &memberDepths) );

  return TU_OKAY;
}

static
TU_ERROR initializeReducedMemberEdgeLists(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new column. */
  int* entryRows,               /**< Array of rows of new column's enries. */
  int numEntries                /**< Length of \p entryRows. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("    Initializing edge lists for members of reduced t-decomposition.\n");
#endif /* TU_DEBUG_TDEC */

  /* (Re)allocate memory for edge lists. */
  assert(newcolumn->usedReducedEdgeStorage == 0);
  int requiredMemReducedEdgeStorage = numEntries;
  if (newcolumn->memReducedEdgeStorage < requiredMemReducedEdgeStorage)
  {
    newcolumn->memReducedEdgeStorage = 2 * requiredMemReducedEdgeStorage;
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->reducedEdgeStorage,
      newcolumn->memReducedEdgeStorage) );
  }

  /* Start with empty lists. */
  for (int i = 0; i < newcolumn->numReducedMembers; ++i)
    newcolumn->reducedMembers[i].firstReducedEdge = NULL;

  /* Fill edge lists. */
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    TU_TDEC_EDGE edge = (row < tdec->numRows) ? tdec->rowEdges[row].edge : -1;
    if (edge >= 0)
    {
      TU_TDEC_MEMBER member = tdec->edges[edge].member;
      assert(member >= 0);
      ReducedMember* reducedMember = newcolumn->membersToReducedMembers[member];
      newcolumn->reducedEdgeStorage[newcolumn->usedReducedEdgeStorage].next = reducedMember->firstReducedEdge;
      newcolumn->reducedEdgeStorage[newcolumn->usedReducedEdgeStorage].edge = edge;
      reducedMember->firstReducedEdge = &newcolumn->reducedEdgeStorage[newcolumn->usedReducedEdgeStorage];
      ++newcolumn->usedReducedEdgeStorage;

#if defined(TU_DEBUG_TDEC)
      printf("      Edge %d <%d> belongs to reduced member %ld which is member %d.\n", edge,
        tdec->edges[edge].name,
        (reducedMember - newcolumn->reducedMembers) / sizeof(ReducedMember),
        reducedMember->member);
#endif /* TU_DEBUG_TDEC */
    }
  }

  return TU_OKAY;
}

static
TU_ERROR computeReducedMemberChildren(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn  /**< new column. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("    Computing children of members of reduced t-decomposition.\n");
#endif /* TU_DEBUG_TDEC */

  /* Allocate memory for children of reduced members. */
  int requiredMemChildren = newcolumn->numReducedMembers;
  if (newcolumn->memChildrenStorage < requiredMemChildren)
  {
    newcolumn->memChildrenStorage = 2 * requiredMemChildren;
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->childrenStorage, newcolumn->memChildrenStorage) );
  }

  /* Initialize numChildren to zero for all reduced members. */
  newcolumn->usedChildrenStorage = 0;
  for (int m = 0; m < newcolumn->numReducedMembers; ++m)
    newcolumn->reducedMembers[m].numChildren = 0;

  /* Count children of each reduced member. */
  for (int m = 1; m < newcolumn->numReducedMembers; ++m)
  {
    ReducedMember* parentReducedMember = newcolumn->membersToReducedMembers[
      tdec->members[newcolumn->reducedMembers[m].member].parentMember];

#if defined(TU_DEBUG_TDEC)
    printf("Reduced member %d (= member %d) has parent %ld (= member %d).\n",
      m, newcolumn->reducedMembers[m].member,
      (parentReducedMember - newcolumn->reducedMembers) / sizeof(ReducedMember),
      tdec->members[newcolumn->reducedMembers[m].member].parentMember);
#endif /* TU_DEBUG_TDEC */
    parentReducedMember->numChildren++;
  }

  /* Set memory pointer of each reduced member. */
  for (int m = 0; m < newcolumn->numReducedMembers; ++m)
  {
    newcolumn->reducedMembers[m].children = &newcolumn->childrenStorage[newcolumn->usedChildrenStorage];
    newcolumn->usedChildrenStorage += newcolumn->reducedMembers[m].numChildren;
    newcolumn->reducedMembers[m].numChildren = 0;
  }

  /* Set children of each reduced member. */
  for (int m = 1; m < newcolumn->numReducedMembers; ++m)
  {
    ReducedMember* parentReducedMember = newcolumn->membersToReducedMembers[
      tdec->members[newcolumn->reducedMembers[m].member].parentMember];

#if defined(TU_DEBUG_TDEC)
    printf("Reduced member %ld (= member %d) has %d (= member %d) as child %d.\n",
      (parentReducedMember - newcolumn->reducedMembers) / sizeof(ReducedMember),
      tdec->members[newcolumn->reducedMembers[m].member].parentMember,
      m, newcolumn->reducedMembers[m].member, parentReducedMember->numChildren);
#endif /* TU_DEBUG_TDEC */
    parentReducedMember->children[parentReducedMember->numChildren] = &newcolumn->reducedMembers[m];
    parentReducedMember->numChildren++;
  }

  return TU_OKAY;
}

/**
 * \brief Count the number of children of a reduced member having certain types.
 */

static
TU_ERROR countChildrenTypes(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new column. */
  ReducedMember* reducedMember, /**< Reduced member. */
  int* pNumClosing,             /**< Number of children of type \ref TYPE_SHORTCUTS_PATH. */
  int* pNumShortcuts,           /**< Number of children of type \ref TYPE_SHORTCUTS_PATH. */
  int* pNumExtending,           /**< Number of children of type \ref TYPE_EXTENDS_PATH. */
  int* pNumConnecting           /**< Number of children of type \ref TYPE_CONNECTS_PATHS. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedMember);

  if (pNumClosing)
    *pNumClosing = 0;
  if (pNumShortcuts)
    *pNumShortcuts = 0;
  if (pNumExtending)
    *pNumExtending = 0;
  if (pNumConnecting)
    *pNumConnecting = 0;

  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    ReducedMember* child = reducedMember->children[c];
    assert(child);
    if (pNumClosing && child->type == TYPE_CLOSES_CYCLE)
      (*pNumClosing)++;
    else if (pNumShortcuts && child->type == TYPE_SHORTCUTS_PATH)
      (*pNumShortcuts)++;
    else if (pNumExtending && child->type == TYPE_EXTENDS_PATH)
      (*pNumExtending)++;
    else if (pNumConnecting && child->type == TYPE_CONNECTS_PATHS)
      (*pNumConnecting)++;
  }

  return TU_OKAY;
}

static
TU_ERROR determineBondType(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< New column of \p tdec. */
  TU_TDEC_MEMBER member         /**< Member of \p tdec. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("    Determining type of bond.\n");
#endif /* TU_DEBUG_TDEC */
  

  assert(false);
}

static
TU_ERROR determinePolygonType(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< New column of \p tdec. */
  TU_TDEC_MEMBER member         /**< Member of \p tdec. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("    Determining type of polygon.\n");
#endif /* TU_DEBUG_TDEC */

  assert(false);

  return TU_OKAY;
}

static
TU_ERROR determinePrimeType(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< New column of \p tdec. */
  TU_TDEC_MEMBER member         /**< Member of \p tdec. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("    Determining type of prime.\n");
#endif /* TU_DEBUG_TDEC */

  assert(false);

  return TU_OKAY;
}

static
TU_ERROR checkBondRoot(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn  /**< New column of \p tdec. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("    Checking bond root.\n");
#endif /* TU_DEBUG_TDEC */

  ReducedMember* root = &newcolumn->reducedMembers[0];
  if (root->numChildren == 0)
  {
    assert(root->firstReducedEdge != NULL); /* No children, no edge contradicts connected. */
    assert(root->firstReducedEdge->next == NULL); /* Two parallel tree edges. */

    TU_TDEC_EDGE edge = root->firstReducedEdge->edge;
    newcolumn->terminalNode1 = findEdgeHead(tdec, edge);
    newcolumn->terminalNode2 = findEdgeTail(tdec, edge);
    newcolumn->terminalMember1 = findEdgeMember(tdec, edge);
    newcolumn->terminalMember2 = newcolumn->terminalMember1;
    newcolumn->remainsGraphic = true;
#if defined(TU_DEBUG_TDEC)
    printf("        No children\n");
#endif /* TU_DEBUG_TDEC */
    return TU_OKAY;
  }

  assert(false);

  return TU_OKAY;
}

static
TU_ERROR checkPolygonRoot(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn  /**< new column. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("    Checking polygon root.\n");
#endif /* TU_DEBUG_TDEC */

  ReducedMember* root = &newcolumn->reducedMembers[0];
  int numClosing;
  int numShortcuts;
  int numExtending;
  TU_CALL( countChildrenTypes(tu, tdec, newcolumn, root, &numClosing, &numShortcuts, &numExtending,
    NULL) );

#if defined(TU_DEBUG_TDEC)
  printf("      It has %d closing children, %d shortcut children and %d extending children.\n",
    numClosing, numShortcuts, numExtending);
#endif /* TU_DEBUG_TDEC */

  newcolumn->remainsGraphic = (numClosing == 0) && (numExtending + numShortcuts <= 2);

  return TU_OKAY;
}

static
TU_ERROR checkPrimeRoot(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn  /**< New column of \p tdec. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("    Checking prime root.\n");
#endif /* TU_DEBUG_TDEC */

  assert(false);

  return TU_OKAY;
}

static
TU_ERROR addColumnBondSame(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new-column structure. */
  TU_TDEC_EDGE newEdge          /**< Edge. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("    Adding column for reduced t-decomposition with a bond root.\n");
#endif /* TU_DEBUG_TDEC */

  assert(newcolumn->terminalMember1 == newcolumn->terminalMember2);
  assert(findMember(tdec, newcolumn->terminalMember1) == newcolumn->terminalMember1);
  assert(findMember(tdec, newcolumn->terminalMember2) == newcolumn->terminalMember2);

  tdec->edges[newEdge].member = newcolumn->terminalMember1;
  tdec->edges[newEdge].head = newcolumn->terminalNode1;
  tdec->edges[newEdge].tail = newcolumn->terminalNode2;
  tdec->members[newcolumn->terminalMember1].numEdges++;

  return TU_OKAY;
}

static
TU_ERROR addColumnPolygonSame(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new-column structure. */
  TU_TDEC_EDGE newEdge          /**< Edge. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("    Adding column for reduced t-decomposition with a polygon root.\n");
#endif /* TU_DEBUG_TDEC */
  
  assert(false);

  return TU_OKAY;
}

static
TU_ERROR addColumnPrimeSame(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new-column structure. */
  TU_TDEC_EDGE newEdge          /**< Edge. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("    Adding column for reduced t-decomposition with a prime root.\n");
#endif /* TU_DEBUG_TDEC */
  
  assert(false);

  return TU_OKAY;
}

TU_ERROR TUtdecAddColumnCheck(TU* tu, TU_TDEC* tdec, TU_TDEC_NEWCOLUMN* newcolumn, int* entryRows,
  int numEntries)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(entryRows);
  assert(numEntries >= 1);
  
#if defined(TU_DEBUG_TDEC)
  printf("\n  Preparing to add a column with %d 1's.\n", numEntries);
#endif /* TU_DEBUG_TDEC */

  TU_CALL( initializeNewColumn(tu, tdec, newcolumn) );
  TU_CALL( findReducedDecomposition(tu, tdec, newcolumn, entryRows, numEntries) );
  TU_CALL( initializeReducedMemberEdgeLists(tu, tdec, newcolumn, entryRows, numEntries) );
  TU_CALL( computeReducedMemberChildren(tu, tdec, newcolumn) );

  TU_TDEC_MEMBER member;
  for (int i = newcolumn->numReducedMembers-1; i > 0; i--)
  {
    member = newcolumn->reducedMembers[i].member;
    /* Perform checks based on feedback from children (#TYPE1, etc.). */

    /* Determine the type. */
    if (tdec->members[member].type == TDEC_MEMBER_TYPE_BOND)
    {
      TU_CALL( determineBondType(tu, tdec, newcolumn, member) );
    }
    else if (tdec->members[member].type == TDEC_MEMBER_TYPE_POLYGON)
    {
      TU_CALL( determinePolygonType(tu, tdec, newcolumn, member) );
    }
    else
    {
      assert(tdec->members[member].type == TDEC_MEMBER_TYPE_PRIME);
      TU_CALL( determinePrimeType(tu, tdec, newcolumn, member) );
    }
  }

  /* Type check for root. */
  member = newcolumn->reducedMembers[0].member;
  if (tdec->members[member].type == TDEC_MEMBER_TYPE_BOND)
  {
    TU_CALL( checkBondRoot(tu, tdec, newcolumn) );
  }
  else if (tdec->members[member].type == TDEC_MEMBER_TYPE_POLYGON)
  {
    TU_CALL( checkPolygonRoot(tu, tdec, newcolumn) );
  }
  else
  {
    TU_CALL( checkPrimeRoot(tu, tdec, newcolumn) );
  }

  if (newcolumn->remainsGraphic)
  {
#if defined(TU_DEBUG_TDEC)
    printf("    Adding the column would maintain graphicness.\n");
    printf("    Co-tree edge from node %d {member %d} to node %d {member %d}.\n", newcolumn->terminalNode1,
      newcolumn->terminalMember1, newcolumn->terminalNode2, newcolumn->terminalMember2);
#endif /* TU_DEBUG_TDEC */
  }

  return TU_OKAY;
}

static
TU_ERROR createNewRowsPolygon(
  TU* tu,               /**< \ref TU environment. */
  TU_TDEC* tdec,        /**< t-decomposition. */
  TU_TDEC_EDGE* pedge,  /**< Pointer for storing the new edge. */
  TU_TDEC_NODE head,    /**< Head node. */
  TU_TDEC_NODE tail,    /**< Tail node. */
  int column,           /**< Index of new column to be added. */
  int* entryRows,       /**< Array of rows with 1-entry in this column. */
  int numEntries        /**< Number of 1-entries in this column. */
)
{
  assert(tu);
  assert(tdec);
  assert(column >= 0);
  assert(entryRows);
  assert(numEntries >= 0);

  /* Count new rows. */
  int countNewRows = 0;
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    if (row >= tdec->numRows || tdec->rowEdges[row].edge < 0)
      ++countNewRows;
  }
  
  if (countNewRows)
  {
#if defined(TU_DEBUG_TDEC)
    printf("    There are %d new rows.\n", countNewRows);
#endif /* TU_DEBUG_TDEC */

    /*
     * newEdge = parent marker edge.
     *
     * markerEdge ---------------> cotreeEdge
     *       ^                         |
     *       |                         |
     * first tree edge <- ... <- last tree edge
     * 
     * Arrow e --> f means that e->next = f, f->prev = e, e->head = f->tail.
     */

    TU_TDEC_MEMBER newMember;
    TU_CALL( createMember(tu, tdec, TDEC_MEMBER_TYPE_POLYGON, &newMember) );
    TU_TDEC_EDGE parentMarkerEdge;
    TU_CALL( createMarkerEdge(tu, tdec, &parentMarkerEdge, INT_MIN, head, tail, true) );
    tdec->edges[parentMarkerEdge].childMember = newMember;
    
    /* Add child marker edge and link it to marker edges. */
    TU_TDEC_NODE childMarkerHead = createNode(tu, tdec);
    TU_TDEC_NODE childMarkerTail = createNode(tu, tdec);
    TU_TDEC_EDGE childMarkerEdge;
    TU_CALL( createMarkerEdge(tu, tdec, &childMarkerEdge, newMember, childMarkerHead,
      childMarkerTail, false) );
    tdec->members[newMember].markerOfParent = parentMarkerEdge;
    tdec->members[newMember].markerToParent = childMarkerEdge;
    tdec->numMarkers++;
    
    /* Add new tree edges. */
    TU_TDEC_EDGE lastEdge = childMarkerEdge;
    for (int p = 0; p < numEntries; ++p)
    {
      int row = entryRows[p];
      if (row >= tdec->numRows || tdec->rowEdges[row].edge < 0)
      {
        TU_TDEC_NODE newTail = createNode(tu, tdec);
        TU_TDEC_EDGE treeEdge;
        TU_CALL( createRowEdge(tu, tdec, &treeEdge, newMember, tdec->edges[lastEdge].tail,
          newTail, row) );
        tdec->edges[lastEdge].prev = treeEdge;
        tdec->edges[treeEdge].next = lastEdge;
        lastEdge = treeEdge;
        tdec->members[newMember].numEdges++;
      }
    }

    /* Add cotree edge. */
    TU_TDEC_EDGE cotreeEdge;
    TU_CALL( createColumnEdge(tu, tdec, &cotreeEdge, newMember, tdec->edges[lastEdge].tail,
      childMarkerHead, column) );
    tdec->edges[childMarkerEdge].next = cotreeEdge;
    tdec->edges[cotreeEdge].prev = childMarkerEdge;
    tdec->edges[lastEdge].prev = cotreeEdge;
    tdec->edges[cotreeEdge].next = lastEdge;
    tdec->members[newMember].numEdges += 2;

    return parentMarkerEdge;
  }
  else
  {
    TU_CALL( createColumnEdge(tu, tdec, pedge, INT_MIN, head, tail, column) );
  }

  return TU_OKAY;
}

TU_ERROR TUtdecAddColumnApply(TU* tu, TU_TDEC* tdec, TU_TDEC_NEWCOLUMN* newcolumn, int column,
  int* entryRows, int numEntries)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(newcolumn->remainsGraphic);
  assert(newcolumn->numReducedMembers > 0);

#if defined(TU_DEBUG_TDEC)
  printf("  Adding a column with %d 1's.\n", numEntries);
#endif /* TU_DEBUG_TDEC */

  TU_TDEC_EDGE newEdge;
  TU_CALL( createNewRowsPolygon(tu, tdec, &newEdge, newcolumn->terminalNode1,
    newcolumn->terminalNode2, column, entryRows, numEntries) );
#if defined(TU_DEBUG_TDEC)
  printf("    New edge is %d.\n", newEdge);
#endif /* TU_DEBUG_TDEC */

  if (newcolumn->terminalMember1 == newcolumn->terminalMember2)
  {
    if (tdec->members[newcolumn->reducedMembers[0].member].type == TDEC_MEMBER_TYPE_BOND)
    {
      TU_CALL( addColumnBondSame(tu, tdec, newcolumn, newEdge) );
    }
    else if (tdec->members[newcolumn->reducedMembers[0].member].type == TDEC_MEMBER_TYPE_POLYGON)
    {
      TU_CALL( addColumnPolygonSame(tu, tdec, newcolumn, newEdge) );
    }
    else
    {
      assert(tdec->members[newcolumn->reducedMembers[0].member].type == TDEC_MEMBER_TYPE_PRIME);
      TU_CALL( addColumnPrimeSame(tu, tdec, newcolumn, newEdge) );
    }
    return TU_OKAY;
  }
  
  assert(false);

  return TU_OKAY;
}

TU_ERROR testGraphicnessTDecomposition(TU* tu, TU_CHRMAT* matrix, TU_CHRMAT* transpose,
  bool* isGraphic, TU_GRAPH* graph, TU_GRAPH_EDGE* basis, TU_GRAPH_EDGE* cobasis,
  TU_SUBMAT** psubmatrix)
{
  assert(tu);
  assert(matrix);
  assert(transpose);
  assert(!graph || (TUgraphNumNodes(graph) == 0 && TUgraphNumEdges(graph) == 0));
  assert(!psubmatrix || !*psubmatrix);
  assert(!basis || graph);
  assert(!cobasis || graph);

#if defined(TU_DEBUG_TDEC)
  printf("testGraphicnessTDecomposition called for a 1-connected %dx%d matrix.\n",
    matrix->numRows, matrix->numColumns);
  TUchrmatPrintDense(stdout, (TU_CHRMAT*) matrix, ' ', true);
#endif /* TU_DEBUG_TDEC */
  
  if (matrix->numNonzeros == 0)
  {
    *isGraphic = true;
    if (graph)
    {
      /* Construct a path with numRows edges and with numColumns loops at 0. */

      TU_GRAPH_NODE s;
      TU_CALL( TUgraphAddNode(tu, graph, &s) );
      for (int c = 0; c < matrix->numColumns; ++c)
      {
        TU_GRAPH_EDGE e;
        TU_CALL( TUgraphAddEdge(tu, graph, s, s, &e) );
        if (cobasis)
          *cobasis++ = e;
      }
      for (int r = 0; r < matrix->numRows; ++r)
      {
        TU_GRAPH_NODE t;
        TU_CALL( TUgraphAddNode(tu, graph, &t) );
        TU_GRAPH_EDGE e;
        TU_CALL( TUgraphAddEdge(tu, graph, s, t, &e) );
        if (basis)
          *basis++ = e;
        s = t;
      }

#if defined(TU_DEBUG_TDEC)
      printf("Constructed graph with %d nodes and %d edges.\n", TUgraphNumNodes(graph),
        TUgraphNumEdges(graph));
#endif /* TU_DEBUG_TDEC */
    }
    return TU_OKAY;
  }

  int rootRow = transpose->entryColumns[0];
#if defined(TU_DEBUG_TDEC)
  printf("  Root row is %d.\n", rootRow);
#endif /* TU_DEBUG_TDEC */
  TU_TDEC* tdec = NULL;
  TU_CALL( TUtdecCreate(tu, &tdec, rootRow, 0, 0, 0, 0, 0) ); /* TODO: avoid reallocations. */
  
  /* Process each column. */
  TU_TDEC_NEWCOLUMN* newcol = NULL;
  TUtdecnewcolumnCreate(tu, &newcol);
  for (int column = 0; column < matrix->numColumns; ++column)
  {
    TUtdecAddColumnCheck(tu, tdec, newcol,
      &transpose->entryColumns[transpose->rowStarts[column]],
      transpose->rowStarts[column+1] - transpose->rowStarts[column]);

    if (newcol->remainsGraphic)
    {
      TUtdecAddColumnApply(tu, tdec, newcol, column,
        &transpose->entryColumns[transpose->rowStarts[column]],
        transpose->rowStarts[column+1] - transpose->rowStarts[column]);
    }
    else
    {
      *isGraphic = false;
      assert(!"Not implemented");
    }
  }
  TUtdecnewcolumnFree(tu, &newcol);

  if (*isGraphic && graph)
  {
    TUtdecToGraph(tu, tdec, graph, true, basis, cobasis, NULL);
  }

  TUtdecFree(tu, &tdec);

  return TU_OKAY;
}
