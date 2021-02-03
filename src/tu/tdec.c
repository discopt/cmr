#define TU_DEBUG_TDEC /* Uncomment to enable debugging of t-decompositions. */

#include <tu/tdec.h>
#include "env_internal.h"

#include <assert.h>
#include <limits.h>

typedef struct
{
  TU_TDEC_NODE representativeNode;  /**< \brief Next representative of same node towards root, or -1 if root. */
} TU_TDEC_NODE_DATA;

typedef struct
{
  int name;                   /**< \brief Name of this edge.
                               *
                               * 0, 1, ..., m-1 indicate rows, -1,-2, ..., -n indicate columns,
                               * and for (small) k >= 0, MAX_INT-k and -MAX_INT+k indicate
                               * markers of the parent and to the parent, respectively. */
  TU_TDEC_MEMBER member;      /**< \brief Member this edge belongs to or -1 if in free list. */
  TU_TDEC_NODE head;          /**< \brief Head node of this edge for a prime member, -1 otherwise. */
  TU_TDEC_NODE tail;          /**< \brief Tail node of this edge for a prime member, -1 otherwise. */
  TU_TDEC_EDGE prev;          /**< \brief Next edge of this member. */
  TU_TDEC_EDGE next;          /**< \brief Previous edge of this member. */
  TU_TDEC_MEMBER childMember; /**< \brief Child member linked to this edge, or -1. */
} TU_TDEC_EDGE_DATA;

typedef struct
{
  TU_TDEC_MEMBER_TYPE type;             /**< \brief Type of member. Only valid for representative member. */
  TU_TDEC_MEMBER representativeMember;  /**< \brief Representative of member, or -1 if this is a representative member. */
  TU_TDEC_MEMBER parentMember;          /**< \brief Parent member of this member or -1 for a root. Only valid for representative member. */
  int numEdges;                         /**< \brief Number of edges. Only valid for representative member. */
  TU_TDEC_EDGE markerToParent;          /**< \brief Parent marker edge. Only valid for representative member. */
  TU_TDEC_EDGE markerOfParent;          /**< \brief Child marker of parent to which this member is linked. Only valid if root representative. */
  TU_TDEC_EDGE firstEdge;               /**< \brief First edge in doubly-linked edge list of this member. */
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

static inline
bool isRepresentativeMember(
  TU_TDEC* tdec,        /**< t-decomposition. */
  TU_TDEC_MEMBER member /**< Member of \p tdec. */
)
{
  assert(tdec);
  return tdec->members[member].representativeMember < 0;
}

/**
 * \brief Checks whether \p tdec has consistent edge data.
 * 
 * \returns Explanation of inconsistency, or \c NULL.
 */

static
char* consistencyEdges(
  TU* tu,       /**< \ref TU environment. */
  TU_TDEC* tdec /**< t-decomposition. */
)
{
  assert(tu);
  assert(tdec);

  for (TU_TDEC_MEMBER member = 0; member < tdec->numMembers; ++member)
  {
    if (!isRepresentativeMember(tdec, member))
      continue;

    TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
    int countEdges = 0;
    if (edge >= 0)
    {
      do
      {
        if (edge < 0 || edge >= tdec->memEdges)
          return TUconsistencyMessage("edge %d of member %d out of range.", member, edge);
        if (tdec->edges[edge].next < 0 || tdec->edges[edge].next > tdec->memEdges)
          return TUconsistencyMessage("edge %d of member %d has next out of range", member, edge);
        if (tdec->edges[tdec->edges[edge].next].prev != edge)
          return TUconsistencyMessage("member %d has inconsistent edge list", member);
        edge = tdec->edges[edge].next;
        countEdges++;
      }
      while (edge != tdec->members[member].firstEdge);
    }
    if (countEdges != tdec->members[member].numEdges)
    {
      return TUconsistencyMessage("member %d has %d edges, but numEdges %d", member, countEdges,
        tdec->members[member].numEdges);
    }
  }

  return NULL;
}

/**
 * \brief Checks whether \p tdec has consistent member data.
 * 
 * \returns Explanation of inconsistency, or \c NULL.
 */

static
char* consistencyMembers(
  TU* tu,       /**< \ref TU environment. */
  TU_TDEC* tdec /**< t-decomposition. */
)
{
  assert(tu);
  assert(tdec);

  for (TU_TDEC_MEMBER member = 0; member < tdec->numMembers; ++member)
  {
    if (tdec->members[member].type != TDEC_MEMBER_TYPE_BOND
      && tdec->members[member].type != TDEC_MEMBER_TYPE_PRIME
      && tdec->members[member].type != TDEC_MEMBER_TYPE_POLYGON)
      return TUconsistencyMessage("member %d has invalid type", member);
  }

  return NULL;
}

/**
 * \brief Checks whether \p tdec has consistent node data.
 * 
 * \returns Explanation of inconsistency, or \c NULL.
 */

static
char* consistencyNodes(
  TU* tu,       /**< \ref TU environment. */
  TU_TDEC* tdec /**< t-decomposition. */
)
{
  assert(tu);
  assert(tdec);

  for (TU_TDEC_MEMBER member = 0; member < tdec->numMembers; ++member)
  {
    bool isPrime = tdec->members[member].type == TDEC_MEMBER_TYPE_PRIME;
    TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
    if (edge < 0)
      continue;
    do
    {
      TU_TDEC_NODE head = tdec->edges[edge].head;
      TU_TDEC_NODE tail = tdec->edges[edge].tail;
      if (isPrime)
      {
        if (head < 0)
          return TUconsistencyMessage("edge %d of prime member %d has invalid head node", edge, member);
        if (tail < 0)
          return TUconsistencyMessage("edge %d of prime member %d has invalid tail node", edge, member);
        if (head >= tdec->memNodes)
          return TUconsistencyMessage("edge %d of prime member %d has head node out of range", edge, member);
        if (tail >= tdec->memNodes)
          return TUconsistencyMessage("edge %d of prime member %d has tail node out of range", edge, member);
      }
      else
      {
        if (head != -1)
          return TUconsistencyMessage("edge %d of non-prime member %d has a head node", edge, member);
        if (tail != -1)
          return TUconsistencyMessage("edge %d of non-prime member %d has a tail node", edge, member);
      }
      edge = tdec->edges[edge].next;
    }
    while (edge != tdec->members[member].firstEdge);
  }

  return NULL;
}

/**
 * \brief Checks whether the members of \p tdec form a forest.
 * 
 * \returns Explanation of inconsistency, or \c NULL.
 */

static
char* consistencyTree(
  TU* tu,       /**< \ref TU environment. */
  TU_TDEC* tdec /**< t-decomposition. */
)
{
  assert(tu);
  assert(tdec);

  for (TU_TDEC_MEMBER member = 0; member < tdec->numMembers; ++member)
  {
    if (member == 0)
      continue;

    int length = 0;
    TU_TDEC_MEMBER current;
    for (current = tdec->members[member].parentMember; current >= 0;
      current = tdec->members[current].parentMember)
    {
      ++length;
      if (length > tdec->numMembers)
        return "infinite member parent loop";
    }
  }

  return NULL;
}

typedef enum
{
  TYPE_UNKNOWN = 0,
  TYPE_1_HEAD_END_TAIL_END = 1, /**< Edge plus path is a cycle. */
  TYPE_2_HEAD_END_TAIL_IN = 2,  /**< Head is path end and tail is inner node. */
  TYPE_3_HEAD_END_TAIL_OUT = 3, /**< Head is path end and tail does not belong to path. */
  TYPE_4_HEAD_IN_TAIL_IN = 4,   /**< Head and tail are inner nodes such that adding the edge yields a path. */
  TYPE_ROOT = 5,                /**< Root member. */
  TYPE_5_OTHER = 6              /**< All other cases. */
} Type;

/**
 * \brief Additional edge information specific to a path.
 */

typedef struct _ReducedEdge
{
  TU_TDEC_EDGE edge;          /**< \brief The edge in the t-decomposition. */
  struct _ReducedEdge* next;  /**< \brief Next edge of this reduced member, or \c NULL. */
} ReducedEdge;

/**
 * \brief Additional member information specfic to a given path.
 */

typedef struct _ReducedMember
{
  TU_TDEC_MEMBER member;                  /**< \brief The member from the t-decomposition. */
  TU_TDEC_MEMBER rootMember;              /**< \brief The root member of this component of the t-decomposition. */
  int depth;                              /**< \brief Depth of this member in the reduced t-decomposition. */
  Type type;                              /**< \brief Type of this member. */
  int numChildren;                        /**< \brief Number of children in the reduced t-decomposition. */
  struct _ReducedMember** children;       /**< \brief Children in the reduced t-decomposition. */
  ReducedEdge* firstReducedEdge;          /**< \brief First edge in linked list of edges of this reduced member. */
  TU_TDEC_EDGE representativePathEdge;    /**< \brief Edge representing squeezed off path polygon. */
  TU_TDEC_EDGE representativeNonpathEdge; /**< \brief Edge representing squeezed off non-path polygon. */
} ReducedMember;

typedef struct _ReducedComponent
{
  int rootDepth;                  /**< \brief Depth of reduced root member. */
  ReducedMember* root;            /**< \brief Reduced root member. */
  TU_TDEC_NODE terminalNode1;     /**< \brief First terminal node of path. */
  TU_TDEC_NODE terminalNode2;     /**< \brief Second terminal node of path. */
  TU_TDEC_MEMBER terminalMember1; /**< \brief First terminal member of path. */
  TU_TDEC_MEMBER terminalMember2; /**< \brief Second terminal member of path. */
} ReducedComponent;

struct _TU_TDEC_NEWCOLUMN
{
  bool remainsGraphic;                      /**< \brief Indicator whether adding this column maintains graphicness. */
  int memReducedMembers;                    /**< \brief Allocated memory for \c reducedMembers. */
  int numReducedMembers;                    /**< \brief Number of members in \c reducedMembers. */
  ReducedMember* reducedMembers;            /**< \brief Array of reduced members, sorted by increasing depth. */
  ReducedMember** membersToReducedMembers;  /**< \brief Array mapping members to members of the reduced t-decomposition. */

  ReducedComponent* reducedComponents;      /**< \brief Array with reduced root members. */
  int memReducedComponents;                 /**< \brief Allocated memory for \c reducedComponents. */
  int numReducedComponents;                 /**< \brief Number of reduced root members. */

  ReducedEdge* reducedEdgeStorage;          /**< \brief Storage for edge lists of reduced members. */
  int memReducedEdgeStorage;                /**< \brief Allocated memory for \c reducedEdgeStorage. */
  int usedReducedEdgeStorage;               /**< \brief Number of stored edges in \c reducedEdgeStorage. */

  ReducedMember** childrenStorage;          /**< \brief Storage for members' arrays of children in reduced t-decomposition. */
  int usedChildrenStorage;                  /**< \brief Number of stored children in \c childrenStorage. */
  int memChildrenStorage;                   /**< \brief Allocated memory for \c childrenStorage. */

  int* nodesDegree;                         /**< \brief Map from nodes to degree w.r.t. path edges. */
  bool* edgesInPath;                        /**< \brief Map from edges to indicator for being in the path. */
};

int compareMemberDepths(const void* a, const void* b)
{
  const ReducedMember* first = a;
  const ReducedMember* second = b;
  /* Negative depths are moved to the end. */
  if (first->depth < 0)
    return +1;
  if (second->depth < 0)
    return -1;
  return first->depth - second->depth;
}

static TU_TDEC_MEMBER findMember(TU_TDEC* tdec, TU_TDEC_MEMBER start)
{
  TU_TDEC_MEMBER current = start;
  TU_TDEC_MEMBER next;
  while ((next = tdec->members[current].representativeMember) >= 0)
    current = next;
  TU_TDEC_MEMBER root = current;
  current = start;
  while ((next = tdec->members[current].representativeMember) >= 0)
  {
    if (next != root)
      tdec->members[current].representativeMember = root;
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

static // TODO: inline
TU_TDEC_MEMBER findEdgeMember(TU_TDEC* tdec, TU_TDEC_EDGE edge)
{
  return findMember(tdec, tdec->edges[edge].member);
}

/**
 * \brief Checks whether \p tdec has consistent parent/child structure of members.
 * 
 * \returns Explanation of inconsistency, or \c NULL.
 */

static
char* consistencyParentChild(
  TU* tu,       /**< \ref TU environment. */
  TU_TDEC* tdec /**< t-decomposition. */
)
{
  assert(tu);
  assert(tdec);

  if (tdec->memMembers < tdec->numMembers)
    return TUconsistencyMessage("member count and memory are inconsistent");
  if (tdec->numMembers < 0)
    return TUconsistencyMessage("negative member count");

  int* countChildren = NULL;
  if (TUallocStackArray(tu, &countChildren, tdec->memMembers) != TU_OKAY)
    return TUconsistencyMessage("stack allocation in consistencyParentChild() failed");
  for (int m = 0; m < tdec->memMembers; ++m)
    countChildren[m] = 0;

  for (TU_TDEC_MEMBER member = 0; member < tdec->numMembers; ++member)
  {
    if (tdec->members[member].representativeMember != -1)
      continue;

    if (tdec->members[member].parentMember >= tdec->memMembers)
    {
      TUfreeStackArray(tu, &countChildren);
      return TUconsistencyMessage("parent member of %d is out of range", member);
    }
    if (tdec->members[member].parentMember >= 0)
      countChildren[tdec->members[member].parentMember]++;
  }

  for (TU_TDEC_MEMBER member = 0; member < tdec->numMembers; ++member)
  {
    if (!isRepresentativeMember(tdec, member))
      continue;

    TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
    if (edge < 0)
      continue;
    do
    {
      if (tdec->edges[edge].childMember >= 0)
      {
        countChildren[member]--;
        
        if (tdec->members[tdec->edges[edge].childMember].parentMember != member)
        {
          TUfreeStackArray(tu, &countChildren);
          return TUconsistencyMessage("member %d has child edge %d for child %d whose parent member is %d",
            member, edge, tdec->edges[edge].childMember,
            tdec->members[tdec->edges[edge].childMember].parentMember);
        }
        if (tdec->members[tdec->edges[edge].childMember].markerOfParent != edge)
        {
          TUfreeStackArray(tu, &countChildren);
          return TUconsistencyMessage("member %d has child edge %d for child %d whose parent's markerOfParent is %d",
            member, edge, tdec->edges[edge].childMember,
            tdec->members[tdec->edges[edge].childMember].markerOfParent);
        }
        TU_TDEC_EDGE markerChild = tdec->members[tdec->edges[edge].childMember].markerToParent;
        if (tdec->edges[markerChild].name != -tdec->edges[edge].name)
        {
          TUfreeStackArray(tu, &countChildren);
          return TUconsistencyMessage("marker edges %d and %d of members %d (parent) and %d (child) have names %d and %d.",
            edge, markerChild, member, findEdgeMember(tdec, markerChild), tdec->edges[edge].name,
            tdec->edges[markerChild].name);
        }
      }
      edge = tdec->edges[edge].next;
    }
    while (edge != tdec->members[member].firstEdge);
  }

  if (TUfreeStackArray(tu, &countChildren) != TU_OKAY)
    return "stack deallocation in consistencyParentChild() failed";

  return NULL;
}

char* TUtdecConsistency(TU* tu, TU_TDEC* tdec)
{
  char* message = NULL;
  if ((message = consistencyMembers(tu, tdec)))
    return message;
  if ((message = consistencyEdges(tu, tdec)))
    return message;
  if ((message = consistencyNodes(tu, tdec)))
    return message;
  if ((message = consistencyParentChild(tu, tdec)))
    return message;
  if ((message = consistencyTree(tu, tdec)))
    return message;

  return NULL;
}

static
TU_TDEC_NODE findNode(TU_TDEC* tdec, TU_TDEC_NODE start)
{
  TU_TDEC_NODE current = start;
  TU_TDEC_NODE next;
  while ((next = tdec->nodes[current].representativeNode) >= 0)
    current = next;
  TU_TDEC_NODE root = current;
  current = start;
  while ((next = tdec->nodes[current].representativeNode) >= 0)
  {
    if (next != root)
      tdec->nodes[current].representativeNode = root;
    current = next;
  }
  return root;
}

static // TODO: inline
TU_TDEC_NODE findEdgeHead(TU_TDEC* tdec, TU_TDEC_EDGE edge)
{
  assert(edge >= 0);
  assert(edge < tdec->memEdges);
  assert(tdec->edges[edge].head >= 0);
  assert(tdec->edges[edge].head < tdec->memNodes);
  return findNode(tdec, tdec->edges[edge].head);
}

static // TODO: inline
TU_TDEC_NODE findEdgeTail(TU_TDEC* tdec, TU_TDEC_EDGE edge)
{
  return findNode(tdec, tdec->edges[edge].tail);
}

static
TU_ERROR createNode(
  TU* tu,             /**< \ref TU environment . */
  TU_TDEC* tdec,      /**< t-decomposition. */
  TU_TDEC_NODE* pnode /**< Pointer for storing new node. */
)
{
  assert(tu);
  assert(tdec);
  assert(pnode);

  TU_TDEC_NODE node = tdec->firstFreeNode;
  if (node >= 0)
  {
#if defined(TU_DEBUG_TDEC)
    printf("          createNode returns free node %d.\n", node);
#endif /* TU_DEBUG_TDEC */
    tdec->firstFreeNode = tdec->nodes[node].representativeNode;
  }
  else /* No member in free list, so we enlarge the array. */
  {
    int newSize = 2 * tdec->memNodes + 16;
    TU_CALL( TUreallocBlockArray(tu, &tdec->nodes, newSize) );
    for (int v = tdec->memNodes + 1; v < newSize; ++v)
      tdec->nodes[v].representativeNode = v+1;
    tdec->nodes[newSize-1].representativeNode = -1;
    tdec->firstFreeNode = tdec->memNodes + 1;
    node = tdec->memNodes;
    tdec->memNodes = newSize;
#if defined(TU_DEBUG_TDEC)
    printf("          createNode enlarges node array to %d and returns node %d.\n", newSize, node);
#endif /* TU_DEBUG_TDEC */
  }
  tdec->nodes[node].representativeNode = -1;
  tdec->numNodes++;

  *pnode = node;

  return TU_OKAY;
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

static
TU_ERROR setColumnEdge(
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
    TU_CALL( TUreallocBlockArray(tu, &tdec->columnEdges, newSize) );
    for (int c = tdec->memColumns; c < newSize; ++c)
      tdec->columnEdges[c].edge = -1;
    tdec->memColumns = newSize;
  }

  assert(tdec->columnEdges[column].edge == -1);
  tdec->columnEdges[column].edge = edge;
  if (column >= tdec->numColumns)
    tdec->numColumns = column + 1;

  return TU_OKAY;
}

/**
 * \brief Adds \p edge to the edge list of \p member.
 */

static
TU_ERROR addEdgeToMembersEdgeList(
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
  assert(isRepresentativeMember(tdec, member));

  TU_TDEC_EDGE first = tdec->members[member].firstEdge;
  if (first >= 0)
  {
    assert(tdec->members[member].numEdges > 0);
    TU_TDEC_EDGE last = tdec->edges[first].prev;
    tdec->edges[edge].next = first;
    tdec->edges[edge].prev = last;
    tdec->edges[first].prev = edge;
    tdec->edges[last].next =  edge;
  }
  else
  {
    assert(tdec->members[member].numEdges == 0);
    tdec->edges[edge].next = edge;
    tdec->edges[edge].prev = edge;
  }
  tdec->members[member].firstEdge = edge;
  tdec->members[member].numEdges++;

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
  assert(member < 0 || isRepresentativeMember(tdec, member));

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

#if defined(TU_DEBUG_TDEC)
  printf("        Created row edge %d <%d> of member %d.\n", edge, row, member);
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
  TU_CALL( setColumnEdge(tu, tdec, column, edge) );

#if defined(TU_DEBUG_TDEC)
  printf("        Created column edge {%d,%d} of member %d for column %d.\n", head, tail, member,
    column);
#endif /* TU_DEBUG_TDEC */

  return TU_OKAY;
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
    data->name = -INT_MAX + tdec->numMarkers;

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

  if (tdec->numMembers == tdec->memMembers)
  {
    tdec->memMembers = 16 + 2 * tdec->memMembers;
    TU_CALL( TUreallocBlockArray(tu, &tdec->members, tdec->memMembers) );
  }

  TU_TDEC_MEMBER_DATA* data = &tdec->members[tdec->numMembers];
  data->markerOfParent = -1;
  data->markerToParent = -1;
  data->firstEdge = -1;
  data->representativeMember = -1;
  data->numEdges = 0;
  data->parentMember = -1;
  data->type = type;
  *pmember = tdec->numMembers;
  tdec->numMembers++;

#if defined(TU_DEBUG_TDEC)
  printf("        Creating %s member %d\n", type == TDEC_MEMBER_TYPE_BOND ? "bond" :
    (type == TDEC_MEMBER_TYPE_PRIME ? "prime" : "polygon"), *pmember);
#endif /* TU_DEBUG_TDEC */

  return TU_OKAY;
}

TU_ERROR TUtdecCreate(TU* tu, TU_TDEC** ptdec, int memEdges, int memNodes, int memMembers,
  int numRows, int numColumns)
{
  assert(tu);
  assert(ptdec);
  assert(!*ptdec);

  TU_CALL( TUallocBlock(tu, ptdec) );
  TU_TDEC* tdec = *ptdec;
  if (memMembers < numRows)
    memMembers = numRows;
  tdec->memMembers = memMembers;
  tdec->numMembers = 0;
  tdec->members = NULL;
  TU_CALL( TUallocBlockArray(tu, &tdec->members, tdec->memMembers) );

  if (memNodes < 1)
    memNodes = 1;
  tdec->memNodes = memNodes;
  tdec->nodes = NULL;
  TU_CALL( TUallocBlockArray(tu, &tdec->nodes, memNodes) );
  tdec->numNodes = 0;
  for (int v = 0; v < memNodes; ++v)
    tdec->nodes[v].representativeNode = v+1;
  tdec->nodes[memNodes-1].representativeNode = -1;
  tdec->firstFreeNode = 0;

  if (memEdges < 1)
    memEdges = 1;
  tdec->memEdges = memEdges;
  tdec->edges = NULL;
  TU_CALL( TUallocBlockArray(tu, &tdec->edges, memEdges) );
  tdec->numEdges = 0;
  tdec->numMarkers = 0;

  /* Initialize free list with unused edges. */
  if (memEdges > tdec->numEdges)
  {
    for (int e = tdec->numEdges; e < memEdges; ++e)
    {
      tdec->edges[e].next = e+1;
      tdec->edges[e].member = -1;
    }
    tdec->edges[memEdges-1].next = -1;
    tdec->firstFreeEdge = tdec->numEdges;
  }
  else
    tdec->firstFreeEdge = -1;

  tdec->numRows = numRows;
  tdec->memRows = numRows;
  tdec->rowEdges = NULL;
  TU_CALL( TUallocBlockArray(tu, &tdec->rowEdges, tdec->numRows) );
  for (int r = 0; r < tdec->numRows; ++r)
    tdec->rowEdges[r].edge = -1;

  tdec->numColumns = numColumns;
  tdec->memColumns = tdec->numColumns;
  tdec->columnEdges = NULL;
  TU_CALL( TUallocBlockArray(tu, &tdec->columnEdges, tdec->numColumns) );
  for (int c = 0; c < tdec->numColumns; ++c)
    tdec->columnEdges[c].edge = -1;

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

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

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

#if defined(TU_DEBUG_TDEC)
  printf("TUtdecToGraph for t-decomposition.\n");
#endif /* TU_DEBUG_TDEC */

  TU_CALL( TUgraphClear(tu, graph) );

  TU_GRAPH_EDGE* localEdgeElements = NULL;
  if (edgeElements)
    localEdgeElements = edgeElements;
  else if (basis || cobasis)
    TU_CALL( TUallocStackArray(tu, &localEdgeElements, tdec->memEdges) );
  TU_GRAPH_NODE* tdecNodesToGraphNodes = NULL;
  TU_CALL( TUallocStackArray(tu, &tdecNodesToGraphNodes, tdec->numNodes) );
  TU_GRAPH_EDGE* tdecEdgesToGraphEdges = NULL;
  TU_CALL( TUallocStackArray(tu, &tdecEdgesToGraphEdges, tdec->memEdges) );

  for (int v = 0; v < tdec->memNodes; ++v)
  {
    if (tdec->nodes[v].representativeNode < 0)
    {
      TU_CALL( TUgraphAddNode(tu, graph, &tdecNodesToGraphNodes[v]) );
    }
    else
      tdecNodesToGraphNodes[v] = -1;
  }

  for (TU_TDEC_MEMBER member = 0; member < tdec->numMembers; ++member)
  {
    if (!isRepresentativeMember(tdec, member))
      continue;

    TU_TDEC_MEMBER_TYPE type = tdec->members[member].type;
#if defined(TU_DEBUG_TDEC)
    printf("  Member %d is %s with %d edges.\n", member, type == TDEC_MEMBER_TYPE_BOND ?
      "a bond" : (type == TDEC_MEMBER_TYPE_POLYGON ? "a polygon" : "prime"),
       tdec->members[member].numEdges);
#endif /* TU_DEBUG_TDEC */

    TU_GRAPH_EDGE graphEdge;
    TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
    if (type == TDEC_MEMBER_TYPE_PRIME)
    {
      do
      {
        TU_TDEC_NODE head = findEdgeHead(tdec, edge);
        TU_TDEC_NODE tail = findEdgeTail(tdec, edge);
        TU_CALL( TUgraphAddEdge(tu, graph, tdecNodesToGraphNodes[head], tdecNodesToGraphNodes[tail],
          &graphEdge) );
        tdecEdgesToGraphEdges[edge] = graphEdge;
        if (localEdgeElements)
          localEdgeElements[graphEdge] = tdec->edges[edge].name;
        edge = tdec->edges[edge].next;
      }
      while (edge != tdec->members[member].firstEdge);
    }
    else if (type == TDEC_MEMBER_TYPE_BOND)
    {
      TU_GRAPH_NODE graphHead, graphTail;
      TU_CALL( TUgraphAddNode(tu, graph, &graphHead) );
      TU_CALL( TUgraphAddNode(tu, graph, &graphTail) );
      do
      {
        TU_CALL( TUgraphAddEdge(tu, graph, graphHead, graphTail, &graphEdge) );
        tdecEdgesToGraphEdges[edge] = graphEdge;
        if (localEdgeElements)
          localEdgeElements[graphEdge] = tdec->edges[edge].name;
        edge = tdec->edges[edge].next;
      }
      while (edge != tdec->members[member].firstEdge);
    }
    else
    {
      assert(type == TDEC_MEMBER_TYPE_POLYGON);

      TU_GRAPH_NODE firstNode, v;
      TU_CALL( TUgraphAddNode(tu, graph, &firstNode) );
      v = firstNode;
      edge = tdec->edges[edge].next;
      while (edge != tdec->members[member].firstEdge)
      {
        
        TU_GRAPH_NODE w;
        TU_CALL( TUgraphAddNode(tu, graph, &w) );
        TU_CALL( TUgraphAddEdge(tu, graph, v, w, &graphEdge) );
        tdecEdgesToGraphEdges[edge] = graphEdge;
        if (localEdgeElements)
          localEdgeElements[graphEdge] = tdec->edges[edge].name;

        edge = tdec->edges[edge].next;
        v = w;
      }
      TU_CALL( TUgraphAddEdge(tu, graph, v, firstNode, &graphEdge) );
      tdecEdgesToGraphEdges[edge] = graphEdge;
      if (localEdgeElements)
        localEdgeElements[graphEdge] = tdec->edges[edge].name;
    }
  }

  /* Merge respective parent and child edges. */

  if (merge)
  {
#if defined(TU_DEBUG_TDEC)
    printf("  Before merging, the graph has %d nodes and %d edges.\n", TUgraphNumNodes(graph),
      TUgraphNumEdges(graph));
    fflush(stdout);
#endif /* TU_DEBUG_TDEC */

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

#if defined(TU_DEBUG_TDEC)
      printf("  Merging edges %d = {%d,%d} <%d> and %d = {%d,%d} <%d>.\n", parent, parentU, parentV,
        tdec->edges[tdec->members[m].markerOfParent].name, child, childU, childV,
        tdec->edges[tdec->members[m].markerToParent].name);
      fflush(stdout);
#endif /* TU_DEBUG_TDEC */

      TU_CALL( TUgraphMergeNodes(tu, graph, parentU, childU) );
      TU_CALL( TUgraphDeleteNode(tu, graph, childU) );
      TU_CALL( TUgraphMergeNodes(tu, graph, parentV, childV) );
      TU_CALL( TUgraphDeleteNode(tu, graph, childV) );

      TU_CALL( TUgraphDeleteEdge(tu, graph, parent) );
      TU_CALL( TUgraphDeleteEdge(tu, graph, child) );
    }
  }

  // TODO: Remove nodes with degree 0 or 1?!

  /* Construct (co)basis. */

  if (basis || cobasis)
  {
#if !defined(NDEBUG)
    /* This is only relevant if a 1-separation exists. */
    for (int r = 0; r < tdec->numRows; ++r)
      basis[r] = INT_MIN;
    for (int c = 0; c < tdec->numColumns; ++c)
      cobasis[c] = INT_MIN;
#endif /* !NDEBUG */

    for (TU_GRAPH_ITER i = TUgraphEdgesFirst(graph); TUgraphEdgesValid(graph, i);
      i = TUgraphEdgesNext(graph, i))
    {
      TU_GRAPH_EDGE e = TUgraphEdgesEdge(graph, i);

#if defined(TU_DEBUG_TDEC)
      printf("  Graph edge %d = {%d,%d} <%d>\n", e, TUgraphEdgeU(graph, e), TUgraphEdgeV(graph, e),
        localEdgeElements[e]);
      fflush(stdout);
#endif /* TU_DEBUG_TDEC */

      int element = localEdgeElements[e];
      if (element >= 0 && basis)
        basis[element] = e;
      else if (element < 0 && cobasis)
        cobasis[-1-element] = e;
    }

#if !defined(NDEBUG)
    /* These assertions indicate a 1-separable input matrix. */
    for (int r = 0; r < tdec->numRows; ++r)
      assert(basis[r] >= 0);
    for (int c = 0; c < tdec->numColumns; ++c)
      assert(cobasis[c] >= 0);
#endif /* !NDEBUG */
  }

  TU_CALL( TUfreeStackArray(tu, &tdecEdgesToGraphEdges) );
  TU_CALL( TUfreeStackArray(tu, &tdecNodesToGraphNodes) );
  if (localEdgeElements != edgeElements)
    TU_CALL( TUfreeStackArray(tu, &localEdgeElements) );

  return TU_OKAY;
}

static
void edgeToDot(
  FILE* stream,
  TU_TDEC* tdec,
  TU_TDEC_MEMBER member,
  TU_TDEC_EDGE edge,
  int u,
  int v,
  bool red
)
{
  assert(stream);
  assert(member >= 0);
  assert(edge >= 0);

  const char* redStyle = red ? ",color=red" : "";
  if (tdec->members[member].markerToParent == edge)
  {
    fprintf(stream, "    %d.%d -- p%d [style=dashed%s];\n", member, u, member, redStyle);
    fprintf(stream, "    p%d -- %d.%d [style=dashed%s];\n", member, member, v, redStyle);
    fprintf(stream, "    %d.%d [shape=box];\n", member, u);
    fprintf(stream, "    %d.%d [shape=box];\n", member, v);
    fprintf(stream, "    p%d [style=dashed];\n", member);
  }
  else if (tdec->edges[edge].childMember >= 0)
  {
    TU_TDEC_MEMBER child = tdec->edges[edge].childMember;
    fprintf(stream, "    %d.%d -- c%d [style=dotted%s];\n", member, u, child, redStyle);
    fprintf(stream, "    c%d -- %d.%d [style=dotted%s];\n", child, member, v, redStyle);
    fprintf(stream, "    %d.%d [shape=box];\n", member, u);
    fprintf(stream, "    %d.%d [shape=box];\n", member, v);
    fprintf(stream, "    c%d [style=dotted];\n", child);

    fprintf(stream, "    p%d -- c%d [style=dashed,dir=forward];\n", child, child);
  }
  else
  {
    fprintf(stream, "    %d.%d -- %d.%d [label=\"%d <%d>\",style=bold%s];\n", member, u, member, v,
      edge, tdec->edges[edge].name, redStyle);
    fprintf(stream, "    %d.%d [shape=box];\n", member, u);
    fprintf(stream, "    %d.%d [shape=box];\n", member, v);
  }
}

TU_ERROR TUtdecToDot(TU* tu, TU_TDEC* tdec, FILE* stream, bool* edgesHighlighted)
{
  assert(tu);
  assert(tdec);
  assert(stream);

  fprintf(stream, "// t-decomposition\n");
  fprintf(stream, "graph tdec {\n");
  fprintf(stream, "  compound = true;\n");
  for (TU_TDEC_MEMBER member = 0; member < tdec->numMembers; ++member)
  {
    fprintf(stream, "  subgraph member%d {\n", member);
    if (tdec->members[member].type == TDEC_MEMBER_TYPE_BOND)
    {
      TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
      do
      {
        edgeToDot(stream, tdec, member, edge, 0, 1, edgesHighlighted ? edgesHighlighted[edge] : false);
        edge = tdec->edges[edge].next;
      }
      while (edge != tdec->members[member].firstEdge);
    }
    else if (tdec->members[member].type == TDEC_MEMBER_TYPE_PRIME)
    {
      TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
      do
      {
        TU_TDEC_NODE u = findEdgeHead(tdec, edge);
        TU_TDEC_NODE v = findEdgeTail(tdec, edge);
        edgeToDot(stream, tdec, member, edge, u, v,
          edgesHighlighted ? edgesHighlighted[edge] : false);
        edge = tdec->edges[edge].next;
      }
      while (edge != tdec->members[member].firstEdge);
    }
    else
    {
      assert(tdec->members[member].type == TDEC_MEMBER_TYPE_POLYGON);
      TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
      int i = 0;
      do
      {
        edgeToDot(stream, tdec, member, edge, i, (i+1) % tdec->members[member].numEdges,
          edgesHighlighted ? edgesHighlighted[edge] : false);
        edge = tdec->edges[edge].next;
        i++;
      }
      while (edge != tdec->members[member].firstEdge);
    }
    fprintf(stream, "  }\n");
  }
  fprintf(stream, "}\n");

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

  newcolumn->numReducedComponents = 0;
  newcolumn->memReducedComponents = 0;
  newcolumn->reducedComponents = NULL;

  newcolumn->reducedEdgeStorage = NULL;
  newcolumn->memReducedEdgeStorage = 0;
  newcolumn->usedReducedEdgeStorage = 0;

  newcolumn->memChildrenStorage = 0;
  newcolumn->usedChildrenStorage = 0;
  newcolumn->childrenStorage = NULL;

  newcolumn->nodesDegree = NULL;
  newcolumn->edgesInPath = NULL;

  return TU_OKAY;
}

TU_ERROR TUtdecnewcolumnFree(TU* tu, TU_TDEC_NEWCOLUMN** pnewcolumn)
{
  assert(tu);
  assert(*pnewcolumn);

  TU_TDEC_NEWCOLUMN* newcolumn = *pnewcolumn;
  
  if (newcolumn->edgesInPath)
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->edgesInPath) );

  if (newcolumn->nodesDegree)
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->nodesDegree) );

  if (newcolumn->reducedComponents)
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->reducedComponents) );
  if (newcolumn->reducedMembers)
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->reducedMembers) );

  if (newcolumn->membersToReducedMembers)
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->membersToReducedMembers) );
  if (newcolumn->reducedEdgeStorage)
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->reducedEdgeStorage) );
  if (newcolumn->childrenStorage)
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->childrenStorage) );

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
  assert(newcolumn->numReducedMembers == 0);

  newcolumn->remainsGraphic = true;
  newcolumn->usedReducedEdgeStorage = 0;

  // TODO: Remember sizes of these arrays.

  /* memEdges does not suffice since new edges can be created by squeezing off.
   * Each squeezing off introduces 4 new edges, and we might apply this twice for each polygon member. */
  TU_CALL( TUreallocBlockArray(tu, &newcolumn->edgesInPath, tdec->memEdges + 8*tdec->numMembers) );
  TU_CALL( TUreallocBlockArray(tu, &newcolumn->nodesDegree, tdec->memNodes) );

  return TU_OKAY;
}

/**
 * \brief Creates, if necessary, the reduced member for \p member and calls itself for the parent.
 */

static
ReducedMember* createReducedMember(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new column. */
  TU_TDEC_MEMBER member,              /**< Member to create reduced member for. */
  ReducedMember** rootDepthMinimizer  /**< Array mapping root members to the depth minimizer. */
)
{
#if defined(TU_DEBUG_TDEC)
  printf("      Attempting to create reduced member %d.\n", member);
#endif /* TU_DEBUG_TDEC */

  ReducedMember* reducedMember = newcolumn->membersToReducedMembers[member];
  if (reducedMember)
  {
    /* This member is a known reduced member. If we meet an existing path of low depth, we remember
     * that. */

#if defined(TU_DEBUG_TDEC)
    printf("      Reduced member exists.\n");
#endif /* TU_DEBUG_TDEC */

    if (!rootDepthMinimizer[reducedMember->rootMember] ||
      reducedMember->depth < rootDepthMinimizer[reducedMember->rootMember]->depth)
    {
#if defined(TU_DEBUG_TDEC)
      printf("      Updating depth to %d.\n", reducedMember->depth);
#endif /* TU_DEBUG_TDEC */
      rootDepthMinimizer[reducedMember->rootMember] = reducedMember;
    }
  }
  else
  {
    reducedMember = &newcolumn->reducedMembers[newcolumn->numReducedMembers];
    newcolumn->numReducedMembers++;
    newcolumn->membersToReducedMembers[member] = reducedMember;
    reducedMember->member = member;
    reducedMember->numChildren = 0;
#if defined(TU_DEBUG_TDEC)
      printf("      Reduced member is new.\n");
#endif /* TU_DEBUG_TDEC */

    TU_TDEC_MEMBER parentMember = findMemberParent(tdec, member);
    if (parentMember >= 0)
    {
      ReducedMember* parentReducedMember = createReducedMember(tu, tdec, newcolumn, parentMember,
        rootDepthMinimizer);
      reducedMember->depth = parentReducedMember->depth + 1;
      reducedMember->rootMember = parentReducedMember->rootMember;
      parentReducedMember->numChildren++;
    }
    else
    {
      reducedMember->depth = 0;
      reducedMember->rootMember = member;

      /* We found a new component. We temporarily store the root member and later compute the actual
         reduced root. */
      if (newcolumn->memReducedComponents == newcolumn->numReducedComponents)
      {
        newcolumn->memReducedComponents = 2 * newcolumn->memReducedComponents + 16;
        TU_CALL( TUreallocBlockArray(tu, &newcolumn->reducedComponents,
          newcolumn->memReducedComponents) );
      }
#if defined(TU_DEBUG_TDEC)
      printf("      Initializing the new reduced component %d.\n", newcolumn->numReducedComponents);
#endif /* TU_DEBUG_TDEC */

      newcolumn->reducedComponents[newcolumn->numReducedComponents].root = reducedMember;
      newcolumn->numReducedComponents++;
    }

#if defined(TU_DEBUG_TDEC)
      printf("      The root member of %d is %d.\n", reducedMember->member, reducedMember->rootMember);
#endif /* TU_DEBUG_TDEC */
  }
  return reducedMember;
}

static
TU_ERROR computeReducedDecomposition(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new column. */
  int* entryRows,               /**< Array of rows of new column's enries. */
  int numEntries                /**< Length of \p entryRows. */
)
{
  /* Identify all members on the path. For the induced sub-arborescence we also compute the
   * depths. After the computation, its root has depth pathRootDepth. */
#if defined(TU_DEBUG_TDEC)
  printf("    Computing reduced t-decomposition.\n");
#endif /* TU_DEBUG_TDEC */

  /* Enlarge members array. */
  if (newcolumn->memReducedMembers < tdec->numMembers + numEntries)
  {
    newcolumn->memReducedMembers = tdec->memMembers + numEntries;
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->reducedMembers, newcolumn->memReducedMembers) );
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->membersToReducedMembers,
      newcolumn->memReducedMembers) );
  }

  /* Initialize the mapping from members to reduced members. */
  for (int m = 0; m < tdec->numMembers; ++m)
    newcolumn->membersToReducedMembers[m] = NULL;

  ReducedMember** rootDepthMinimizer = NULL;
  TU_CALL( TUallocStackArray(tu, &rootDepthMinimizer, tdec->numMembers) );
  for (int m = 0; m < tdec->numMembers; ++m)
    rootDepthMinimizer[m] = NULL;
  newcolumn->numReducedMembers = 0;
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    TU_TDEC_EDGE edge = (row < tdec->numRows) ? tdec->rowEdges[row].edge : -1;
    printf("      Entry %d is row %d of %d and corresponds to edge %d.\n", p, row, tdec->numRows, edge);
    if (edge >= 0)
    {
      TU_TDEC_MEMBER member = findEdgeMember(tdec, edge);
#if defined(TU_DEBUG_TDEC)
      printf("      Edge %d exists and belongs to member %d.\n", edge, member);
#endif /* TU_DEBUG_TDEC */
      ReducedMember* reducedMember = createReducedMember(tu, tdec, newcolumn, member,
        rootDepthMinimizer);

      /* For the first edge of this member, we set the depth minimizer to the new reduced member. */
      if (!rootDepthMinimizer[reducedMember->rootMember])
      {
        rootDepthMinimizer[reducedMember->rootMember] = reducedMember;
      }
    }
  }

  /* We set the reduced roots according to the minimizers. */
  for (int i = 0; i < newcolumn->numReducedComponents; ++i)
  {
#if defined(TU_DEBUG_TDEC)
    printf("      Considering reduced component %d with root member %d.\n", i,
      newcolumn->reducedComponents[i].root->member);
    printf("      The minimizer is %d.\n", rootDepthMinimizer[newcolumn->reducedComponents[i].root->member]->member);
#endif /* TU_DEBUG_TDEC */
    newcolumn->reducedComponents[i].root =
      rootDepthMinimizer[newcolumn->reducedComponents[i].root->member];
//     newcolumn->reducedComponents[i].rootDepth =
//       rootDepthMinimizer[newcolumn->reducedComponents[i].root->member]->depth;
#if defined(TU_DEBUG_TDEC)
    printf("      Member %d is a reduced root.\n", newcolumn->reducedComponents[i].root->member);
#endif /* TU_DEBUG_TDEC */
  }

  /* Allocate memory for children. */
  if (newcolumn->memChildrenStorage < newcolumn->numReducedMembers)
  {
    newcolumn->memChildrenStorage = 2*newcolumn->numReducedMembers;
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->childrenStorage, newcolumn->memChildrenStorage) );
  }

  /* Set memory pointer of each reduced member. */
  for (int m = 0; m < newcolumn->numReducedMembers; ++m)
  {
    ReducedMember* reducedMember = &newcolumn->reducedMembers[m];
    if (reducedMember->depth < rootDepthMinimizer[reducedMember->rootMember]->depth)
    {
#if defined(TU_DEBUG_TDEC)
      printf("      Member %d's depth is smaller than its reduced root.\n", reducedMember->member);
#endif /* TU_DEBUG_TDEC */
      continue;
    }

    reducedMember->children = &newcolumn->childrenStorage[newcolumn->usedChildrenStorage];
    newcolumn->usedChildrenStorage += reducedMember->numChildren;
    reducedMember->numChildren = 0;
  }

#if defined(TU_DEBUG_TDEC)
  printf("    Total number of children is %d / %d.\n", newcolumn->usedChildrenStorage,
    newcolumn->memChildrenStorage);
  fflush(stdout);
#endif /* TU_DEBUG_TDEC */

  /* Set children of each reduced member. */
  for (int m = 0; m < newcolumn->numReducedMembers; ++m)
  {
    ReducedMember* reducedMember = &newcolumn->reducedMembers[m];
    if (reducedMember->depth <= rootDepthMinimizer[reducedMember->rootMember]->depth)
    {
#if defined(TU_DEBUG_TDEC)
      printf("      Member %d's depth is smaller than or equal to that of its reduced root.\n",
        reducedMember->member);
      fflush(stdout);
#endif /* TU_DEBUG_TDEC */
      continue;
    }

    TU_TDEC_MEMBER parentMember = tdec->members[newcolumn->reducedMembers[m].member].parentMember;
    ReducedMember* parentReducedMember = parentMember >= 0 ?
      newcolumn->membersToReducedMembers[parentMember] : NULL;
    if (parentReducedMember)
    {

#if defined(TU_DEBUG_TDEC)
      printf("      Reduced member %ld (= member %d) has %d (= member %d) as child %d.\n",
        (parentReducedMember - newcolumn->reducedMembers),
        parentReducedMember->member, m, newcolumn->reducedMembers[m].member, parentReducedMember->numChildren);
#endif /* TU_DEBUG_TDEC */
      parentReducedMember->children[parentReducedMember->numChildren] = &newcolumn->reducedMembers[m];
      parentReducedMember->numChildren++;
    }
  }

  TU_CALL( TUfreeStackArray(tu, &rootDepthMinimizer) );

  return TU_OKAY;
}

/**
 * \brief Creates members and reduced members of new edges.
 */

static
TU_ERROR completeReducedDecomposition(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new column. */
  int* entryRows,               /**< Array of rows of new column's enries. */
  int numEntries                /**< Length of \p entryRows. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);

  /* Check if we need new rows. */

  int newNumRows = tdec->numRows-1;
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    TU_TDEC_EDGE edge = (row < tdec->numRows) ? tdec->rowEdges[row].edge : -1;
    if (edge < 0)
    {
      if (row > newNumRows)
        newNumRows = row;
    }
  }
  newNumRows++;

#if defined(TU_DEBUG_TDEC)
  printf("    Completing reduced decomposition: increasing #rows from %d to %d.\n", tdec->numRows,
    newNumRows);
#endif /* TU_DEBUG_TDEC */  

  /* Create single-edge bond members for all new rows. */

  if (newNumRows > tdec->numRows)
  {
    if (newNumRows > tdec->memRows)
    {
      tdec->memRows = 2*newNumRows;
      TU_CALL( TUreallocBlockArray(tu, &tdec->rowEdges, tdec->memRows) );
    }

    for (int r = tdec->numRows; r < newNumRows; ++r)
    {
      TU_TDEC_MEMBER member;
      TU_CALL( createMember(tu, tdec, TDEC_MEMBER_TYPE_BOND, &member) );
      
      TU_TDEC_EDGE edge;
      TU_CALL( createEdge(tu, tdec, member, &edge) );
      TU_CALL( addEdgeToMembersEdgeList(tu, tdec, edge, member) );
      tdec->edges[edge].name = r;
      tdec->edges[edge].head = -1;
      tdec->edges[edge].tail = -1;
      tdec->edges[edge].childMember = -1;

#if defined(TU_DEBUG_TDEC)
      printf("    New row %d is edge %d of member %d.\n", r, edge, member);
#endif /* TU_DEBUG_TDEC */

      tdec->rowEdges[r].edge = edge;
    }
  }

  /* Create reduced members. */

  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    if (row >= tdec->numRows)
    {
      /* Edge and member for this row were created. We create the reduced member. */
      TU_TDEC_EDGE edge = tdec->rowEdges[row].edge;
      TU_TDEC_MEMBER member = findEdgeMember(tdec, edge);

#if defined(TU_DEBUG_TDEC)
      printf("    Creating reduced member for edge %d of member %d.\n", edge, member);
#endif /* TU_DEBUG_TDEC */

      assert(newcolumn->numReducedMembers < newcolumn->memReducedMembers);
      ReducedMember* reducedMember = &newcolumn->reducedMembers[newcolumn->numReducedMembers];
      newcolumn->numReducedMembers++;
      reducedMember->numChildren = 0;
      reducedMember->member = member;
      reducedMember->depth = 0;
      reducedMember->rootMember = -1;
      reducedMember->type = TYPE_ROOT;

      assert(newcolumn->usedReducedEdgeStorage + 1 < newcolumn->memReducedEdgeStorage);
      ReducedEdge* reducedEdge = &newcolumn->reducedEdgeStorage[newcolumn->usedReducedEdgeStorage];
      newcolumn->usedReducedEdgeStorage++;
      reducedEdge->next = NULL;
      reducedEdge->edge = edge;
      reducedMember->firstReducedEdge = reducedEdge;

      if (newcolumn->numReducedComponents == newcolumn->memReducedComponents)
      {
        newcolumn->memReducedComponents = 2 * newcolumn->memReducedComponents + 16;
        TU_CALL( TUreallocBlockArray(tu, &newcolumn->reducedComponents,
          newcolumn->memReducedComponents) );
      }

      newcolumn->membersToReducedMembers[member] = reducedMember;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].root = reducedMember;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].rootDepth = 0;
      newcolumn->numReducedComponents++;
    }
  }

  tdec->numRows = newNumRows;

#if defined(TU_DEBUG_TDEC)
  printf("    Number of reduced components is %d.\n", newcolumn->numReducedComponents);
#endif /* TU_DEBUG_TDEC */

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
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  
#if defined(TU_DEBUG_TDEC)
  printf("    Initializing edge lists for members of reduced t-decomposition.\n");
#endif /* TU_DEBUG_TDEC */

  for (int v = 0; v < tdec->memNodes; ++v)
    newcolumn->nodesDegree[v] = 0;
  for (int e = 0; e < tdec->memEdges; ++e)
    newcolumn->edgesInPath[e] = false;

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

      newcolumn->edgesInPath[edge] = true;
      if (tdec->members[member].type == TDEC_MEMBER_TYPE_PRIME)
      {
        newcolumn->nodesDegree[findEdgeHead(tdec, edge)]++;
        newcolumn->nodesDegree[findEdgeTail(tdec, edge)]++;
      }

#if defined(TU_DEBUG_TDEC)
      printf("      Edge %d <%d> belongs to reduced member %ld which is member %d.\n", edge,
        tdec->edges[edge].name, (reducedMember - newcolumn->reducedMembers), reducedMember->member);
#endif /* TU_DEBUG_TDEC */
    }
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
  ReducedMember* reducedMember, /**< Reduced member. */
  int* pNumOneEnd,              /**< Number of children that (recursively) must contain one path end. */
  int* pNumTwoEnds              /**< Number of children that (recursively) must contain two path ends. */
)
{
  assert(tu);
  assert(tdec);
  assert(reducedMember);

  if (pNumOneEnd)
    *pNumOneEnd = 0;
  if (pNumTwoEnds)
    *pNumTwoEnds = 0;

  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    ReducedMember* child = reducedMember->children[c];
    assert(child);
    if (pNumOneEnd && (child->type == TYPE_2_HEAD_END_TAIL_IN
      || child->type == TYPE_3_HEAD_END_TAIL_OUT
      || 0)) // child->type == TYPE_3_HEAD_OUT_TAIL_END))
    {
      (*pNumOneEnd)++;
    }
    else if (pNumTwoEnds && child->type == TYPE_4_HEAD_IN_TAIL_IN)
      (*pNumTwoEnds)++;
  }

  return TU_OKAY;
}

static
TU_ERROR determineTypes(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember        /**< Reduced member. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);

#if defined(TU_DEBUG_TDEC)
  printf("  determineTypes(reduced member %ld = member %d)\n",
    (reducedMember - &newcolumn->reducedMembers[0]), reducedMember->member);
#endif /* TU_DEBUG_TDEC */

  /* First handle children recursively. */
  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    TU_CALL( determineTypes(tu, tdec, newcolumn, reducedComponent, reducedMember->children[c]) );

#if defined(TU_DEBUG_TDEC)
    if (newcolumn->remainsGraphic)
      printf("    Child has type %d\n", reducedMember->children[c]->type);
    else
      printf("    Child prohibits graphicness.\n");
#endif /* TU_DEBUG_TDEC */

    /* Abort if some part indicates non-graphicness. */
    if (!newcolumn->remainsGraphic)
      return TU_OKAY;
  }

  int numOneEnd;
  int numTwoEnds;
  TU_CALL( countChildrenTypes(tu, tdec, reducedMember, &numOneEnd, &numTwoEnds) );

#if defined(TU_DEBUG_TDEC)
  printf("    It has %d children with one end and %d with two ends.\n", numOneEnd,
    numTwoEnds);
#endif /* TU_DEBUG_TDEC */
  
  if (2*numTwoEnds + numOneEnd > 2)
  {
    newcolumn->remainsGraphic = false;
    return TU_OKAY;
  }

  bool isRoot = reducedMember == reducedComponent->root;

  /* Different behavior for bonds, polygons and prime components. */
  TU_TDEC_MEMBER member = reducedMember->member;
  if (tdec->members[member].type == TDEC_MEMBER_TYPE_BOND)
  {
    if (isRoot)
    {
      /* A bond root should always work. */
      reducedMember->type = TYPE_ROOT;
    }
    else
    {
      // No children, but an reduced edge.
      if (2*numTwoEnds + numOneEnd == 0 && reducedMember->firstReducedEdge)
        reducedMember->type = TYPE_1_HEAD_END_TAIL_END;
      else
      {
        assert(0 == "Typing of non-root bond not fully implemented.");
      }
    }
  }
  else if (tdec->members[member].type == TDEC_MEMBER_TYPE_POLYGON)
  {
    if (isRoot)
    {
      newcolumn->remainsGraphic = (numTwoEnds == 0) && (numOneEnd <= 2);
    }
    else
    {
      int countReducedEdges = 0;
      for (ReducedEdge* edge = reducedMember->firstReducedEdge; edge != NULL; edge = edge->next)
        ++countReducedEdges;
      int numEdges = tdec->members[member].numEdges;
      if (countReducedEdges == numEdges - 1)
      {
        reducedMember->type = TYPE_1_HEAD_END_TAIL_END;
        return TU_OKAY;
      }
      else if (countReducedEdges + numTwoEnds == numEdges - 1)
      {
        assert(numTwoEnds == 1);
        reducedMember->type = TYPE_4_HEAD_IN_TAIL_IN;
        return TU_OKAY;
      }
      else if (numTwoEnds == 1)
      {
        newcolumn->remainsGraphic = false;
        return TU_OKAY;
      }
      else if (numOneEnd == 1)
      {
        reducedMember->type = TYPE_3_HEAD_END_TAIL_OUT;
        return TU_OKAY;
      }
      else if (numOneEnd == 2)
      {
        reducedMember->type = TYPE_4_HEAD_IN_TAIL_IN;
        return TU_OKAY;
      }
      else
      {
        assert(numOneEnd == 0);
        assert(numTwoEnds == 0);
        reducedMember->type = TYPE_3_HEAD_END_TAIL_OUT;
        return TU_OKAY;
      }
    }
  }
  else
  {
    assert(tdec->members[member].type == TDEC_MEMBER_TYPE_PRIME);

    assert(0 == "Typing of prime not implemented.");
  }

  if (!isRoot && reducedMember->type == TYPE_1_HEAD_END_TAIL_END)
  {
    TU_TDEC_MEMBER parentMember = tdec->members[reducedMember->member].parentMember;
    ReducedMember* reducedParent = newcolumn->membersToReducedMembers[parentMember];
    TU_TDEC_EDGE markerOfParent = tdec->members[member].markerOfParent;

#if defined(TU_DEBUG_TDEC)
    printf("    Marker edge closes cycle.\n");
    printf("    Parent member %d is reduced member %ld.\n", parentMember,
      (reducedParent - newcolumn->reducedMembers));
#endif /* TU_DEBUG_TDEC */

    /* Add marker edge of parent to reduced parent's reduced edges. */

    assert(newcolumn->usedReducedEdgeStorage < newcolumn->memReducedEdgeStorage);
    ReducedEdge* reducedEdge = &newcolumn->reducedEdgeStorage[newcolumn->usedReducedEdgeStorage];
    ++newcolumn->usedReducedEdgeStorage;
    reducedEdge->edge = markerOfParent;
    reducedEdge->next = reducedParent->firstReducedEdge;
    reducedParent->firstReducedEdge = reducedEdge;

    /* Indicate that marker edge of parent belongs to path. */
    newcolumn->edgesInPath[markerOfParent] = true;

    /* Increase node degrees of nodes in parent. */
    if (tdec->members[reducedParent->member].type == TDEC_MEMBER_TYPE_PRIME)
    {
      newcolumn->nodesDegree[findEdgeHead(tdec, markerOfParent)]++;
      newcolumn->nodesDegree[findEdgeTail(tdec, markerOfParent)]++;
    }

#if defined(TU_DEBUG_TDEC)
    printf("    Added marker edge of parent to list of reduced edges.\n");
#endif /* TU_DEBUG_TDEC */
  }

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

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

  TU_CALL( initializeNewColumn(tu, tdec, newcolumn) );
  TU_CALL( computeReducedDecomposition(tu, tdec, newcolumn, entryRows, numEntries) );
  TU_CALL( initializeReducedMemberEdgeLists(tu, tdec, newcolumn, entryRows, numEntries) );

  for (int i = 0; i < newcolumn->numReducedComponents; ++i)
  {
    TU_CALL( determineTypes(tu, tdec, newcolumn, &newcolumn->reducedComponents[i],
      newcolumn->reducedComponents[i].root) );
  }

  if (newcolumn->remainsGraphic)
  {
#if defined(TU_DEBUG_TDEC)
    printf("    Adding the column would maintain graphicness.\n");
#endif /* TU_DEBUG_TDEC */
  }

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

  return TU_OKAY;
}

static
TU_ERROR addColumnPreprocessBond(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  bool isRoot,                        /**< Whether \p reducedMember is the reduced root. */
  int* pNumAssignedTerminals          /**< Pointer to number of assigned terminal nodes. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedMember);
  assert(pNumAssignedTerminals);

  int numOneEnd;
  int numTwoEnds;
  TU_CALL( countChildrenTypes(tu, tdec, reducedMember, &numOneEnd, &numTwoEnds) );
  
#if defined(TU_DEBUG_TDEC)
  printf("        addColumnPreprocessBond for reduced%s member %ld (member %d), #one-ends = %d, #two-ends = %d.\n",
    isRoot ? " root" : "", (reducedMember - newcolumn->reducedMembers),
    reducedMember->member, numOneEnd, numTwoEnds);
#endif /* TU_DEBUG_TDEC */

  if (isRoot && numOneEnd == 0 && numTwoEnds == 0)
  {
    assert(reducedMember->firstReducedEdge);
    assert(*pNumAssignedTerminals == 0);
    reducedComponent->terminalMember1 = reducedMember->member;
    reducedComponent->terminalMember2 = reducedMember->member;
    *pNumAssignedTerminals = 2;

    return TU_OKAY;
  }
  else if (isRoot && numOneEnd == 2)
  {
    assert(*pNumAssignedTerminals == 2);

    if (tdec->members[reducedMember->member].numEdges >= 4)
    {
      assert(0 == "addColumnPreprocessBond for root with two 1-ends and >= 4 edges is not implemented.");
      
      /* Do not forget to mark new bond as reduced root! */
    }

    return TU_OKAY;
  }

  assert(0 == "addColumnPreprocessBond is not implemented.");

  return TU_OKAY;
}

static
TU_ERROR addColumnPreprocessPrime(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  bool isRoot,                        /**< Whether \p reducedMember is the reduced root. */
  int* pNumAssignedTerminals          /**< Pointer to number of assigned terminal nodes. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedMember);
  assert(pNumAssignedTerminals);

  int numOneEnd;
  int numTwoEnds;
  TU_CALL( countChildrenTypes(tu, tdec, reducedMember, &numOneEnd, &numTwoEnds) );
  
#if defined(TU_DEBUG_TDEC)
  printf("        addColumnPreprocessPrime for reduced%s member %ld (member %d), #one-ends = %d, #two-ends = %d.\n",
    isRoot ? " root" : "", (reducedMember - newcolumn->reducedMembers),
    reducedMember->member, numOneEnd, numTwoEnds);
#endif /* TU_DEBUG_TDEC */

  assert(0 == "addColumnPreprocessPrime is not implemented.");

  return TU_OKAY;
}

/**
 * \brief Replaces an edge by a bond containing it.
 * 
 * The given member should have at least two edges.
 */

static
TU_ERROR createEdgeBond(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new column. */
  TU_TDEC_MEMBER member,        /**< Polygon member to be squeezed. */
  TU_TDEC_EDGE edge,            /**< Edge. */
  TU_TDEC_EDGE* pChildEdge      /**< Pointer for storing the child marker edge to the new bond. */
)
{
  assert(tu);
  assert(tdec);
  assert(member >= 0);
  assert(member < tdec->memMembers);
  assert(pChildEdge);

#if defined(TU_DEBUG_TDEC)
  printf("    Creating bond for edge %d in member %d.\n", edge, member);
#endif /* TU_DEBUG_TDEC */

  TU_TDEC_MEMBER bond = -1;
  TU_CALL( createMember(tu, tdec, TDEC_MEMBER_TYPE_BOND, &bond) );
  tdec->members[bond].parentMember = member;
  tdec->members[bond].markerOfParent = edge;

  TU_TDEC_EDGE markerOfParent;
  TU_CALL( createMarkerEdge(tu, tdec, &markerOfParent, member, tdec->edges[edge].head,
    tdec->edges[edge].tail, true) );
  tdec->edges[markerOfParent].childMember = bond;
  tdec->edges[markerOfParent].next = tdec->edges[edge].next;
  tdec->edges[markerOfParent].prev = tdec->edges[edge].prev;
  assert(tdec->edges[markerOfParent].next != markerOfParent);
  tdec->edges[tdec->edges[markerOfParent].next].prev = markerOfParent;
  tdec->edges[tdec->edges[markerOfParent].prev].next = markerOfParent;
  if (tdec->members[member].firstEdge == edge)
    tdec->members[member].firstEdge = markerOfParent;
  tdec->members[bond].markerOfParent = markerOfParent;

  TU_TDEC_EDGE markerToParent;
  TU_CALL( createMarkerEdge(tu, tdec, &markerToParent, bond, -1, -1, false) );
  TU_CALL( addEdgeToMembersEdgeList(tu, tdec, markerToParent, bond) );
  tdec->members[bond].markerToParent = markerToParent;
  tdec->numMarkers++;

  TU_CALL( addEdgeToMembersEdgeList(tu, tdec, edge, bond) );

  *pChildEdge = markerOfParent;

  return TU_OKAY;
}

/**
 * \brief Squeezes subset of polygon edges into a new polygon connected via a bond.
 * 
 * Takes all edges of the polygon \p member for which \p edgesPredicate is the same as
 * \p predicateValue.
 */

static
TU_ERROR squeezePolygonEdges(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new column. */
  TU_TDEC_MEMBER member,        /**< Polygon member to be squeezed. */
  bool* edgesPredicate,         /**< Map from edges to predicate. */
  bool predicateValue,          /**< Value of predicate. */
  TU_TDEC_EDGE* pChildEdge      /**< Pointer for storing the parent marker edge for the new bond. */
)
{
  assert(tu);
  assert(tdec);
  assert(member >= 0);
  assert(member < tdec->memMembers);
  assert(edgesPredicate);
  assert(tdec->members[member].type == TDEC_MEMBER_TYPE_POLYGON);
  assert(pChildEdge);

#if defined(TU_DEBUG_TDEC)
  printf("    Squeezing polygon %d.\n", member);
#endif /* TU_DEBUG_TDEC */

  /* Initialize new polygon. */
  TU_TDEC_MEMBER polygon;
  TU_CALL( createMember(tu, tdec, TDEC_MEMBER_TYPE_POLYGON, &polygon) );
  TU_TDEC_EDGE polygonParentMarker;
  TU_CALL( createMarkerEdge(tu, tdec, &polygonParentMarker, polygon, -1, -1, false) );
  TU_CALL( addEdgeToMembersEdgeList(tu, tdec, polygonParentMarker, polygon) );
  tdec->members[polygon].markerToParent = polygonParentMarker;

  /* Initialize new bond. */
  TU_TDEC_MEMBER bond;
  TU_CALL( createMember(tu, tdec, TDEC_MEMBER_TYPE_BOND, &bond) );
  TU_TDEC_EDGE bondChildMarker;
  TU_CALL( createMarkerEdge(tu, tdec, &bondChildMarker, bond, -1, -1, true) );
  TU_CALL( addEdgeToMembersEdgeList(tu, tdec, bondChildMarker, bond) );
  tdec->numMarkers++;
  TU_TDEC_EDGE bondParentMarker;
  TU_CALL( createMarkerEdge(tu, tdec, &bondParentMarker, bond, -1, -1, false) );
  TU_CALL( addEdgeToMembersEdgeList(tu, tdec, bondParentMarker, bond) );
  tdec->members[polygon].markerOfParent = bondChildMarker;
  tdec->members[bond].markerToParent = bondParentMarker;

  /* Go through old polygon. */

  TU_TDEC_EDGE firstEdge = tdec->members[member].firstEdge;
  TU_TDEC_EDGE edge = firstEdge;
  bool encounteredStayingEdge = false;
  do
  {
#if defined(TU_DEBUG_TDEC)
    printf("        Edge %d <%d>", edge, tdec->edges[edge].name);
    if (tdec->edges[edge].childMember >= 0)
      printf(" (with child %d)", tdec->edges[edge].childMember);
    if (edge == tdec->members[member].markerToParent)
      printf(" (with parent %d)", tdec->members[member].parentMember);
    printf(" (prev = %d, next = %d)", tdec->edges[edge].prev, tdec->edges[edge].next);
    fflush(stdout);
#endif /* TU_DEBUG_TDEC */
    /* Evaluate predicate. */
    bool value = edgesPredicate[edge];
    if ((value && !predicateValue) || (!value && predicateValue))
    {
#if defined(TU_DEBUG_TDEC)
      printf(" does not satisfy the predicate.\n");
      fflush(stdout);
#endif /* TU_DEBUG_TDEC */
      edge = tdec->edges[edge].next;
      encounteredStayingEdge = true;
      continue;
    }

#if defined(TU_DEBUG_TDEC)
    printf(" satisfies the predicate.\n");
    fflush(stdout);
#endif /* TU_DEBUG_TDEC */

    assert(edge != tdec->members[member].markerToParent);

    /* Remove edge from old edge list. */
    TU_TDEC_EDGE oldPrev = tdec->edges[edge].prev;
    TU_TDEC_EDGE oldNext = tdec->edges[edge].next;
    tdec->edges[oldPrev].next = oldNext;
    tdec->edges[oldNext].prev = oldPrev;
    tdec->members[member].numEdges--;

    /* Add edge to new edge list. */
    TU_TDEC_EDGE newPrev = tdec->edges[polygonParentMarker].prev;
    tdec->edges[newPrev].next = edge;
    tdec->edges[polygonParentMarker].prev = edge;
    tdec->edges[edge].prev = newPrev;
    tdec->edges[edge].next = polygonParentMarker;
    tdec->edges[edge].member = polygon;
    if (tdec->edges[edge].childMember >= 0)
    {
      assert( tdec->members[tdec->edges[edge].childMember].parentMember == member);
      tdec->members[tdec->edges[edge].childMember].parentMember = polygon;
    }
    tdec->members[polygon].numEdges++;

    /* Did we move the first edge of this member? */
    if (edge == firstEdge)
    {
      tdec->members[member].firstEdge = oldNext;
      firstEdge = oldNext;
      edge = oldNext;
      continue;
    }

    edge = oldNext;
  }
  while (edge != firstEdge || !encounteredStayingEdge);

  /* Add child marker edge from old polygon to bond and add it to edge list. */
  TU_TDEC_EDGE memberChildMarker;
  TU_CALL( createMarkerEdge(tu, tdec, &memberChildMarker, member, -1, -1, true) );
  tdec->numMarkers++;
  tdec->members[bond].markerOfParent = memberChildMarker;
  TU_TDEC_EDGE oldPrev = tdec->edges[firstEdge].prev;
  tdec->edges[memberChildMarker].next = firstEdge;
  tdec->edges[memberChildMarker].prev = oldPrev;
  tdec->edges[oldPrev].next = memberChildMarker;
  tdec->edges[firstEdge].prev = memberChildMarker;
  tdec->members[member].numEdges++;

  /* Link all. */
  tdec->members[polygon].parentMember = bond;
  tdec->edges[bondChildMarker].childMember = polygon;
  tdec->members[bond].parentMember = member;
  tdec->edges[memberChildMarker].childMember = bond;

#if defined(TU_DEBUG_TDEC)
  printf("        Updated old polygon:\n");
  edge = firstEdge;
  do
  {
    printf("          Edge %d <%d>.\n", edge, tdec->edges[edge].name);
    edge = tdec->edges[edge].next;
  }
  while (edge != firstEdge);

  printf("        New polygon:\n");
  edge = polygonParentMarker;
  do
  {
    printf("          Edge %d <%d>.\n", edge, tdec->edges[edge].name);
    edge = tdec->edges[edge].next;
  }
  while (edge != polygonParentMarker);
#endif /* TU_DEBUG_TDEC */

  *pChildEdge = memberChildMarker;

  return TU_OKAY;
}


static
TU_ERROR addColumnPreprocessPolygon(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  bool isRoot,                        /**< Whether \p reducedMember is the reduced root. */
  int* pNumAssignedTerminals          /**< Pointer to number of assigned terminal nodes. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedMember);
  assert(pNumAssignedTerminals);

  int numOneEnd;
  int numTwoEnds;
  TU_CALL( countChildrenTypes(tu, tdec, reducedMember, &numOneEnd, &numTwoEnds) );

#if defined(TU_DEBUG_TDEC)
  printf("        addColumnPreprocessPolygon for reduced%s member %ld (member %d), #one-ends = %d, #two-ends = %d.\n",
    isRoot ? " root" : "", (reducedMember - newcolumn->reducedMembers), reducedMember->member,
    numOneEnd, numTwoEnds);
#endif /* TU_DEBUG_TDEC */

  if (isRoot && numOneEnd == 0 && numTwoEnds == 0)
  {
    /* Root polygon containing both ends. */

    assert(reducedMember->firstReducedEdge);
    if (reducedMember->firstReducedEdge->next == NULL)
    {
      /* There is only one path edge, so we create a bond for that edge. */
      TU_TDEC_EDGE bondChildMarker;
      TU_CALL( createEdgeBond(tu, tdec, newcolumn, reducedMember->member,
        reducedMember->firstReducedEdge->edge, &bondChildMarker) );

      TU_TDEC_MEMBER bond = tdec->edges[bondChildMarker].childMember;
      reducedComponent->terminalMember1 = bond;
      reducedComponent->terminalMember2 = bond;
      *pNumAssignedTerminals = 2;
      return TU_OKAY;
    }
    else
    {
      /* Squeeze off all path edges by moving them to a new polygon and creating a bond to connect
       * it to the remaining polygon. */

      TU_TDEC_EDGE bondChildMarker;
      TU_CALL( squeezePolygonEdges(tu, tdec, newcolumn, reducedMember->member,
        newcolumn->edgesInPath, true, &bondChildMarker) );

      TU_TDEC_MEMBER bond = tdec->edges[bondChildMarker].childMember;
      reducedComponent->terminalMember1 = bond;
      reducedComponent->terminalMember2 = bond;
      *pNumAssignedTerminals = 2;
      return TU_OKAY;
    }
  }
  else if (!isRoot && reducedMember->type == TYPE_3_HEAD_END_TAIL_OUT
    && numOneEnd + numTwoEnds == 0)
  {
    /* Squeeze off all path edges by moving them to a new polygon and creating a bond to connect
     * it to the remaining polygon. */

    assert(*pNumAssignedTerminals < 2);
    assert(reducedMember->firstReducedEdge);

    /* Squeeze off path edges (if more than one). */
    if (reducedMember->firstReducedEdge->next)
    {
      TU_CALL( squeezePolygonEdges(tu, tdec, newcolumn, reducedMember->member,
        newcolumn->edgesInPath, true, &reducedMember->representativePathEdge) );
      newcolumn->edgesInPath[reducedMember->representativePathEdge] = true;
    }
    else
    {
      reducedMember->representativePathEdge = reducedMember->firstReducedEdge->edge;
    }
    printf("squeezedPathEdge = %d. Old polygon has length %d\n", reducedMember->representativePathEdge, tdec->members[reducedMember->member].numEdges);

    assert(tdec->members[reducedMember->member].numEdges >= 3);
    if (tdec->members[reducedMember->member].numEdges == 3)
    {
      reducedMember->representativeNonpathEdge
        = tdec->edges[tdec->members[reducedMember->member].markerToParent].next;
      if (reducedMember->representativeNonpathEdge == reducedMember->representativePathEdge)
      {
        reducedMember->representativeNonpathEdge
          = tdec->edges[reducedMember->representativeNonpathEdge].next;
      }
    }
    else
    {
      /* We temporarily mark the parent edge to belong to the path. */
      TU_TDEC_EDGE markerToParent = tdec->members[reducedMember->member].markerToParent;
      newcolumn->edgesInPath[markerToParent] = true;
      TU_CALL( squeezePolygonEdges(tu, tdec, newcolumn, reducedMember->member,
        newcolumn->edgesInPath, false, &reducedMember->representativeNonpathEdge) );
      newcolumn->edgesInPath[markerToParent] = false;
    }
    assert(tdec->members[reducedMember->member].numEdges == 3);

    printf("squeezedNonpathEdge = %d. Old polygon has length %d\n", reducedMember->representativeNonpathEdge, tdec->members[reducedMember->member].numEdges);

    if (*pNumAssignedTerminals == 0)
    {
      *pNumAssignedTerminals = 1;
      reducedComponent->terminalMember1 = reducedMember->member;
      return TU_OKAY;
    }
    else
    {
      *pNumAssignedTerminals = 2;
      reducedComponent->terminalMember2 = reducedMember->member;
      return TU_OKAY;
    }
  }

  assert(0 == "addColumnPreprocessPolygon is not implemented.");

  return TU_OKAY;
}

/**
 * \brief Preprocessing of reduced t-decomposition before the actual modification.
 * 
 * Processes the reduced members in depth-first search manner and does the following:
 * - Polygons are squeezed.
 * - Terminal nodes and (reduced) members are detected.
 */

static
TU_ERROR addColumnPreprocess(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int* pNumAssignedTerminals          /**< Pointer to number of assigned terminal nodes. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedComponent);

#if defined(TU_DEBUG_TDEC)
  printf("      addColumnPreprocess(reduced member %ld = member %d)\n",
    (reducedMember - &newcolumn->reducedMembers[0]), reducedMember->member);
#endif /* TU_DEBUG_TDEC */

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );
  
  /* If we are type 1, then we don't need to do anything. */
  if (reducedMember->type == TYPE_1_HEAD_END_TAIL_END)
  {
    return TU_OKAY;
  }

  /* Handle children recursively. */
  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    TU_CALL( addColumnPreprocess(tu, tdec, newcolumn, reducedComponent, reducedMember->children[c],
      pNumAssignedTerminals) );
  }

  bool isRoot = reducedMember == reducedComponent->root;

  /* Different behavior for bonds, polygons and prime components. */

  if (tdec->members[reducedMember->member].type == TDEC_MEMBER_TYPE_BOND)
  {
    TU_CALL( addColumnPreprocessBond(tu, tdec, newcolumn, reducedComponent, reducedMember, isRoot,
      pNumAssignedTerminals) );
  }
  else if (tdec->members[reducedMember->member].type == TDEC_MEMBER_TYPE_PRIME)
  {
    TU_CALL( addColumnPreprocessPrime(tu, tdec, newcolumn, reducedComponent, reducedMember, isRoot,
      pNumAssignedTerminals) );
  }
  else
  {
    assert(tdec->members[reducedMember->member].type == TDEC_MEMBER_TYPE_POLYGON);

    TU_CALL( addColumnPreprocessPolygon(tu, tdec, newcolumn, reducedComponent, reducedMember,
      isRoot, pNumAssignedTerminals) );
  }

  if (reducedMember->type == TYPE_1_HEAD_END_TAIL_END)
  {
    ReducedMember* parent = newcolumn->membersToReducedMembers[findMember(tdec,
      tdec->members[reducedMember->member].markerOfParent)];
    assert(parent);

    if (newcolumn->usedReducedEdgeStorage == newcolumn->memReducedEdgeStorage)
    {
      newcolumn->memReducedEdgeStorage *= 2;
      TU_CALL( TUreallocBlockArray(tu, &newcolumn->reducedEdgeStorage,
        newcolumn->memReducedEdgeStorage) );
      for (int e = newcolumn->usedReducedEdgeStorage; e < newcolumn->memReducedEdgeStorage; ++e)
        newcolumn->reducedEdgeStorage[e].next = &newcolumn->reducedEdgeStorage[e+1];
      
    }

    assert(0 == "Not tested.");
  }

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

#if 0
  
    if (isRoot)
    {
      if (*pNumAssignedTerminals == 0)
      {
        newcolumn->terminalMember1 = reducedMember;
        newcolumn->terminalNode1 = findEdgeHead(tdec, reducedMember->firstReducedEdge->edge);
        (*pNumAssignedTerminals)++;
      }
      if (*pNumAssignedTerminals == 1)
      {
        newcolumn->terminalMember2 = reducedMember;
        newcolumn->terminalNode2 = findEdgeTail(tdec, reducedMember->firstReducedEdge->edge);
        (*pNumAssignedTerminals)++;
      }
    }
    else
    {
      assert(0 == "addColumnPreprocess for non-root bond not yet implemented.");
    }
  }
  else if (tdec->members[reducedMember->member].type == TDEC_MEMBER_TYPE_POLYGON)
  {
#if defined(TU_DEBUG_TDEC)
  printf("    Member is a polygon.\n");
#endif /* TU_DEBUG_TDEC */

    /* The parent marker edge p (unless root). */
    TU_TDEC_EDGE markerToParent = tdec->members[reducedMember->member].markerToParent;
    TU_TDEC_EDGE start, edge;

    /* The child marker edges c1, c2 for the ends of the path. */
//     ReducedMember* endChild[2] = {NULL, NULL};
    TU_TDEC_EDGE endChildEdge[2] = {-1, -1};
//     bool endChildHead[2] = {false, false};
    for (int c = 0; c < reducedMember->numChildren; ++c)
    {
      Type type = reducedMember->children[c]->type;
      if (type == TYPE_2_HEAD_END_TAIL_IN || type == TYPE_3_HEAD_END_TAIL_OUT)
      {
        int childNum = endChildEdge[0] >= 0 ? 1 : 0;
        assert(endChildEdge[childNum] < 0);
//         endChild[childNum] = reducedMember->children[c];
        endChildEdge[childNum] = tdec->members[reducedMember->children[c]->member].markerOfParent;
//         endChildHead[childNum] = (type == TYPE_2_HEAD_END_TAIL_IN || type == TYPE_3_HEAD_END_TAIL_OUT );
      }
    }

    assert(endChildEdge[0] == -1 || 0 == "Not implemented: addColumnPreprocess with path-end children.");

    /* Reorder the edge list of the member as follows:
     *
     * [child-marker 0]
     * regular edge
     * regular edge
     * ...
     * regular edge
     * [parent marker]
     * [child-marker 1]
     * 
     * the last two are interchanged in case of a root.
     */

    /* Swap child-marker 1 with edge after parent marker. */
    if (endChildEdge[1] >= 0)
    {
      TU_CALL( swapEdges(tu, tdec, reducedMember->member, tdec->edges[markerToParent].next,
        endChildEdge[1]) );
    }

    /* Swap child marker 1 with parent marker if root. */
    if (endChildEdge[1] >= 0 && isRoot)
    {
      TU_CALL( swapEdges(tu, tdec, reducedMember->member, markerToParent, endChildEdge[1]) );
      edge = endChildEdge[1];
    }
    else
    {
      edge = markerToParent;
    }

    /* Swap regular edges. */
    for (ReducedEdge* reducedEdge = reducedMember->firstReducedEdge; reducedEdge;
      reducedEdge = reducedEdge->next)
    {
      TU_CALL( swapEdges(tu, tdec, reducedMember->member, tdec->edges[edge].prev,
        reducedEdge->edge) );
      edge = reducedEdge->edge;
    }

    /* Swap child-marker 0 with edge before last reduced edge. */
    if (endChildEdge[0] >= 0)
    {
      edge = tdec->edges[edge].prev;
      TU_CALL( swapEdges(tu, tdec, reducedMember->member, edge,
        endChildEdge[0]) );
    }

    /* Make this edge the first in the member's list. */
    tdec->members[reducedMember->member].firstEdge = edge;

#if defined(TU_DEBUG_TDEC)
    start = (endChildEdge[1] >= 0) ? tdec->edges[markerToParent].next : markerToParent;
    edge = start;
    do
    {
      if (edge == endChildEdge[0])
        printf("      1st end ");
      else if (edge == endChildEdge[1])
        printf("      2nd end ");
      else if (edge == markerToParent)
        printf("      Parent marker ");
      else
        printf("      Regular ");
      if (edge == tdec->members[reducedMember->member].firstEdge)
        printf("(first of this member) ");
      printf("edge %d <%d> {%d,%d}\n", edge, tdec->edges[edge].name, findEdgeHead(tdec, edge),
        findEdgeTail(tdec, edge));
      edge = tdec->edges[edge].prev;
    }
    while  (edge != start);
#endif /* TU_DEBUG_TDEC */

    if (isRoot)
    {
      if (*pNumAssignedTerminals == 0)
      {
        /* Ends of path in this polygon. */

        TU_TDEC_EDGE firstEdge = tdec->members[reducedMember->member].firstEdge;
        TU_TDEC_NODE head1 = findEdgeHead(tdec, firstEdge);
        TU_TDEC_NODE tail1 = findEdgeTail(tdec, firstEdge);
        TU_TDEC_NODE head2 = findEdgeHead(tdec, tdec->edges[firstEdge].next);
        TU_TDEC_NODE tail2 = findEdgeTail(tdec, tdec->edges[firstEdge].next);
        newcolumn->terminalNode1 = (head1 == head2 || head1 == tail2) ? tail1 : head1;

#if defined(TU_DEBUG_TDEC)
        printf("      First node of path is %d.\n", newcolumn->terminalNode1);
#endif /* TU_DEBUG_TDEC */

        TU_TDEC_EDGE lastEdge = tdec->edges[tdec->members[reducedMember->member].markerToParent].prev;
        head1 = findEdgeHead(tdec, lastEdge);
        tail1 = findEdgeTail(tdec, lastEdge);
        head2 = findEdgeHead(tdec, tdec->edges[lastEdge].prev);
        tail2 = findEdgeTail(tdec, tdec->edges[lastEdge].prev);
        newcolumn->terminalNode2 = (head1 == head2 || head1 == tail2) ? tail1 : head1;

#if defined(TU_DEBUG_TDEC)
        printf("      Last node of path is %d.\n", newcolumn->terminalNode2);
#endif /* TU_DEBUG_TDEC */

        newcolumn->terminalMember1 = reducedMember;
        newcolumn->terminalMember2 = reducedMember;
        (*pNumAssignedTerminals) = 2;
      }
      else
      {
        assert(0 == "addColumnPreprocess for root polygon with end-children not yet implemented.");
      }
    }
    else
    {
      assert(0 == "addColumnPreprocess for non-root polygon not yet implemented.");
    }
  }
  else
  {
    assert(tdec->members[reducedMember->member].type == TDEC_MEMBER_TYPE_PRIME);

#if defined(TU_DEBUG_TDEC)
    printf("    Member is prime.\n");
#endif /* TU_DEBUG_TDEC */

    if (isRoot)
    {
      assert(0 == "addColumnPreprocess for root prime not yet implemented.");
    }
    else
    {
      assert(0 == "addColumnPreprocess for non-root prime not yet implemented.");
    }
  }

#endif

  return TU_OKAY;
}

static
TU_ERROR createNewRowsPolygon(
  TU* tu,                   /**< \ref TU environment. */
  TU_TDEC* tdec,            /**< t-decomposition. */
  TU_TDEC_MEMBER* pmember,  /**< Pointer for storing the new member or -1 if there is none. */
  TU_TDEC_EDGE* pedge,      /**< Pointer for storing the new edge. */
  TU_TDEC_NODE head,        /**< Head node. */
  TU_TDEC_NODE tail,        /**< Tail node. */
  int column,               /**< Index of new column to be added. */
  int* entryRows,           /**< Array of rows with 1-entry in this column. */
  int numEntries            /**< Number of 1-entries in this column. */
)
{
  assert(tu);
  assert(tdec);
  assert(column >= 0);
  assert(entryRows);
  assert(numEntries >= 0);
  assert(pmember);
  assert(pedge);

#if defined(TU_DEBUG_TDEC)
  printf("      Creating polygon for new rows.\n");
#endif /* TU_DEBUG_TDEC */
  
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
    printf("      There are %d new rows.\n", countNewRows);
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

    TU_CALL( createMember(tu, tdec, TDEC_MEMBER_TYPE_POLYGON, pmember) );
    TU_TDEC_MEMBER newMember = *pmember;
    TU_TDEC_EDGE parentMarkerEdge;
    TU_CALL( createMarkerEdge(tu, tdec, &parentMarkerEdge, INT_MIN, head, tail, true) );
    tdec->edges[parentMarkerEdge].childMember = newMember;

    /* Add child marker edge and link it to marker edges. */
    TU_TDEC_EDGE childMarkerEdge;
    TU_CALL( createMarkerEdge(tu, tdec, &childMarkerEdge, newMember, -1, -1, false) );
    TU_CALL( addEdgeToMembersEdgeList(tu, tdec, childMarkerEdge, newMember) );
    tdec->members[newMember].markerOfParent = parentMarkerEdge;
    tdec->members[newMember].markerToParent = childMarkerEdge;
    tdec->numMarkers++;

    /* Add new tree edges. */
    for (int p = 0; p < numEntries; ++p)
    {
      int row = entryRows[p];
      if (row >= tdec->numRows || tdec->rowEdges[row].edge < 0)
      {
        TU_TDEC_EDGE treeEdge;
        TU_CALL( createRowEdge(tu, tdec, &treeEdge, newMember, -1, -1, row) );
        TU_CALL( addEdgeToMembersEdgeList(tu, tdec, treeEdge, newMember) );
      }
    }

    /* Add cotree edge. */
    TU_TDEC_EDGE cotreeEdge;
    TU_CALL( createColumnEdge(tu, tdec, &cotreeEdge, newMember, -1, -1, column) );
    TU_CALL( addEdgeToMembersEdgeList(tu, tdec, cotreeEdge, newMember) );

    *pedge = parentMarkerEdge;
  }
  else
  {
    *pmember = -1;
    TU_CALL( createColumnEdge(tu, tdec, pedge, INT_MIN, head, tail, column) );
  }

  return TU_OKAY;
}

/**
 * \brief Creates nodes for \p member if necessary.
 */

static
TU_ERROR createMemberNodes(
  TU* tu,               /**< \ref TU environment. */
  TU_TDEC* tdec,        /**< t-decomposition. */
  TU_TDEC_MEMBER member /**< member. */
)
{
  assert(tu);
  assert(tdec);
  assert(member >= 0 && member < tdec->memMembers);

  if (tdec->members[member].type == TDEC_MEMBER_TYPE_PRIME)
    return TU_OKAY;

#if defined(TU_DEBUG_TDEC)
  printf("        Creating nodes for member %d.\n", member);
#endif /* TU_DEBUG_TDEC */

  if (tdec->members[member].type == TDEC_MEMBER_TYPE_BOND)
  {
    assert(tdec->members[member].firstEdge >= 0);

    TU_TDEC_NODE head, tail;
    TU_CALL( createNode(tu, tdec, &head) );
    TU_CALL( createNode(tu, tdec, &tail) );

    TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
    do
    {
      assert(tdec->edges[edge].head < 0);
      tdec->edges[edge].head = head;
      assert(tdec->edges[edge].tail < 0);
      tdec->edges[edge].tail = tail;
      edge = tdec->edges[edge].next;
    }
    while (edge != tdec->members[member].firstEdge);
  }
  else if (tdec->members[member].type == TDEC_MEMBER_TYPE_POLYGON)
  {
    assert(tdec->members[member].firstEdge >= 0);

    TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
    do
    {
      TU_TDEC_NODE v;
      TU_CALL( createNode(tu, tdec, &v) );

      assert(tdec->edges[edge].head < 0);
      tdec->edges[edge].head = v;

      edge = tdec->edges[edge].next;

      assert(tdec->edges[edge].tail < 0);
      tdec->edges[edge].tail = v;
    }
    while (edge != tdec->members[member].firstEdge);
  }

  return TU_OKAY;
}

/**
 * \brief Finds terminal node in \p member.
 */

static
TU_TDEC_NODE findTerminalNode(
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new column. */
  TU_TDEC_MEMBER member         /**< member. */
)
{
  assert(tdec);
  assert(newcolumn);
  assert(member >= 0 && member < tdec->memMembers);

  

  return -1;
}

static
TU_ERROR addColumnUpdate(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced member. */
  TU_TDEC_EDGE newEdge                /**< Edge to be added appropriately. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);

#if defined(TU_DEBUG_TDEC)
  printf("      Updating reduced decomposition with component %ld, adding edge %d.\n",
    (reducedComponent - &newcolumn->reducedComponents[0]), newEdge);
#endif /* TU_DEBUG_TDEC */
  
  if (reducedComponent->terminalMember1 == reducedComponent->terminalMember2)
  {
#if defined(TU_DEBUG_TDEC)
    printf("      Unique terminal member %d is %s.\n", reducedComponent->terminalMember1,
      tdec->members[reducedComponent->terminalMember1].type == TDEC_MEMBER_TYPE_BOND ? "a bond" :
      (tdec->members[reducedComponent->terminalMember1].type == TDEC_MEMBER_TYPE_PRIME ? "prime" :
      "a polygon"));
    fflush(stdout);
#endif /* TU_DEBUG_TDEC */
    
    if (tdec->members[reducedComponent->terminalMember1].type == TDEC_MEMBER_TYPE_BOND)
    {
      /* Add edge to the bond.  */

      tdec->edges[newEdge].member = reducedComponent->terminalMember1;
      tdec->edges[newEdge].head = -1;
      tdec->edges[newEdge].tail = -1;
      TU_CALL( addEdgeToMembersEdgeList(tu, tdec, newEdge, reducedComponent->terminalMember1) );
    }
    else
    {
      assert(0 == "Adding of column with same end members not implemented.");
    }
  }
  else
  {
    assert(0 == "Adding of column with different end members not implemented.");
  }
  
  return TU_OKAY;
}

static
TU_ERROR doReorderComponent(
  TU* tu,                           /**< \ref TU environment. */
  TU_TDEC* tdec,                    /**< t-decomposition. */
  TU_TDEC_MEMBER member,            /**< Member to be processed. */
  TU_TDEC_MEMBER newParent,         /**< New parent of \p member. */
  TU_TDEC_MEMBER newMarkerToParent, /**< New marker edge linked to new parent of \p member. */
  TU_TDEC_MEMBER markerOfNewParent  /**< Counterpart to \p newMarkerToParent. */
)
{
  assert(tu);
  assert(tdec);
  assert(member >= 0);
  assert(newParent >= 0);

  TU_TDEC_MEMBER oldParent = findMemberParent(tdec, member);
  TU_TDEC_EDGE oldMarkerToParent = tdec->members[member].markerToParent;
  TU_TDEC_EDGE oldMarkerOfParent = tdec->members[member].markerOfParent;

#if defined(TU_DEBUG_TDEC)
  printf("        Flipping parenting of member %d with old parent %d and new parent %d.\n", member,
    oldParent, newParent);
  fflush(stdout);
#endif /* TU_DEBUG_TDEC */
  
  tdec->members[member].markerToParent = newMarkerToParent;
  tdec->members[member].markerOfParent = markerOfNewParent;
  tdec->edges[markerOfNewParent].childMember = member;

  if (oldMarkerToParent >= 0)
    TU_CALL( doReorderComponent(tu, tdec, oldParent, member, oldMarkerOfParent, oldMarkerToParent) );

  return TU_OKAY;
}

static
TU_ERROR reorderComponent(
  TU* tu,                 /**< \ref TU environment. */
  TU_TDEC* tdec,          /**< t-decomposition. */
  TU_TDEC_MEMBER newRoot  /**< The new root of the component. */
)
{
  assert(tu);
  assert(tdec);
  assert(newRoot >= 0 && newRoot < tdec->memMembers);
  assert(isRepresentativeMember(tdec, newRoot));

#if defined(TU_DEBUG_TDEC)
  printf("      Making member %d the new root of its component.\n", newRoot);
  fflush(stdout);
#endif /* TU_DEBUG_TDEC */


  if (tdec->members[newRoot].parentMember >= 0)
  {
    TU_CALL( doReorderComponent(tu, tdec, findMemberParent(tdec, newRoot), newRoot,
      tdec->members[newRoot].markerOfParent, tdec->members[newRoot].markerToParent) );
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

#if defined(TU_DEBUG_TDEC)
  printf("  Adding a column with %d 1's.\n", numEntries);
#endif /* TU_DEBUG_TDEC */

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

  /* Create reduced components for new edges. */
  TU_CALL( completeReducedDecomposition(tu, tdec, newcolumn, entryRows, numEntries) );
  
  /* Preprocess all reduced components individually. */

  TU_TDEC_EDGE* componentNewEdges = NULL;
  TU_CALL( TUallocStackArray(tu, &componentNewEdges, newcolumn->numReducedComponents) );

  int maxDepthComponent = -1;
  for (int i = 0; i < newcolumn->numReducedComponents; ++i)
  {
#if defined(TU_DEBUG_TDEC)
    printf("      Considering reduced component %d of depth %d.\n", i,
      newcolumn->reducedComponents[i].rootDepth);
#endif /* TU_DEBUG_TDEC */

    if (maxDepthComponent < 0 || newcolumn->reducedComponents[i].rootDepth
      > newcolumn->reducedComponents[maxDepthComponent].rootDepth)
    {
      maxDepthComponent = i;
    }

    ReducedComponent* reducedComponent = &newcolumn->reducedComponents[i];
    int numAssignedTerminals = 0;
    TU_CALL( addColumnPreprocess(tu, tdec, newcolumn, reducedComponent, reducedComponent->root,
      &numAssignedTerminals) );
    assert(numAssignedTerminals == 2);

#if defined(TU_DEBUG_TDEC)
    printf("      Preprocessing done. Reduced terminal members are %d and %d.\n",
      reducedComponent->terminalMember1, reducedComponent->terminalMember2);
    fflush(stdout);
#endif /* TU_DEBUG_TDEC */

    /* Create new edge for this component. If there is one component, this is a column edge, and
     * otherwise it is a marker edge that will be linked to a new polygon consisting of all these
     * marker edges and the column edge. */
    TU_TDEC_EDGE newEdge;
    TU_CALL( createEdge(tu, tdec, -1, &newEdge) );
    componentNewEdges[i] = newEdge;
    tdec->edges[newEdge].childMember = -1;
    tdec->edges[newEdge].head = -1;
    tdec->edges[newEdge].tail = -1;

    TU_CALL( addColumnUpdate(tu, tdec, newcolumn, reducedComponent, newEdge) );
  }

  for (int c = 0; c < newcolumn->numReducedComponents; ++c)
  {
    if (c == maxDepthComponent)
    {
#if defined(TU_DEBUG_TDEC)
      printf("      Reduced component %d has maximum depth and will remain a root.\n", c);
#endif /* TU_DEBUG_TDEC */
    }
    else
    {
      TU_CALL( reorderComponent(tu, tdec, findMember(tdec,
        tdec->edges[componentNewEdges[c]].member)) );
    }
  }
  
  if (newcolumn->numReducedComponents == 0)
  {
    assert(0 == "Adding a zero column not implemented.");
  }
  else if (newcolumn->numReducedComponents == 1)
  {
    TU_TDEC_EDGE columnEdge = componentNewEdges[0];
    tdec->edges[columnEdge].name = -1 - column;
    tdec->edges[columnEdge].childMember = -1;
  }
  else
  {
    /*
     * We create another edge for the column as well as a polygon containing all new edges and this one.
     */

    TU_TDEC_MEMBER polygon;
    TU_CALL( createMember(tu, tdec, TDEC_MEMBER_TYPE_POLYGON, &polygon) );

    TU_TDEC_EDGE columnEdge;
    TU_CALL( createEdge(tu, tdec, polygon, &columnEdge) );
    tdec->edges[columnEdge].childMember = -1;
    tdec->edges[columnEdge].head = -1;
    tdec->edges[columnEdge].tail = -1;
    tdec->edges[columnEdge].name = -1-column;
    TU_CALL( addEdgeToMembersEdgeList(tu, tdec, columnEdge, polygon) );

    for (int i = 0; i < newcolumn->numReducedComponents; ++i)
    {
      TU_TDEC_EDGE newEdge = componentNewEdges[i];
      TU_TDEC_EDGE markerEdge;
      TU_CALL( createEdge(tu, tdec, polygon, &markerEdge) );
      TU_CALL( addEdgeToMembersEdgeList(tu, tdec, markerEdge, polygon) );
      tdec->edges[markerEdge].head = -1;
      tdec->edges[markerEdge].tail = -1;

      TU_TDEC_MEMBER partnerMember = findEdgeMember(tdec, newEdge);

      if (i == maxDepthComponent)
      {
        tdec->edges[markerEdge].childMember = -1;
        tdec->edges[markerEdge].name = INT_MAX - tdec->numMarkers;
        tdec->members[polygon].parentMember = partnerMember;
        tdec->members[polygon].markerToParent = markerEdge;
        tdec->members[polygon].markerOfParent = newEdge;
        tdec->edges[newEdge].name = -INT_MAX + tdec->numMarkers;
        tdec->edges[newEdge].childMember = polygon;
      }
      else
      {
        tdec->edges[markerEdge].childMember = partnerMember;
        tdec->members[partnerMember].markerOfParent = markerEdge;
        tdec->members[partnerMember].markerToParent = newEdge;
        tdec->members[partnerMember].parentMember = polygon;
        tdec->edges[markerEdge].name = INT_MAX - tdec->numMarkers;
        tdec->edges[newEdge].name = -INT_MAX + tdec->numMarkers;
      }

      tdec->numMarkers++;
    }
  }

  TU_CALL( TUfreeStackArray(tu, &componentNewEdges) );

  newcolumn->numReducedMembers = 0;
  newcolumn->numReducedComponents = 0;

//   TU_TDEC_MEMBER newMember;
//   TU_TDEC_EDGE newEdge;
//   TU_CALL( createNewRowsPolygon(tu, tdec, &newMember, &newEdge, newcolumn->terminalNode1,
//     newcolumn->terminalNode2, column, entryRows, numEntries) );
// #if defined(TU_DEBUG_TDEC)
//   printf("    New edge is %d of %d. Terminal members are %d and %d.\n", newEdge, newMember,
//     newcolumn->terminalMember1, newcolumn->terminalMember2);
//   fflush(stdout);
// #endif /* TU_DEBUG_TDEC */


//     if (tdec->members[newcolumn->terminalMember1].type == TDEC_MEMBER_TYPE_BOND)
//     {
      /* Add edge to the bond.  */

//       tdec->edges[newEdge].member = newcolumn->terminalMember1;
//       tdec->edges[newEdge].head = newcolumn->terminalNode1;
//       tdec->edges[newEdge].tail = newcolumn->terminalNode2;
//       TU_CALL( addEdgeToMembersEdgeList(tu, tdec, newEdge, newcolumn->terminalMember1) );
//       if (newMember >= 0)
//         tdec->members[newMember].parentMember = newcolumn->terminalMember1;
//     }
//     else if (tdec->members[newcolumn->terminalMember1].type == TDEC_MEMBER_TYPE_POLYGON)
//     {
//       assert(0 == "Adding of column with same end component of type polygon not implemented.");
//     }
//     else
//     {
//       assert(tdec->members[newcolumn->terminalMember1].type == TDEC_MEMBER_TYPE_PRIME);

//       assert(0 == "Adding of column with same end component of type prime not implemented.");
//     }
//   }
//   else
//   {
    /* All modifications to the members of the path were done, so we only have to merge all along 
     * it. */

#if defined(TU_DEBUG_TDEC)
//     printf("      Creating new prime member from path between terminals.\n");
//     fflush(stdout);
#endif /* TU_DEBUG_TDEC */

//     TU_CALL( createMemberNodes(tu, tdec, newcolumn->terminalMember1) );
//     TU_TDEC_NODE terminalNode1 = findTerminalNode(tdec, newcolumn, newcolumn->terminalMember1);
//     assert(terminalNode1 >= 0);
//     TU_CALL( createMemberNodes(tu, tdec, newcolumn->terminalMember2) );
//     TU_TDEC_NODE terminalNode2 = findTerminalNode(tdec, newcolumn, newcolumn->terminalMember2);
//     assert(terminalNode2 >= 0);

//     TU_TDEC_MEMBER root = newcolumn->reducedMembers[0].member;
//     TU_TDEC_MEMBER member = newcolumn->terminalMember1;
//     while (member != root)
//     {
//       TU_TDEC_MEMBER next = findMemberParent(tdec, member);
//       TU_CALL( createMemberNodes(tu, tdec, next) );
//       member = next;
//     }
//     member = newcolumn->terminalMember2;
//     while (member != root)
//     {
//       TU_TDEC_MEMBER next = findMemberParent(tdec, member);
//       if (next != root)
//         TU_CALL( createMemberNodes(tu, tdec, next) );
//       member = next;
//     }

//     assert(0 == "Adding of column with different end components not implemented.");
//   }

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

  return TU_OKAY;
}

TU_ERROR testGraphicnessTDecomposition(TU* tu, TU_CHRMAT* matrix, TU_CHRMAT* transpose,
  bool* pisGraphic, TU_GRAPH* graph, TU_GRAPH_EDGE* basis, TU_GRAPH_EDGE* cobasis,
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
    *pisGraphic = true;
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

  TU_TDEC* tdec = NULL;
  TU_CALL( TUtdecCreate(tu, &tdec, 0, 0, 0, 0, 0) ); /* TODO: avoid reallocations. */

  /* Process each column. */
  TU_TDEC_NEWCOLUMN* newcol = NULL;
  TUtdecnewcolumnCreate(tu, &newcol);
  *pisGraphic = true;
  for (int column = 0; column < matrix->numColumns; ++column)
  {
#if defined(TU_DEBUG_TDEC)
    char name[256];
    snprintf(name, 256, "tdec-before-column-%02d.dot", column);
    FILE* dotFile = fopen(name, "w");
    bool* edgesHighlighted = NULL;
    TU_CALL( TUallocStackArray(tu, &edgesHighlighted, tdec->memEdges) );
    for (int e = 0; e < tdec->memEdges; ++e)
      edgesHighlighted[e] = false;
    for (int p = transpose->rowStarts[column]; p < transpose->rowStarts[column+1]; ++p)
    {
      int row = transpose->entryColumns[p];
      if (row < tdec->numRows && tdec->rowEdges[row].edge >= 0)
        edgesHighlighted[tdec->rowEdges[row].edge] =  true;
    }
    TU_CALL( TUtdecToDot(tu, tdec, dotFile, edgesHighlighted) );
    TU_CALL( TUfreeStackArray(tu, &edgesHighlighted) );
    fclose(dotFile);
#endif /* TU_DEBUG_TDEC */

    TU_CALL( TUtdecAddColumnCheck(tu, tdec, newcol,
      &transpose->entryColumns[transpose->rowStarts[column]],
      transpose->rowStarts[column+1] - transpose->rowStarts[column]) );

    if (newcol->remainsGraphic)
    {
      TUtdecAddColumnApply(tu, tdec, newcol, column,
        &transpose->entryColumns[transpose->rowStarts[column]],
        transpose->rowStarts[column+1] - transpose->rowStarts[column]);
    }
    else
    {
      *pisGraphic = false;
      assert(!"Not implemented");
    }
  }
  TUtdecnewcolumnFree(tu, &newcol);

#if defined(TU_DEBUG_TDEC)
  char name[256];
  snprintf(name, 256, "tdec-final.dot");
  FILE* dotFile = fopen(name, "w");
  bool* edgesHighlighted = NULL;
  TU_CALL( TUallocStackArray(tu, &edgesHighlighted, tdec->memEdges) );
  for (int e = 0; e < tdec->memEdges; ++e)
    edgesHighlighted[e] = false;
  TU_CALL( TUtdecToDot(tu, tdec, dotFile, edgesHighlighted) );
  TU_CALL( TUfreeStackArray(tu, &edgesHighlighted) );
  fclose(dotFile);
#endif /* TU_DEBUG_TDEC */

  if (*pisGraphic && graph)
  {
    TUtdecToGraph(tu, tdec, graph, true, basis, cobasis, NULL);
  }

  TUtdecFree(tu, &tdec);

  return TU_OKAY;
}

