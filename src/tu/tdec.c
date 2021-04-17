#define TU_DEBUG /* Uncomment to enable general debugging. */
// #define TU_DEBUG_SPLITTING /* Uncomment to enable debug output for splitting of polygons. */
// #define TU_DEBUG_DOT /* Uncomment to output dot files after modifications of the t-decomposition. */

// TODO: Refactor replacement of an edge by another one.
// TODO: Refactor creation of a pair of marker edges instead of one.

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

#define SWAP_INTS(a, b) \
  do \
  { \
    int tmp = a; \
    a = b; \
    b = tmp; \
  } \
  while (false)

static inline
bool isRepresentativeMember(
  TU_TDEC* tdec,        /**< t-decomposition. */
  TU_TDEC_MEMBER member /**< Member of \p tdec. */
)
{
  assert(tdec);
  return tdec->members[member].representativeMember < 0;
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
        if (findEdgeMember(tdec, edge) != member)
          return TUconsistencyMessage("edge %d belongs to member %d but is in member %d's edge list.", edge,
            findEdgeMember(tdec, edge), member);
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
      && tdec->members[member].type != TDEC_MEMBER_TYPE_POLYGON
      && tdec->members[member].type != TDEC_MEMBER_TYPE_LOOP)
    {
      return TUconsistencyMessage("member %d has invalid type", member);
    }
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
    if (!isRepresentativeMember(tdec, member))
      continue;

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
    if (!isRepresentativeMember(tdec, member))
      continue;

    int length = 0;
    TU_TDEC_MEMBER current;
    for (current = tdec->members[member].parentMember; current >= 0; current = tdec->members[current].parentMember)
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
  TYPE_1_CLOSES_CYCLE = 1,        /**< Parent marker edge plus path is a cycle. */
  TYPE_2_SHORTCUT = 2,            /**< One node of the parent marker edge is a path end and the other is an inner node. */
  TYPE_3_EXTENSION = 3,           /**< One node of the parent marker edge is a path end and the other does not belong to the path. */
  TYPE_4_CONNECTS_TWO_PATHS = 4,  /**< Adding the parent marker edge connects two paths to a single one. */
  TYPE_5_ROOT = 5,                /**< Root member. */
  TYPE_6_OTHER = 6                /**< All other cases. */
} Type;

/**
 * \brief Additional edge information specific to a path.
 */

typedef struct _PathEdge
{
  TU_TDEC_EDGE edge;          /**< \brief The edge in the t-decomposition. */
  struct _PathEdge* next;  /**< \brief Next edge of this reduced member, or \c NULL. */
} PathEdge;

/**
 * \brief Additional member information specfic to a given path.
 *
 * @TODO: Maybe add parent reduced member as well, so we don't have to go via the membersToReducedMembers array.
 */

typedef struct _ReducedMember
{
  TU_TDEC_MEMBER member;            /**< \brief The member from the t-decomposition. */
  TU_TDEC_MEMBER rootMember;        /**< \brief The root member of this component of the t-decomposition. */
  int depth;                        /**< \brief Depth of this member in the reduced t-decomposition. */
  Type type;                        /**< \brief Type of this member. */
  struct _ReducedMember* parent;    /**< \brief Parent in the reduced t-decomposition. */
  int numChildren;                  /**< \brief Number of children in the reduced t-decomposition. */
  struct _ReducedMember** children; /**< \brief Children in the reduced t-decomposition. */
  PathEdge* firstPathEdge;          /**< \brief First edge in linked list of path edges of \p member. */
  TU_TDEC_NODE primeEndNodes[4];    /**< \brief For primes, the end nodes of the paths inside the member (or -1). */
} ReducedMember;

typedef struct _ReducedComponent
{
  int rootDepth;                    /**< \brief Depth of reduced root member. */
  ReducedMember* root;              /**< \brief Reduced root member. */
  TU_TDEC_NODE terminalNode[2];     /**< \brief Terminal nodes of path. */
  TU_TDEC_MEMBER terminalMember[2]; /**< \brief Terminal members of path. */
  int numTerminals;
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

  PathEdge* pathEdges;                      /**< \brief Storage for edge lists of path edges. */
  int memPathEdges;                         /**< \brief Allocated memory for \c pathEdges. */
  int numPathEdges;                         /**< \brief Number of stored edges in \c pathEdges. */

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
    if (!isRepresentativeMember(tdec, member))
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

        if (findMember(tdec, tdec->members[findMember(tdec, tdec->edges[edge].childMember)].parentMember) != findMember(tdec, member))
        {
          TUfreeStackArray(tu, &countChildren);
          return TUconsistencyMessage("member %d has child edge %d for child %d whose parent member is %d",
            member, edge, findMember(tdec, tdec->edges[edge].childMember),
            findMember(tdec, tdec->members[findMember(tdec, tdec->edges[edge].childMember)].parentMember));
        }
        if (tdec->members[findMember(tdec, tdec->edges[edge].childMember)].markerOfParent != edge)
        {
          TUfreeStackArray(tu, &countChildren);
          return TUconsistencyMessage("member %d has child edge %d for child %d whose parent's markerOfParent is %d",
            member, edge, findMember(tdec, tdec->edges[edge].childMember),
            tdec->members[findMember(tdec, tdec->edges[edge].childMember)].markerOfParent);
        }
        TU_TDEC_EDGE markerChild = tdec->members[findMember(tdec, tdec->edges[edge].childMember)].markerToParent;
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
  assert(tdec);
  assert(edge >= 0);
  assert(edge < tdec->memEdges);
  assert(tdec->edges[edge].head >= 0);
  assert(tdec->edges[edge].head < tdec->memNodes);
  return findNode(tdec, tdec->edges[edge].head);
}

static // TODO: inline
TU_TDEC_NODE findEdgeTail(TU_TDEC* tdec, TU_TDEC_EDGE edge)
{
  assert(tdec);
  assert(edge >= 0);
  assert(edge < tdec->memEdges);
  assert(tdec->edges[edge].tail >= 0);
  assert(tdec->edges[edge].tail < tdec->memNodes);
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
    TUdbgMsg(10, "createNode returns free node %d.\n", node);
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
    TUdbgMsg(12, "createNode enlarges node array to %d and returns node %d.\n", newSize, node);
  }
  tdec->nodes[node].representativeNode = -1;
  tdec->numNodes++;

  *pnode = node;

  return TU_OKAY;
}

/**
 * \brief Adds \p edge to the edge list of \p member.
 *
 * Requires that \p edge already has \p member as its member.
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
  assert(tdec->edges[edge].member == member);

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
 * \brief Removes \p edge from the edge list of \p member.
 */

static
TU_ERROR removeEdgeFromMembersEdgeList(
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
  assert(findMember(tdec, tdec->edges[edge].member) == member);

  if (tdec->members[member].numEdges == 1)
    tdec->members[member].firstEdge = -1;
  else
  {
    if (tdec->members[member].firstEdge == edge)
      tdec->members[member].firstEdge = tdec->edges[edge].next;

    assert(tdec->members[member].firstEdge != edge);

    tdec->edges[tdec->edges[edge].prev].next = tdec->edges[edge].next;
    tdec->edges[tdec->edges[edge].next].prev = tdec->edges[edge].prev;
  }

  tdec->members[member].numEdges--;

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
    TUdbgMsg(12, "Creating edge %d by using a free edge.\n", edge);
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
    TUdbgMsg(12, "Creating edge %d and reallocating edge to array %d elements.\n", edge, newSize);
  }

  tdec->edges[edge].tail = -1;
  tdec->edges[edge].head = -1;
  tdec->edges[edge].name = INT_MAX/2;
  tdec->edges[edge].member = member;
  tdec->numEdges++;

  *pedge = edge;

  return TU_OKAY;
}

static
TU_ERROR createMarkerEdge(
  TU* tu,                 /**< \ref TU environment. */
  TU_TDEC* tdec,          /**< t-decomposition. */
  TU_TDEC_EDGE* pedge,    /**< Pointer for storing the new edge. */
  TU_TDEC_MEMBER member,  /**< Member this edge belongs to. */
  TU_TDEC_NODE tail,      /**< Tail node of this edge. */
  TU_TDEC_NODE head,      /**< Head node of this edge. */
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
  TUdbgMsg(12, "Created %s marker edge {%d,%d} of member %d.\n", isParent ? "parent" : "child", head, tail, member);

  return TU_OKAY;
}

static
TU_ERROR createMember(
  TU* tu,                   /**< \ref TU environment. */
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

  TUdbgMsg(10, "Creating %s member %d.\n",
    type == TDEC_MEMBER_TYPE_BOND ? "bond" : (type == TDEC_MEMBER_TYPE_PRIME ? "prime" : "polygon"), *pmember);

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

  TUdbgMsg(0, "TUtdecToGraph for t-decomposition.\n");

  TU_CALL( TUgraphClear(tu, graph) );

  TU_GRAPH_EDGE* localEdgeElements = NULL;
  if (edgeElements)
    localEdgeElements = edgeElements;
  else if (basis || cobasis)
    TU_CALL( TUallocStackArray(tu, &localEdgeElements, tdec->memEdges) );
  TU_GRAPH_NODE* tdecNodesToGraphNodes = NULL;
  TU_CALL( TUallocStackArray(tu, &tdecNodesToGraphNodes, tdec->memNodes) );
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
    TUdbgMsg(2, "Member %d is %s with %d edges.\n", member,
      type == TDEC_MEMBER_TYPE_BOND ? "a bond" : (type == TDEC_MEMBER_TYPE_POLYGON ? "a polygon" : "prime"),
       tdec->members[member].numEdges);

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
    else if (type == TDEC_MEMBER_TYPE_POLYGON)
    {
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
    else
    {
      assert(type == TDEC_MEMBER_TYPE_LOOP);

      TU_GRAPH_NODE v;
      TU_CALL( TUgraphAddNode(tu, graph, &v) );
      TU_CALL( TUgraphAddEdge(tu, graph, v, v, &graphEdge) );
      tdecEdgesToGraphEdges[edge] = graphEdge;
      if (localEdgeElements)
        localEdgeElements[graphEdge] = tdec->edges[edge].name;
    }
  }

  /* Merge respective parent and child edges. */

  if (merge)
  {
    TUdbgMsg(2, "Before merging, the graph has %d nodes and %d edges.\n", TUgraphNumNodes(graph),
      TUgraphNumEdges(graph));

    for (int m = 0; m < tdec->numMembers; ++m)
    {
      if (!isRepresentativeMember(tdec, m) || tdec->members[m].parentMember < 0)
        continue;

      TU_GRAPH_EDGE parent = tdecEdgesToGraphEdges[tdec->members[m].markerOfParent];
      TU_GRAPH_EDGE child = tdecEdgesToGraphEdges[tdec->members[m].markerToParent];
      TU_GRAPH_NODE parentU = TUgraphEdgeU(graph, parent);
      TU_GRAPH_NODE parentV = TUgraphEdgeV(graph, parent);
      TU_GRAPH_NODE childU = TUgraphEdgeU(graph, child);
      TU_GRAPH_NODE childV = TUgraphEdgeV(graph, child);

      TUdbgMsg(2, "Merging edges %d = {%d,%d} <%d> and %d = {%d,%d} <%d>.\n", parent, parentU, parentV,
        tdec->edges[tdec->members[m].markerOfParent].name, child, childU, childV,
        tdec->edges[tdec->members[m].markerToParent].name);

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

      TUdbgMsg(2, "Graph edge %d = {%d,%d} <%d>\n", e, TUgraphEdgeU(graph, e), TUgraphEdgeV(graph, e),
        localEdgeElements[e]);

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
  FILE* stream,           /**< File stream. */
  TU_TDEC* tdec,          /**< t-decomposition. */
  TU_TDEC_MEMBER member,  /**< Member this edge belongs to. */
  TU_TDEC_EDGE edge,      /**< Edge. */
  int u,                  /**< First node. */
  int v,                  /**< Second node. */
  bool red                /**< Whether to color it red. */
)
{
  assert(stream);
  assert(member >= 0);
  assert(edge >= 0);

  member = findMember(tdec, member);

  char type = (tdec->members[member].type == TDEC_MEMBER_TYPE_BOND) ?
    'P' : (tdec->members[member].type == TDEC_MEMBER_TYPE_POLYGON ? 'S' : 'R');
  const char* redStyle = red ? ",color=red" : "";
  if (tdec->members[member].markerToParent == edge)
  {
    fprintf(stream, "    %c_%d_%d -- %c_p_%d [label=\"%d\",style=dashed%s];\n", type, member, u, type, member, edge, redStyle);
    fprintf(stream, "    %c_p_%d -- %c_%d_%d [label=\"%d\",style=dashed%s];\n", type, member, type, member, v, edge, redStyle);
    fprintf(stream, "    %c_%d_%d [shape=box];\n", type, member, u);
    fprintf(stream, "    %c_%d_%d [shape=box];\n", type, member, v);
    fprintf(stream, "    %c_p_%d [style=dashed];\n", type, member);
  }
  else if (tdec->edges[edge].childMember >= 0)
  {
    TU_TDEC_MEMBER child = findMember(tdec, tdec->edges[edge].childMember);
    char childType = (tdec->members[child].type == TDEC_MEMBER_TYPE_BOND) ?
      'P' : (tdec->members[child].type == TDEC_MEMBER_TYPE_POLYGON ? 'S' : 'R');
    fprintf(stream, "    %c_%d_%d -- %c_c_%d [label=\"%d\",style=dotted%s];\n", type, member, u, type, child, edge, redStyle);
    fprintf(stream, "    %c_c_%d -- %c_%d_%d [label=\"%d\",style=dotted%s];\n", type, child, type, member, v, edge, redStyle);
    fprintf(stream, "    %c_%d_%d [shape=box];\n", type, member, u);
    fprintf(stream, "    %c_%d_%d [shape=box];\n", type, member, v);
    fprintf(stream, "    %c_c_%d [style=dotted];\n", type, child);

    fprintf(stream, "    %c_p_%d -- %c_c_%d [style=dashed,dir=forward];\n", childType, child, type, child);
  }
  else
  {
    fflush(stdout);
    fprintf(stream, "    %c_%d_%d -- %c_%d_%d [label=\"%d <%d>\",style=bold%s];\n", type, member, u, type, member, v,
      edge, tdec->edges[edge].name, redStyle);
    fprintf(stream, "    %c_%d_%d [shape=box];\n", type, member, u);
    fprintf(stream, "    %c_%d_%d [shape=box];\n", type, member, v);
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
    if (!isRepresentativeMember(tdec, member))
      continue;

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
    else if (tdec->members[member].type == TDEC_MEMBER_TYPE_POLYGON)
    {
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
    else
    {
      assert(tdec->members[member].type == TDEC_MEMBER_TYPE_LOOP);

      edgeToDot(stream, tdec, member, tdec->members[member].firstEdge, 0, 0, false);
    }
    fprintf(stream, "  }\n");
  }
  fprintf(stream, "}\n");

  return TU_OKAY;
}

static inline
void debugDot(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC* tdec,                /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn  /**< new column. */
)
{
#if defined(TU_DEBUG_DOT)
  static int dotFileCounter = 1;
  char name[256];
  snprintf(name, 256, "tdec-%03d.dot", dotFileCounter);
  TUdbgMsg(0, "Writing <%s>...", name);
  FILE* dotFile = fopen(name, "w");
  TU_CALL( TUtdecToDot(tu, tdec, dotFile, newcolumn ? newcolumn->edgesInPath : NULL) );
  fclose(dotFile);
  TUdbgMsg(0, " done.\n");

  dotFileCounter++;
#endif /* TU_DEBUG_DOT */
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

  newcolumn->pathEdges = NULL;
  newcolumn->memPathEdges = 0;
  newcolumn->numPathEdges = 0;

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
  if (newcolumn->pathEdges)
    TU_CALL( TUfreeBlockArray(tu, &newcolumn->pathEdges) );
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
  newcolumn->numPathEdges = 0;

  // TODO: Remember sizes of these arrays.

  /* memEdges does not suffice since new edges can be created by squeezing off.
   * Each squeezing off introduces 4 new edges, and we might apply this twice for each polygon member. */
  TU_CALL( TUreallocBlockArray(tu, &newcolumn->edgesInPath, tdec->memEdges + 8*tdec->numMembers + 32) );
  TU_CALL( TUreallocBlockArray(tu, &newcolumn->nodesDegree, tdec->memNodes) );

#if defined(TU_DEBUG_DOT)
  for (int e = 0; e < tdec->memEdges + 8*tdec->numMembers + 32; ++e)
    newcolumn->edgesInPath[e] = false;
#endif /* TU_DEBUG_DOT */

  return TU_OKAY;
}

/**
 * \brief Ensures that the child marker in \p member for \p childMember is not parallel to the parent marker of
 * \p member.
 * 
 * Note that this can only happen when a component is reordered (receiving a new root) because reduced components shall
 * be joined.
 */

static
TU_ERROR ensureChildParentMarkersNotParallel(
  TU* tu,                     /**< \ref TU environment. */
  TU_TDEC* tdec,              /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new column. */
  TU_TDEC_MEMBER childMember, /**< Child member of \p member. */
  TU_TDEC_MEMBER member       /**< Parent of \p childMember. */
)
{
  assert(tu);
  assert(tdec);
  assert(childMember >= 0);
  assert(member == findMemberParent(tdec, childMember));

  TU_TDEC_MEMBER parentMember = findMemberParent(tdec, member);
  if (tdec->members[member].type != TDEC_MEMBER_TYPE_PRIME || parentMember < 0)
    return TU_OKAY;

  TUdbgMsg(10, "Checking if the child marker of %d for child member %d is not parallel to %d's parent marker.\n",
    member, childMember, member);

  TU_TDEC_EDGE childMarkerEdge = tdec->members[childMember].markerOfParent;
  TU_TDEC_NODE nodes[4] = {
    findEdgeTail(tdec, childMarkerEdge),
    findEdgeHead(tdec, childMarkerEdge),
    findEdgeTail(tdec, tdec->members[member].markerToParent),
    findEdgeHead(tdec, tdec->members[member].markerToParent),
  };

  if ((nodes[0] == nodes[2] && nodes[1] == nodes[3]) || (nodes[0] == nodes[3] && nodes[1] == nodes[2]))
  {
    TUdbgMsg(12, "Child marker edge %d is parallel to parent marker edge %d.\n",
      tdec->members[childMember].markerOfParent, tdec->members[member].markerToParent);

    if (tdec->members[parentMember].type != TDEC_MEMBER_TYPE_BOND)
    {
      TUdbgMsg(12, "Parent is not a bond, so we create one between edges %d and %d.\n",
        tdec->members[member].markerToParent, tdec->members[member].markerOfParent);
      
      TU_TDEC_MEMBER newBond;
      TU_CALL( createMember(tu, tdec, TDEC_MEMBER_TYPE_BOND, &newBond) );
      TU_TDEC_EDGE markerOfParent = tdec->members[member].markerOfParent;
      TU_TDEC_EDGE newMarkerOfParent;
      TU_TDEC_EDGE newMarkerToParent;
      TU_CALL( createMarkerEdge(tu, tdec, &newMarkerOfParent, parentMember, tdec->edges[markerOfParent].tail,
        tdec->edges[markerOfParent].head, false) );
      TU_CALL( createMarkerEdge(tu, tdec, &newMarkerToParent, newBond, -1, -1, true) );
      tdec->numMarkers++;

      /* Replace markerOfParent by newMarkerOfParent. */
      tdec->edges[newMarkerOfParent].next = tdec->edges[markerOfParent].next;
      tdec->edges[newMarkerOfParent].prev = tdec->edges[markerOfParent].prev;
      tdec->edges[tdec->edges[markerOfParent].next].prev = newMarkerOfParent;
      tdec->edges[tdec->edges[markerOfParent].prev].next = newMarkerOfParent;
      if (tdec->members[parentMember].firstEdge == markerOfParent)
        tdec->members[parentMember].firstEdge = newMarkerOfParent;
      tdec->edges[newMarkerOfParent].childMember = newBond;
      tdec->edges[newMarkerOfParent].tail = tdec->edges[markerOfParent].tail;
      tdec->edges[newMarkerOfParent].head = tdec->edges[markerOfParent].head;
      tdec->edges[markerOfParent].childMember = member;
      tdec->edges[markerOfParent].member = newBond;
      tdec->edges[markerOfParent].tail = -1;
      tdec->edges[markerOfParent].head = -1;
      TU_CALL( addEdgeToMembersEdgeList(tu, tdec, markerOfParent, newBond) );
      TU_CALL( addEdgeToMembersEdgeList(tu, tdec, newMarkerToParent, newBond) );
      tdec->members[newBond].markerOfParent = newMarkerOfParent;
      tdec->members[newBond].markerToParent = newMarkerToParent;
      tdec->members[newBond].parentMember = parentMember;
      tdec->members[member].parentMember = newBond;
      parentMember = newBond;
      newcolumn->membersToReducedMembers[newBond] = NULL;
    }

    assert(tdec->members[parentMember].type == TDEC_MEMBER_TYPE_BOND);

    TUdbgMsg(12, "Removing child marker edge %d from member %d and adding it to the new bond %d.\n", childMarkerEdge,
      member, parentMember);

    TU_CALL( removeEdgeFromMembersEdgeList(tu, tdec, childMarkerEdge, member) );
    tdec->edges[childMarkerEdge].member = parentMember;
    TU_CALL( addEdgeToMembersEdgeList(tu, tdec, childMarkerEdge, parentMember) );
    tdec->edges[childMarkerEdge].tail = -1;
    tdec->edges[childMarkerEdge].head = -1;
    tdec->members[childMember].parentMember = parentMember;

    debugDot(tu, tdec, NULL);
  }
  
  return TU_OKAY;
}

/**
 * \brief Creates, if necessary, the reduced member for \p member and calls itself for the parent.
 */

static
TU_ERROR createReducedMembers(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new column. */
  TU_TDEC_MEMBER member,              /**< Member to create reduced member for. */
  ReducedMember** rootDepthMinimizer, /**< Array mapping root members to the depth minimizer. */
  ReducedMember** pReducedMember      /**< Pointer for storing the created reduced member. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(member >= 0);
  assert(pReducedMember);

  TUdbgMsg(8, "Attempting to create reduced member %d.\n", member);

  ReducedMember* reducedMember = newcolumn->membersToReducedMembers[member];
  if (reducedMember)
  {
    /* This member is a known reduced member. If we meet an existing path of low depth, we remember
     * that. */

    TUdbgMsg(10, "Reduced member exists.\n");

    if (!rootDepthMinimizer[reducedMember->rootMember] ||
      reducedMember->depth < rootDepthMinimizer[reducedMember->rootMember]->depth)
    {
      TUdbgMsg(8, "Updating depth to %d.\n", reducedMember->depth);
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
    reducedMember->primeEndNodes[0] = -1;
    reducedMember->primeEndNodes[1] = -1;
    reducedMember->primeEndNodes[2] = -1;
    reducedMember->primeEndNodes[3] = -1;

    TUdbgMsg(10, "Reduced member is new.\n");

    TU_TDEC_MEMBER parentMember = findMemberParent(tdec, member);
    if (parentMember >= 0)
    {
      TU_CALL( ensureChildParentMarkersNotParallel(tu, tdec, newcolumn, member, parentMember) );
      /* The parent member might have changed. */
      parentMember = findMemberParent(tdec, member);

      ReducedMember* parentReducedMember;
      TU_CALL( createReducedMembers(tu, tdec, newcolumn, parentMember, rootDepthMinimizer, &parentReducedMember) );
      reducedMember->parent = parentReducedMember;
      reducedMember->depth = parentReducedMember->depth + 1;
      reducedMember->rootMember = parentReducedMember->rootMember;
      parentReducedMember->numChildren++;
    }
    else
    {
      reducedMember->parent = NULL;
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
      TUdbgMsg(8, "Initializing the new reduced component %d.\n", newcolumn->numReducedComponents);

      newcolumn->reducedComponents[newcolumn->numReducedComponents].root = reducedMember;
      newcolumn->numReducedComponents++;
    }

    TUdbgMsg(8, "The root member of %d is %d.\n", reducedMember->member, reducedMember->rootMember);
  }
  
  *pReducedMember = reducedMember;

  return TU_OKAY;
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
  TUdbgMsg(4, "Computing reduced t-decomposition.\n");

  /* Enlarge members array. */
  int maxRow = 0;
  for (int p = 0; p < numEntries; ++p)
  {
    if (entryRows[p] > maxRow)
      maxRow = entryRows[p];
  }
  if (newcolumn->memReducedMembers < tdec->numMembers + numEntries)
  {
    newcolumn->memReducedMembers = tdec->memMembers + maxRow + 1;
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->reducedMembers, newcolumn->memReducedMembers) );
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->membersToReducedMembers,
      newcolumn->memReducedMembers) );
  }

  /* Initialize the mapping from members to reduced members. */
  for (int m = 0; m < tdec->numMembers; ++m)
    newcolumn->membersToReducedMembers[m] = NULL;

  ReducedMember** rootDepthMinimizer = NULL;
  /* Factor 2 because of possible new bonds due to ensureChildParentMarkersNotParallel */
  TU_CALL( TUallocStackArray(tu, &rootDepthMinimizer, 2*tdec->numMembers) ); 
  for (int m = 0; m < 2*tdec->numMembers; ++m)
    rootDepthMinimizer[m] = NULL;
  newcolumn->numReducedMembers = 0;
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    TU_TDEC_EDGE edge = (row < tdec->numRows) ? tdec->rowEdges[row].edge : -1;
    TUdbgMsg(6, "Entry %d is row %d of %d and corresponds to edge %d.\n", p, row, tdec->numRows, edge);
    if (edge >= 0)
    {
      TU_TDEC_MEMBER member = findEdgeMember(tdec, edge);
      TUdbgMsg(8, "Edge %d exists and belongs to member %d.\n", edge, member);
      ReducedMember* reducedMember;
      TU_CALL( createReducedMembers(tu, tdec, newcolumn, member, rootDepthMinimizer, &reducedMember) );

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
    TUdbgMsg(6, "Considering reduced component %d with initial root member %d.\n", i,
      newcolumn->reducedComponents[i].root->member);
    TUdbgMsg(8, "The minimizer is %d.\n",
      rootDepthMinimizer[newcolumn->reducedComponents[i].root->member]->member);

    newcolumn->reducedComponents[i].rootDepth =
      rootDepthMinimizer[newcolumn->reducedComponents[i].root->member]->depth;
    newcolumn->reducedComponents[i].root =
      rootDepthMinimizer[newcolumn->reducedComponents[i].root->member];
    newcolumn->reducedComponents[i].numTerminals = 0;
    TUdbgMsg(8, "Member %d is the new root of the reduced decomposition of this component.\n",
      newcolumn->reducedComponents[i].root->member);
  }

  /* Allocate memory for children. */
  if (newcolumn->memChildrenStorage < newcolumn->numReducedMembers)
  {
    newcolumn->memChildrenStorage = 2*newcolumn->numReducedMembers;
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->childrenStorage, newcolumn->memChildrenStorage) );
  }
  newcolumn->usedChildrenStorage = 0;

  /* Set memory pointer of each reduced member. */
  for (int m = 0; m < newcolumn->numReducedMembers; ++m)
  {
    ReducedMember* reducedMember = &newcolumn->reducedMembers[m];
    if (reducedMember->depth >= rootDepthMinimizer[reducedMember->rootMember]->depth)
    {
      TUdbgMsg(6, "Member %d's depth is greater than that of the root, and it has %d children.\n",
        reducedMember->member, reducedMember->numChildren);
      reducedMember->children = &newcolumn->childrenStorage[newcolumn->usedChildrenStorage];
      newcolumn->usedChildrenStorage += reducedMember->numChildren;
      reducedMember->numChildren = 0;
    }
    else
    {
      TUdbgMsg(6, "Member %d's depth is smaller than that of its new root.\n", reducedMember->member);
      continue;
    }
  }

  TUdbgMsg(4, "Total number of children is %d / %d.\n", newcolumn->usedChildrenStorage, newcolumn->memChildrenStorage);

  /* Set children of each reduced member. */
  for (int m = 0; m < newcolumn->numReducedMembers; ++m)
  {
    ReducedMember* reducedMember = &newcolumn->reducedMembers[m];
    if (reducedMember->depth <= rootDepthMinimizer[reducedMember->rootMember]->depth)
    {
      TUdbgMsg(6, "Member %d's depth is smaller than or equal to that of its reduced root.\n", reducedMember->member);
      continue;
    }

    TU_TDEC_MEMBER parentMember = findMemberParent(tdec, newcolumn->reducedMembers[m].member);
    ReducedMember* parentReducedMember = parentMember >= 0 ? newcolumn->membersToReducedMembers[parentMember] : NULL;
    TUdbgMsg(6, "Member %d's depth is greater than that of its reduced root. Its parent is %d, and reduced parent %p.\n",
      reducedMember->member, parentMember, parentReducedMember);

    if (parentReducedMember)
    {
      TUdbgMsg(6, "Reduced member %ld (= member %d) has %d (= member %d) as child %d.\n",
        (parentReducedMember - newcolumn->reducedMembers),
        parentReducedMember->member, m, newcolumn->reducedMembers[m].member, parentReducedMember->numChildren);
      parentReducedMember->children[parentReducedMember->numChildren] = &newcolumn->reducedMembers[m];
      parentReducedMember->numChildren++;
    }
  }

  TU_CALL( TUfreeStackArray(tu, &rootDepthMinimizer) );

  return TU_OKAY;
}

/**
 * Allocates a path edge structure.
 */

static
TU_ERROR createPathEdge(
  TU* tu,                       /**< \ref TU environment. */
  TU_TDEC_NEWCOLUMN* newcolumn, /**< new column. */
  TU_TDEC_EDGE edge,            /**< The edge it refers to. */
  PathEdge* next,               /**< The next edge in the singly linked list of this reduced member. */
  PathEdge** pNewPathEdge       /**< Pointer for storing the new path edge. */
)
{
  assert(tu);
  assert(newcolumn);

  if (newcolumn->numPathEdges == newcolumn->memPathEdges)
  {
    newcolumn->memPathEdges = 16 + 2 * newcolumn->memPathEdges;
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->pathEdges, newcolumn->memPathEdges) );
  }

  *pNewPathEdge = &newcolumn->pathEdges[newcolumn->numPathEdges];
  newcolumn->pathEdges[newcolumn->numPathEdges].edge = edge;
  newcolumn->pathEdges[newcolumn->numPathEdges].next = next;
  newcolumn->numPathEdges++;

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

  TUdbgMsg(4, "Completing reduced decomposition: increasing #rows from %d to %d.\n", tdec->numRows, newNumRows);

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

      TUdbgMsg(8, "New row %d is edge %d of member %d.\n", r, edge, member);

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

      TUdbgMsg(4, "Creating reduced member for edge %d of member %d.\n", edge, member);

      assert(newcolumn->numReducedMembers < newcolumn->memReducedMembers);
      ReducedMember* reducedMember = &newcolumn->reducedMembers[newcolumn->numReducedMembers];
      newcolumn->numReducedMembers++;
      reducedMember->numChildren = 0;
      reducedMember->member = member;
      reducedMember->depth = 0;
      reducedMember->rootMember = -1;
      reducedMember->type = TYPE_5_ROOT;

      PathEdge* reducedEdge;
      TU_CALL( createPathEdge(tu, newcolumn, edge, NULL, &reducedEdge) );
      reducedMember->firstPathEdge = reducedEdge;

      if (newcolumn->numReducedComponents == newcolumn->memReducedComponents)
      {
        newcolumn->memReducedComponents = 2 * newcolumn->memReducedComponents + 16;
        TU_CALL( TUreallocBlockArray(tu, &newcolumn->reducedComponents,
          newcolumn->memReducedComponents) );
      }

      newcolumn->membersToReducedMembers[member] = reducedMember;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].root = reducedMember;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].rootDepth = 0;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].numTerminals = 0;
#if !defined(NDEBUG)
      newcolumn->reducedComponents[newcolumn->numReducedComponents].terminalMember[0] = INT_MIN;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].terminalMember[1] = INT_MIN;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].terminalNode[0] = INT_MIN;
      newcolumn->reducedComponents[newcolumn->numReducedComponents].terminalNode[1] = INT_MIN;
#endif /* !NDEBUG */
      newcolumn->numReducedComponents++;
    }
  }

  tdec->numRows = newNumRows;

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

  TUdbgMsg(4, "Initializing edge lists for members of reduced t-decomposition.\n");

  for (int v = 0; v < tdec->memNodes; ++v)
    newcolumn->nodesDegree[v] = 0;
  for (int e = 0; e < tdec->memEdges; ++e)
    newcolumn->edgesInPath[e] = false;

  assert(newcolumn->numPathEdges == 0);

  /* Start with empty lists. */
  for (int i = 0; i < newcolumn->numReducedMembers; ++i)
    newcolumn->reducedMembers[i].firstPathEdge = NULL;

  /* Fill edge lists. */
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    TU_TDEC_EDGE edge = (row < tdec->numRows) ? tdec->rowEdges[row].edge : -1;
    if (edge >= 0)
    {
      TU_TDEC_MEMBER member = findEdgeMember(tdec, edge);
      assert(member >= 0);
      ReducedMember* reducedMember = newcolumn->membersToReducedMembers[member];
      assert(reducedMember);
      PathEdge* pathEdge;
      TU_CALL( createPathEdge(tu, newcolumn, edge, reducedMember->firstPathEdge, &pathEdge) );
      reducedMember->firstPathEdge = pathEdge;
      newcolumn->edgesInPath[edge] = true;
      if (tdec->members[member].type == TDEC_MEMBER_TYPE_PRIME)
      {
        newcolumn->nodesDegree[findEdgeHead(tdec, edge)]++;
        newcolumn->nodesDegree[findEdgeTail(tdec, edge)]++;
      }

      TUdbgMsg(6, "Edge %d <%d> belongs to reduced member %ld which is member %d.\n", edge, tdec->edges[edge].name,
        (reducedMember - newcolumn->reducedMembers), reducedMember->member);
    }
  }

  return TU_OKAY;
}

/**
 * \brief Count the number of children of a reduced member having certain types.
 */

static
TU_ERROR countChildrenTypes(
  TU* tu,                           /**< \ref TU environment. */
  TU_TDEC* tdec,                    /**< t-decomposition. */
  ReducedMember* reducedMember,     /**< Reduced member. */
  int* pNumOneEnd,                  /**< Number of children that (recursively) must contain one path end. */
  int* pNumTwoEnds,                 /**< Number of children that (recursively) must contain two path ends. */
  TU_TDEC_EDGE childMarkerEdges[2]  /**< Array for storing a child marker edges containing one/two path ends. */
)
{
  assert(tu);
  assert(tdec);
  assert(reducedMember);
  assert(pNumOneEnd);
  assert(pNumTwoEnds);
  assert(childMarkerEdges);

  *pNumOneEnd = 0;
  *pNumTwoEnds = 0;
  childMarkerEdges[0] = -1;
  childMarkerEdges[1] = -1;
  int nextChildMarker = 0;

  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    ReducedMember* child = reducedMember->children[c];
    assert(child);

    if (child->type == TYPE_2_SHORTCUT || child->type == TYPE_3_EXTENSION)
    {
      if (nextChildMarker < 2)
      {
        childMarkerEdges[nextChildMarker] = tdec->members[findMember(tdec, child->member)].markerOfParent;
        nextChildMarker++;
      }
      (*pNumOneEnd)++;
    }
    else if (child->type == TYPE_4_CONNECTS_TWO_PATHS)
    {
      if (nextChildMarker < 2)
      {
        childMarkerEdges[nextChildMarker] = tdec->members[findMember(tdec, child->member)].markerOfParent;
        nextChildMarker++;
      }
      (*pNumTwoEnds)++;
    }
  }

  return TU_OKAY;
}

static
TU_ERROR determineTypeBond(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int numOneEnd,                      /**< Number of child markers containing one path end. */
  int numTwoEnds,                     /**< Number of child markers containing two path ends. */
  TU_TDEC_EDGE childMarkerEdges[2],   /**< Marker edges of children containing one/two path ends. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedComponent);
  assert(reducedMember);
  assert(numOneEnd >= 0);
  assert(numTwoEnds >= 0);
  assert(numOneEnd + 2*numTwoEnds <= 2);
  assert(childMarkerEdges);
  assert(tdec->members[findMember(tdec, reducedMember->member)].type == TDEC_MEMBER_TYPE_BOND);

  if (depth == 0)
  {
    /* A bond root always works. */
    reducedMember->type = TYPE_5_ROOT;
    return TU_OKAY;
  }

  /* No children, but a path edge. */
  if (2*numTwoEnds + numOneEnd == 0 && reducedMember->firstPathEdge)
    reducedMember->type = TYPE_1_CLOSES_CYCLE;
  else if (numOneEnd == 1)
    reducedMember->type = reducedMember->firstPathEdge ? TYPE_2_SHORTCUT : TYPE_3_EXTENSION;
  else if (numOneEnd + 2 * numTwoEnds == 2)
  {
    if (reducedMember->firstPathEdge)
    {
      /* If there is a path edge then this should have been the root, but it is not! */
      newcolumn->remainsGraphic = false;
    }
    else
      reducedMember->type = TYPE_4_CONNECTS_TWO_PATHS;
  }
  else
  {
    /* Since no children contains path edges, the bond must be a leaf of the reduced decomposition and contains one. */
    assert(reducedMember->firstPathEdge);
    reducedMember->type = TYPE_1_CLOSES_CYCLE;
  }

  return TU_OKAY;
}

static
TU_ERROR determineTypePolygon(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int numOneEnd,                      /**< Number of child markers containing one path end. */
  int numTwoEnds,                     /**< Number of child markers containing two path ends. */
  TU_TDEC_EDGE childMarkerEdges[2],   /**< Marker edges of children containing one/two path ends. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedComponent);
  assert(reducedMember);
  assert(numOneEnd >= 0);
  assert(numTwoEnds >= 0);
  assert(numOneEnd + 2*numTwoEnds <= 2);
  assert(childMarkerEdges);
  assert(tdec->members[findMember(tdec, reducedMember->member)].type == TDEC_MEMBER_TYPE_POLYGON);

  TU_TDEC_MEMBER member = findMember(tdec, reducedMember->member);

  int countPathEdges = 0;
  for (PathEdge* edge = reducedMember->firstPathEdge; edge != NULL; edge = edge->next)
    ++countPathEdges;
  int numEdges = tdec->members[member].numEdges;

  TUdbgMsg(8+2*depth,
    "Determining type of polygon with %d edges, %d path edges, %d 1-end children and %d 2-end children.\n", numEdges,
    countPathEdges, numOneEnd, numTwoEnds);

  if (depth == 0)
  {
    /* We assume that we are not the root of the whole decomposition. */
    assert(tdec->members[member].parentMember >= 0);

    newcolumn->remainsGraphic = (numTwoEnds == 0);
    reducedMember->type = (countPathEdges == numEdges - 1) ? TYPE_1_CLOSES_CYCLE : TYPE_5_ROOT;
    return TU_OKAY;
  }

  if (countPathEdges == numEdges - 1)
  {
    reducedMember->type = TYPE_1_CLOSES_CYCLE;
  }
  else if (countPathEdges + numTwoEnds == numEdges - 1)
  {
    assert(numTwoEnds == 1);
    reducedMember->type = TYPE_4_CONNECTS_TWO_PATHS;
  }
  else if (numTwoEnds == 1)
  {
    newcolumn->remainsGraphic = false;
  }
  else if (numOneEnd == 1)
  {
    reducedMember->type = TYPE_3_EXTENSION;
  }
  else if (numOneEnd == 2)
  {
    reducedMember->type = TYPE_4_CONNECTS_TWO_PATHS;
  }
  else
  {
    assert(numOneEnd == 0);
    assert(numTwoEnds == 0);
    reducedMember->type = TYPE_3_EXTENSION;
  }

  return TU_OKAY;
}

static
TU_ERROR determineTypePrime(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int numOneEnd,                      /**< Number of child markers containing one path end. */
  int numTwoEnds,                     /**< Number of child markers containing two path ends. */
  TU_TDEC_EDGE childMarkerEdges[2],   /**< Marker edges of children containing one/two path ends. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedComponent);
  assert(reducedMember);
  assert(numOneEnd >= 0);
  assert(numTwoEnds >= 0);
  assert(numOneEnd + 2*numTwoEnds <= 2);
  assert(childMarkerEdges);
  assert(tdec->members[findMember(tdec, reducedMember->member)].type == TDEC_MEMBER_TYPE_PRIME);

  TU_TDEC_MEMBER member = findMember(tdec, reducedMember->member);

  /* Collect nodes of parent marker and of child markers containing one/two path end nodes. */
  TU_TDEC_NODE parentMarkerNodes[2] = {
    depth == 0 ? -1 : findEdgeTail(tdec, tdec->members[member].markerToParent),
    depth == 0 ? -1 : findEdgeHead(tdec, tdec->members[member].markerToParent)
  };
  TU_TDEC_NODE childMarkerNodes[4] = {
    childMarkerEdges[0] < 0 ? -1 : findEdgeTail(tdec, childMarkerEdges[0]),
    childMarkerEdges[0] < 0 ? -1 : findEdgeHead(tdec, childMarkerEdges[0]),
    childMarkerEdges[1] < 0 ? -1 : findEdgeTail(tdec, childMarkerEdges[1]),
    childMarkerEdges[1] < 0 ? -1 : findEdgeHead(tdec, childMarkerEdges[1])
  };

  TU_TDEC_NODE* endNodes = reducedMember->primeEndNodes;
  int numEndNodes = 0;

  /* Check the node degrees (with respect to path edges) in this component. */
  for (PathEdge* reducedEdge = reducedMember->firstPathEdge; reducedEdge; reducedEdge = reducedEdge->next)
  {
    TU_TDEC_NODE nodes[2] = { findEdgeHead(tdec, reducedEdge->edge), findEdgeTail(tdec, reducedEdge->edge) };
    for (int i = 0; i < 2; ++i)
    {
      TU_TDEC_NODE v = nodes[i];
      if (newcolumn->nodesDegree[v] >= 3)
      {
        TUdbgMsg(6 + 2*depth, "Node %d of prime member %d has path-degree at least 3.\n", v, reducedMember->member);
        newcolumn->remainsGraphic = false;
        return TU_OKAY;
      }

      if (newcolumn->nodesDegree[v] == 1)
      {
        if (numEndNodes == 4)
        {
          TUdbgMsg(6 + 2*depth, "Prime member %d has at least five path end nodes: %d, %d, %d, %d and %d.\n",
            reducedMember->member, endNodes[0], endNodes[1], endNodes[2], endNodes[3], v);
          newcolumn->remainsGraphic = false;
          return TU_OKAY;
        }

        endNodes[numEndNodes] = v;
        ++numEndNodes;
      }
    }
  }

  /* Exchange end nodes array using these rules:
   * If there are 4 edges, then one path connects 0 with 1 and the other 2 with 3.
   * For each path, parent marker nodes should come first (i.e., index 0 and 2).
   */

  if (numEndNodes == 4)
  {
    TU_TDEC_EDGE* nodeEdges = NULL;
    TUallocStackArray(tu, &nodeEdges, 2*tdec->memNodes);

    /* Initialize relevant entries to -1. */
    for (PathEdge* reducedEdge = reducedMember->firstPathEdge; reducedEdge; reducedEdge = reducedEdge->next)
    {
      TU_TDEC_NODE nodes[2] = { findEdgeHead(tdec, reducedEdge->edge), findEdgeTail(tdec, reducedEdge->edge) };
      for (int i = 0; i < 2; ++i)
      {
        TU_TDEC_NODE v = nodes[i];
        nodeEdges[2*v] = -1;
        nodeEdges[2*v + 1] = -1;
      }
    }

    /* Store incident edges for every node. */
    for (PathEdge* reducedEdge = reducedMember->firstPathEdge; reducedEdge; reducedEdge = reducedEdge->next)
    {
      TU_TDEC_NODE nodes[2] = { findEdgeHead(tdec, reducedEdge->edge), findEdgeTail(tdec, reducedEdge->edge) };
      for (int i = 0; i < 2; ++i)
      {
        TU_TDEC_NODE v = nodes[i];
        nodeEdges[2*v + (nodeEdges[2*v] == -1 ? 0 : 1)] = reducedEdge->edge;
      }
    }

    /* Start at end node 0 and see where we end. */
    TU_TDEC_EDGE previousEdge = -1;
    TU_TDEC_NODE currentNode = endNodes[0];
    while (true)
    {
      TU_TDEC_EDGE edge = nodeEdges[2*currentNode];
      if (edge == previousEdge)
        edge = nodeEdges[2*currentNode+1];
      if (edge == -1)
        break;
      previousEdge = edge;
      TU_TDEC_NODE v = findEdgeHead(tdec, edge);
      currentNode = (v != currentNode) ? v : findEdgeTail(tdec, edge);
    }
    TUfreeStackArray(tu, &nodeEdges);

    /* Exchange such that we end nodes 0 and 1 are end nodes of the same path. */
    if (currentNode == endNodes[2])
    {
      endNodes[2] = endNodes[1];
      endNodes[1] = currentNode;
    }
    else if (currentNode == endNodes[3])
    {
      endNodes[3] = endNodes[1];
      endNodes[1] = currentNode;
    }

    /* Exchange such that end node 2 is at the parent marker. */
    if (endNodes[2] != parentMarkerNodes[0] && endNodes[2] != parentMarkerNodes[1])
      SWAP_INTS(endNodes[2], endNodes[3]);
  }

  /* Exchange such that end node 0 is at the parent marker. */
  if (numEndNodes >= 2 && endNodes[0] != parentMarkerNodes[0] && endNodes[0] != parentMarkerNodes[1])
  {
    SWAP_INTS(endNodes[0], endNodes[1]);
  }

  TUdbgMsg(6 + 2*depth, "Prime member %d has %d path end nodes:", member, numEndNodes);
  if (numEndNodes >= 2)
    TUdbgMsg(0, " a path from %d to %d.", endNodes[0], endNodes[1]);
  if (numEndNodes == 4)
    TUdbgMsg(0, " a path from %d to %d.", endNodes[0], endNodes[1]);
  TUdbgMsg(0, "\n");

  if (depth == 0)
  {
    bool hasParentMarker = tdec->members[member].markerToParent >= 0;

    if (numEndNodes == 0)
    {
      /* No path edges, so there should be two adjacent child marker edges. */
      if (numOneEnd == 2 && (childMarkerNodes[0] == childMarkerNodes[2] || childMarkerNodes[0] == childMarkerNodes[3]
          || childMarkerNodes[1] == childMarkerNodes[2] || childMarkerNodes[2] == childMarkerNodes[3]))
      {
        reducedMember->type = TYPE_5_ROOT;
      }
      else
        newcolumn->remainsGraphic = false;
    }
    else if (numEndNodes == 2)
    {
      if (numOneEnd == 1)
      {
        if ((endNodes[0] == parentMarkerNodes[0] || endNodes[0] == parentMarkerNodes[1] || !hasParentMarker)
          && (endNodes[1] == childMarkerNodes[0] || endNodes[1] == childMarkerNodes[1]))
        {
          reducedMember->type = TYPE_5_ROOT;
        }
        else if ((endNodes[1] == parentMarkerNodes[0] || endNodes[1] == parentMarkerNodes[1] || !hasParentMarker)
          && (endNodes[0] == childMarkerNodes[0] || endNodes[0] == childMarkerNodes[1]))
        {
          reducedMember->type = TYPE_5_ROOT;
        }
        else
        {
          newcolumn->remainsGraphic = false;
        }
      }
      else if (numOneEnd == 2)
      {
        bool matched[2] = { false, false };
        for (int i = 0; i < 2; ++i)
        {
          for (int j = 0; j < 4; ++j)
          {
            if (reducedMember->primeEndNodes[i] == childMarkerNodes[j])
              matched[j/2] = true;
          }
        }
        
        if (matched[0] && matched[1])
          reducedMember->type = TYPE_5_ROOT;
        else
          newcolumn->remainsGraphic = false;
      }
      else if (numTwoEnds == 0)
      {
        assert(numOneEnd == 0);
        reducedMember->type = TYPE_5_ROOT;
      }
      else
      {
        assert(numOneEnd == 0);
        assert(numTwoEnds == 1);

        if ((childMarkerNodes[0] == endNodes[0] && childMarkerNodes[1] == endNodes[1])
          || (childMarkerNodes[0] == endNodes[1] && childMarkerNodes[1] == endNodes[0]))
        {
          reducedMember->type = TYPE_5_ROOT;
        }
        else
        {
          newcolumn->remainsGraphic = false;
        }
      }
    }
    else
    {
      assert(numEndNodes == 4);
      newcolumn->remainsGraphic = false;
    }
  }
  else
  {
    /* Non-root prime member. */

    int parentMarkerDegrees[2] = {
      newcolumn->nodesDegree[parentMarkerNodes[0]],
      newcolumn->nodesDegree[parentMarkerNodes[1]]
    };

    if (numEndNodes == 0)
    {
      /* We have no path edges, so there must be at least one child containing one/two path ends. */
      assert(numOneEnd + numTwoEnds > 0);
      /* We should not have a child marker edge parallel to the parent marker edge! */
      assert(!(parentMarkerNodes[0] == childMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[1])
        && !(parentMarkerNodes[0] == childMarkerNodes[1] && parentMarkerNodes[1] == childMarkerNodes[0]));

      if (numOneEnd == 0)
      {
        /* Even if parent and child marker (type 4) are adjacent, this is non-graphic. */
        newcolumn->remainsGraphic = false;
      }
      else if (numOneEnd == 1)
      {
        if (childMarkerNodes[0] == parentMarkerNodes[0] || childMarkerNodes[0] == parentMarkerNodes[1]
          || childMarkerNodes[1] == parentMarkerNodes[0] || childMarkerNodes[1] == parentMarkerNodes[1])
        {
          reducedMember->type = TYPE_3_EXTENSION;
        }
        else
        {
          newcolumn->remainsGraphic = false;
        }
      }
      else
      {
        /* We're not the root but have two proper children. */
        int childMarkerParentNode[2] = { -1, -1 };
        bool isParallel = false;
        for (int i = 0; i < 4; ++i)
        {
          for (int j = 0; j < 2; ++j)
          {
            if (childMarkerNodes[i] == parentMarkerNodes[j])
            {
              if (childMarkerParentNode[i/2] >= 0)
                isParallel = true;
              childMarkerParentNode[i/2] = j;
            }
          }
        }

        TUdbgMsg(6 + 2*depth, "No paths, %s, child[0] incident to %d and child[1] incident to %d.\n",
          isParallel ? "parallel" : "not parallel", childMarkerParentNode[0], childMarkerParentNode[1]);

        if (!isParallel && childMarkerParentNode[0] >= 0 && childMarkerParentNode[1] >= 0
          && childMarkerParentNode[0] != childMarkerParentNode[1])
        {
          reducedMember->type = TYPE_4_CONNECTS_TWO_PATHS;
        }
        else
          newcolumn->remainsGraphic = false;
      }
    }
    else if (numEndNodes == 2)
    {
      /* Exchange such that end node 0 is at the parent marker. */
      if (endNodes[0] != parentMarkerNodes[0] && endNodes[0] != parentMarkerNodes[1])
        SWAP_INTS(endNodes[0], endNodes[0]);

      if (numOneEnd == 1)
      {
        TUdbgMsg(6 + 2*depth, "%d-%d-path, parent {%d,%d} and child {%d,%d}\n", endNodes[0], endNodes[1],
          parentMarkerNodes[0], parentMarkerNodes[1], childMarkerNodes[0], childMarkerNodes[1]);

        if (parentMarkerNodes[0] != endNodes[0])
        {
          SWAP_INTS(parentMarkerNodes[0], parentMarkerNodes[1]);
          SWAP_INTS(parentMarkerDegrees[0], parentMarkerDegrees[1]);
        }
        if (parentMarkerNodes[0] != endNodes[0])
        {
          newcolumn->remainsGraphic = false;
          return TU_OKAY;
        }

        if (parentMarkerNodes[1] == endNodes[1])
        {
          /* Path closes a cycle with parent marker edge. */
          if (childMarkerNodes[0] == parentMarkerNodes[0] || childMarkerNodes[0] == parentMarkerNodes[1]
            || childMarkerNodes[1] == parentMarkerNodes[0] || childMarkerNodes[1] == parentMarkerNodes[1])
          {
            reducedMember->type = TYPE_2_SHORTCUT;
          }
          else
          {
            newcolumn->remainsGraphic = false;
          }
        }
        else
        {
          /* Path end is not incident to parent marker edge. */
          if (childMarkerNodes[0] == endNodes[1] || childMarkerNodes[1] == endNodes[1])
          {
            reducedMember->type = TYPE_3_EXTENSION;
          }
          else
          {
            newcolumn->remainsGraphic = false;
          }
        }
      }
      else if (numOneEnd == 2)
      {
        TU_TDEC_NODE otherParentNode = -1;
        if (endNodes[0] == parentMarkerNodes[0])
          otherParentNode = parentMarkerNodes[1];
        else if (endNodes[0] == parentMarkerNodes[1])
          otherParentNode = parentMarkerNodes[0];
        else
        {
          newcolumn->remainsGraphic = false;
          return TU_OKAY;
        }

        bool matched[4] = {
          childMarkerNodes[0] == endNodes[1] || childMarkerNodes[0] == otherParentNode,
          childMarkerNodes[1] == endNodes[1] || childMarkerNodes[1] == otherParentNode,
          endNodes[1] == childMarkerNodes[0] || endNodes[1] == childMarkerNodes[1]
            || endNodes[1] == childMarkerNodes[2] || endNodes[1] == childMarkerNodes[3],
          otherParentNode == childMarkerNodes[0] || otherParentNode == childMarkerNodes[1]
            || otherParentNode == childMarkerNodes[2] || otherParentNode == childMarkerNodes[3]
        };
        if (matched[0] && matched[1] && matched[2] && matched[3])
        {
          reducedMember->type = TYPE_4_CONNECTS_TWO_PATHS;
        }
        else
        {
          newcolumn->remainsGraphic = false;
        }
      }
      else if (numTwoEnds == 0)
      {
        if (parentMarkerDegrees[0] % 2 == 0 && parentMarkerDegrees[1] == 1)
        {
          reducedMember->type = parentMarkerDegrees[0] == 0 ? TYPE_3_EXTENSION : TYPE_2_SHORTCUT;
        }
        else if (parentMarkerDegrees[0] == 1 && parentMarkerDegrees[1] % 2 == 0)
        {
          reducedMember->type = parentMarkerDegrees[1] == 0 ? TYPE_3_EXTENSION : TYPE_2_SHORTCUT;
        }
        else if (parentMarkerDegrees[0] == 1 && parentMarkerDegrees[1] == 1)
        {
          reducedMember->type = TYPE_1_CLOSES_CYCLE;
        }
        else
        {
          /* Both have degree 0 or 2. */
          newcolumn->remainsGraphic = false;
        }
      }
      else
      {
        assert(numTwoEnds == 1);

        if ((endNodes[0] == parentMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[0]
          && childMarkerNodes[1] == endNodes[1])
          || (endNodes[0] == parentMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[1]
          && childMarkerNodes[0] == endNodes[1])
          || (endNodes[0] == parentMarkerNodes[1] && parentMarkerNodes[0] == childMarkerNodes[0]
          && childMarkerNodes[1] == endNodes[1])
          || (endNodes[0] == parentMarkerNodes[1] && parentMarkerNodes[0] == childMarkerNodes[1]
          && childMarkerNodes[0] == endNodes[1]))
        {
          reducedMember->type = TYPE_4_CONNECTS_TWO_PATHS;
        }
        else
        {
          newcolumn->remainsGraphic = false;
        }
      }
    }
    else if (numEndNodes == 4)
    {
      if (reducedMember->primeEndNodes[0] != parentMarkerNodes[0] && reducedMember->primeEndNodes[0] != parentMarkerNodes[1])
      {
        TUdbgMsg(6 + 2*depth, "First path does not start at parent marker edge.\n");
        newcolumn->remainsGraphic = false;
        return TU_OKAY;
      }
      if (reducedMember->primeEndNodes[2] != parentMarkerNodes[0] && reducedMember->primeEndNodes[2] != parentMarkerNodes[1])
      {
        TUdbgMsg(6 + 2*depth, "Second path does not start at parent marker edge.\n");
        newcolumn->remainsGraphic = false;
        return TU_OKAY;
      }

      if (numOneEnd == 0 && numTwoEnds == 0)
      {
        reducedMember->type = TYPE_4_CONNECTS_TWO_PATHS;
      }
      else if (numOneEnd == 1)
      {
        bool pathConnects[2] = {
          reducedMember->primeEndNodes[1] == childMarkerNodes[0]
            || reducedMember->primeEndNodes[1] == childMarkerNodes[1],
          reducedMember->primeEndNodes[3] == childMarkerNodes[0]
            || reducedMember->primeEndNodes[3] == childMarkerNodes[1]
        };

        if (pathConnects[0] && pathConnects[1])
        {
          TUdbgMsg(6 + 2*depth, "The both paths end at the child marker edge.\n");
          reducedMember->type = TYPE_4_CONNECTS_TWO_PATHS;
        }
        else if (pathConnects[0])
        {
          TUdbgMsg(6 + 2*depth, "Only the first path ends at the child marker edge.\n");
          reducedMember->type = TYPE_4_CONNECTS_TWO_PATHS;
        }
        else if (pathConnects[1])
        {
          TUdbgMsg(6 + 2*depth, "Only the second path ends at the child marker edge.\n");
          reducedMember->type = TYPE_4_CONNECTS_TWO_PATHS;
        }
        else
        {
          TUdbgMsg(6 + 2*depth, "No path ends at the child marker edge.\n");
          newcolumn->remainsGraphic = false;
        }
      }
      else
      {
        assert(numOneEnd == 2);

        bool pathConnected[2] = { false, false };
        bool childConnected[2] = { false, false };
        for (int i = 0; i < 2; ++i)
        {
          for (int j = 0; j < 4; ++j)
          {
            if (reducedMember->primeEndNodes[1 + 2*i] == childMarkerNodes[j])
            {
              pathConnected[i] = true;
              childConnected[j/2] = true;
            }
          }
        }
        if (pathConnected[0] && pathConnected[1] && childConnected[0] && childConnected[1])
        {
          reducedMember->type = TYPE_4_CONNECTS_TWO_PATHS;
        }
        else
        {
          TUdbgMsg(6 + 2*depth, "No pairing of paths to nodes of child marker edges possible.\n");
          newcolumn->remainsGraphic = false;
        }
      }
    }
  }

  return TU_OKAY;
}

static
TU_ERROR determineTypes(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);

  TUdbgMsg(4 + 2*depth, "determineTypes(member %d = reduced member %ld)\n", reducedMember->member,
    reducedMember - &newcolumn->reducedMembers[0]);

  /* First handle children recursively. */
  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    TU_CALL( determineTypes(tu, tdec, newcolumn, reducedComponent, reducedMember->children[c], depth + 1) );

    if (newcolumn->remainsGraphic)
    {
      TUdbgMsg(6 + 2*depth, "Child member %d of %d has type %d\n", reducedMember->children[c]->member, reducedMember->member,
        reducedMember->children[c]->type);
    }
    else
    {
      TUdbgMsg(6 + 2*depth, "Child prohibits graphicness.\n");
    }

    /* Abort if some part indicates non-graphicness. */
    if (!newcolumn->remainsGraphic)
      return TU_OKAY;
  }

  int numOneEnd, numTwoEnds;
  TU_TDEC_EDGE childMarkerEdges[2];
  TU_CALL( countChildrenTypes(tu, tdec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

#if defined(TU_DEBUG)
  int countReducedEdges = 0;
  for (PathEdge* e = reducedMember->firstPathEdge; e; e = e->next)
    ++countReducedEdges;
  TUdbgMsg(6 + 2*depth, "Member %d has %d children with one end and %d with two ends and %d path edges.\n",
    reducedMember->member, numOneEnd, numTwoEnds, countReducedEdges);
#endif /* TU_DEBUG */

  if (2*numTwoEnds + numOneEnd > 2)
  {
    newcolumn->remainsGraphic = false;
    return TU_OKAY;
  }

  bool isRoot = reducedMember == reducedComponent->root;
  TU_TDEC_MEMBER member = findMember(tdec, reducedMember->member);
  if (tdec->members[member].type == TDEC_MEMBER_TYPE_BOND)
  {
    TU_CALL( determineTypeBond(tu, tdec, newcolumn, reducedComponent, reducedMember, numOneEnd, numTwoEnds,
      childMarkerEdges, depth) );
  }
  else if (tdec->members[member].type == TDEC_MEMBER_TYPE_POLYGON)
  {
    TU_CALL( determineTypePolygon(tu, tdec, newcolumn, reducedComponent, reducedMember, numOneEnd, numTwoEnds,
      childMarkerEdges, depth) );
  }
  else
  {
    assert(tdec->members[member].type == TDEC_MEMBER_TYPE_PRIME);
    TU_CALL( determineTypePrime(tu, tdec, newcolumn, reducedComponent, reducedMember, numOneEnd, numTwoEnds,
      childMarkerEdges, depth) );
  }

  TUdbgMsg(6 + 2*depth, "Determined type %d.\n", reducedMember->type);

  /* Parent marker edge closes cycle, so we propagate information to parent. */

  if (!isRoot && reducedMember->type == TYPE_1_CLOSES_CYCLE)
  {
    TU_TDEC_MEMBER parentMember = findMemberParent(tdec, reducedMember->member);
    ReducedMember* reducedParent = newcolumn->membersToReducedMembers[parentMember];
    TU_TDEC_EDGE markerOfParent = tdec->members[member].markerOfParent;

    TUdbgMsg(6 + 2*depth, "Marker edge closes cycle.\n");
    TUdbgMsg(6 + 2*depth, "Parent member %d is reduced member %ld.\n", parentMember,
      (reducedParent - newcolumn->reducedMembers));

    /* Add marker edge of parent to reduced parent's path edges. */
    PathEdge* reducedEdge = NULL;
    TU_CALL( createPathEdge(tu, newcolumn, markerOfParent, reducedParent->firstPathEdge, &reducedEdge) );
    reducedParent->firstPathEdge = reducedEdge;

    /* Indicate that marker edge of parent belongs to path. */
    newcolumn->edgesInPath[markerOfParent] = true;

    /* Increase node degrees of nodes in a prime parent. */
    if (tdec->members[reducedParent->member].type == TDEC_MEMBER_TYPE_PRIME)
    {
      newcolumn->nodesDegree[findEdgeHead(tdec, markerOfParent)]++;
      newcolumn->nodesDegree[findEdgeTail(tdec, markerOfParent)]++;
    }

    TUdbgMsg(6 + 2*depth, "Added marker edge of parent to list of path edges.\n");
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

  TUdbgMsg(0, "\n  Checking whether we can add a column with %d 1's.\n", numEntries);

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

  TU_CALL( initializeNewColumn(tu, tdec, newcolumn) );
  TU_CALL( computeReducedDecomposition(tu, tdec, newcolumn, entryRows, numEntries) );
  TU_CALL( initializeReducedMemberEdgeLists(tu, tdec, newcolumn, entryRows, numEntries) );

  debugDot(tu, tdec, newcolumn);

  for (int i = 0; i < newcolumn->numReducedComponents; ++i)
  {
    TU_CALL( determineTypes(tu, tdec, newcolumn, &newcolumn->reducedComponents[i],
      newcolumn->reducedComponents[i].root, 0) );
  }

  if (newcolumn->remainsGraphic)
    TUdbgMsg(4, "Adding the column would maintain graphicness.\n");

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

  return TU_OKAY;
}

static
TU_ERROR addTerminal(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  TU_TDEC_MEMBER member,              /**< New terminal member. */
  TU_TDEC_NODE node                   /**< New terminal node. */
)
{
  assert(reducedComponent);
  assert(member >= 0);
  assert(isRepresentativeMember(tdec, member));

  /* For bonds we don't need to provide a node. */
  assert(node >= 0 || tdec->members[member].type == TDEC_MEMBER_TYPE_BOND);
  assert(reducedComponent->numTerminals != 1 || node >= 0 || member == reducedComponent->terminalMember[0]);

  TUdbgMsg(12, "Setting terminal node %d of member %d.\n", node, member);

  if (reducedComponent->numTerminals == 2)
  {
    TUdbgMsg(8, "Attempted to add terminal but already 2 known. We should detect non-graphicness soon.\n");
  }
  else
  {
    reducedComponent->terminalMember[reducedComponent->numTerminals] = member;
    reducedComponent->terminalNode[reducedComponent->numTerminals] = node;
    reducedComponent->numTerminals++;
  }

  return TU_OKAY;
}

static
TU_ERROR moveReducedRoot(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent  /**< Reduced component. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedComponent);

  ReducedMember* reducedMember = reducedComponent->root;
  TU_TDEC_MEMBER member = reducedMember->member;
  assert(isRepresentativeMember(tdec, member));

  int numOneEnd, numTwoEnds;
  TU_TDEC_EDGE childMarkerEdges[2];
  TU_CALL( countChildrenTypes(tu, tdec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  bool cycleWithUniqueEndChild;
  if (tdec->members[member].type == TDEC_MEMBER_TYPE_BOND)
  {
    cycleWithUniqueEndChild = (numTwoEnds == 1 || numOneEnd == 1);
  }
  else if (tdec->members[member].type == TDEC_MEMBER_TYPE_PRIME)
  {
    if (numTwoEnds == 1 || numOneEnd == 1)
    {
      /* For root primes, we have to check whether the path edges form a cycle with a two-end child. */
      TU_TDEC_NODE childMarkerNodes[2] = {
        findEdgeTail(tdec, childMarkerEdges[0]),
        findEdgeHead(tdec, childMarkerEdges[0])
      };
      cycleWithUniqueEndChild = (reducedMember->primeEndNodes[2] == -1
        && ((reducedMember->primeEndNodes[0] == childMarkerNodes[0]
        && reducedMember->primeEndNodes[1] == childMarkerNodes[1])
        || (reducedMember->primeEndNodes[0] == childMarkerNodes[1]
        && reducedMember->primeEndNodes[1] == childMarkerNodes[0])));
    }
    else
      cycleWithUniqueEndChild = false;
  }
  else
  {
    /* For root polygons, the parent marker is not a path edge, so we cannot close a cycle with a child marker edge. */
    cycleWithUniqueEndChild = false;
  }

  if (!cycleWithUniqueEndChild)
  {
    TUdbgMsg(6, "Reduced root member does not close a cycle with a 1- or 2-end child marker edge.\n");
    return TU_OKAY;
  }

  while (cycleWithUniqueEndChild)
  {
    TUdbgMsg(6, "Reduced member closes a cycle with a 1- or 2-end child marker edge.\n");

    TU_TDEC_MEMBER childMember = findMember(tdec, tdec->edges[childMarkerEdges[0]].childMember);
    newcolumn->edgesInPath[ tdec->members[childMember].markerToParent ] = true;

    /* Find the unique child member in order to process that. */
    for (int c = 0; c < reducedMember->numChildren; ++c)
    {
      Type type = reducedMember->children[c]->type;
      if (type == TYPE_2_SHORTCUT || type == TYPE_3_EXTENSION || type == TYPE_4_CONNECTS_TWO_PATHS)
      {
        reducedMember = reducedMember->children[c];
        break;
      }
    }

    member = reducedMember->member;
    assert(isRepresentativeMember(tdec, member));
    TU_CALL( countChildrenTypes(tu, tdec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

    if (tdec->members[member].type == TDEC_MEMBER_TYPE_BOND)
    {
      assert(!reducedMember->firstPathEdge);
      cycleWithUniqueEndChild = (numOneEnd == 1 || numTwoEnds == 1);
    }
    else if (tdec->members[member].type == TDEC_MEMBER_TYPE_PRIME)
    {
      cycleWithUniqueEndChild = ((numOneEnd == 1 || numTwoEnds == 1) && reducedMember->primeEndNodes[3] >= 0);
    }
    else
    {
      assert(tdec->members[member].type == TDEC_MEMBER_TYPE_POLYGON);
      if (numOneEnd == 1 || numTwoEnds == 1)
      {
        int countPathEdges = 1;
        for (PathEdge* pathEdge = reducedMember->firstPathEdge; pathEdge != NULL; pathEdge = pathEdge->next)
          ++countPathEdges;
        cycleWithUniqueEndChild = countPathEdges == tdec->members[member].numEdges - 1;
      }
      else
        cycleWithUniqueEndChild = false;
    }
  }

  TUdbgMsg(6, "Updated reduced root member is %d.\n", reducedMember->member);

  reducedComponent->root = reducedMember;

  return TU_OKAY;
}



static
TU_ERROR setEdgeNodes(
  TU* tu,             /**< \ref TU environment. */
  TU_TDEC* tdec,      /**< t-decomposition. */
  TU_TDEC_EDGE edge,  /**< Reduced component. */
  TU_TDEC_NODE tail,  /**< New tail node. */
  TU_TDEC_NODE head   /**< New head node. */
)
{
  assert(tu);
  assert(tdec);
  assert(edge >= 0);
  assert(edge < tdec->memEdges);

  tdec->edges[edge].tail = tail;
  tdec->edges[edge].head = head;

  return TU_OKAY;
}

static
TU_ERROR mergeMemberIntoParent(
  TU* tu,                 /**< \ref TU environment. */
  TU_TDEC* tdec,          /**< t-decomposition. */
  TU_TDEC_MEMBER member,  /**< Reduced component. */
  bool headToHead         /**< Whether the heads of the edges shall be joined. */
)
{
  assert(tu);
  assert(tdec);
  assert(member >= 0);

  member = findMember(tdec, member);
  TU_TDEC_MEMBER parentMember = findMemberParent(tdec, member);
  assert(parentMember >= 0);

  TUdbgMsg(10, "Merging child member %d into its parent member %d.\n", member, parentMember);

#if defined(TU_DEBUG)
  TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
  do
  {
    if (tdec->edges[edge].head < 0 || tdec->edges[edge].tail < 0)
      TUdbgMsg(10, "Edge %d of merge member %d does not have nodes.\n", edge, member);
    assert(tdec->edges[edge].tail >= 0);
    assert(tdec->edges[edge].head >= 0);
    edge = tdec->edges[edge].next;
  }
  while (edge != tdec->members[member].firstEdge);
#endif /* TU_DEBUG */

  TU_TDEC_EDGE parentEdge = tdec->members[member].markerOfParent;
  assert(parentEdge >= 0);
  TU_TDEC_EDGE childEdge = tdec->members[member].markerToParent;
  assert(childEdge >= 0);

  TU_TDEC_NODE parentEdgeNodes[2] = { findEdgeTail(tdec, parentEdge), findEdgeHead(tdec, parentEdge) };
  TU_TDEC_NODE childEdgeNodes[2] = { findEdgeTail(tdec, childEdge), findEdgeHead(tdec, childEdge) };
  TUdbgMsg(10, "Merging member %d into %d, identifying %d = {%d,%d} with %d = {%d,%d}.\n", member, parentMember,
    childEdge, childEdgeNodes[0], childEdgeNodes[1], parentEdge, parentEdgeNodes[0], parentEdgeNodes[1]);

  /* Identify nodes. */

  tdec->nodes[childEdgeNodes[0]].representativeNode = parentEdgeNodes[headToHead ? 0 : 1];
  tdec->nodes[childEdgeNodes[1]].representativeNode = parentEdgeNodes[headToHead ? 1 : 0];

  /* Identify members. */

  tdec->members[member].representativeMember = parentMember;

  /* We add the member's edges to the parent's edge list and thereby remove the two marker edges. */
  if (tdec->members[parentMember].firstEdge == parentEdge)
    tdec->members[parentMember].firstEdge = tdec->edges[parentEdge].next;

  tdec->edges[tdec->edges[parentEdge].next].prev = tdec->edges[childEdge].prev;
  tdec->edges[tdec->edges[parentEdge].prev].next = tdec->edges[childEdge].next;
  tdec->edges[tdec->edges[childEdge].next].prev = tdec->edges[parentEdge].prev;
  tdec->edges[tdec->edges[childEdge].prev].next = tdec->edges[parentEdge].next;
  tdec->members[parentMember].numEdges += tdec->members[member].numEdges - 2;
  tdec->numEdges -= 2;
  tdec->edges[parentEdge].next = tdec->firstFreeEdge;
  tdec->edges[childEdge].next = parentEdge;
  tdec->firstFreeEdge = childEdge;
  tdec->members[parentMember].type = TDEC_MEMBER_TYPE_PRIME;

  return TU_OKAY;
}

static
TU_ERROR createBondNodes(
  TU* tu,               /**< \ref TU environment. */
  TU_TDEC* tdec,        /**< t-decomposition. */
  TU_TDEC_MEMBER member /**< A bond member. */
)
{
  assert(tu);
  assert(tdec);
  assert(member >= 0);
  member = findMember(tdec, member);
  assert(tdec->members[member].type == TDEC_MEMBER_TYPE_BOND);

  TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
  assert(tdec->edges[edge].head < 0);
  assert(tdec->edges[edge].tail < 0);
  if (tdec->edges[edge].head >= 0)
  {
    assert(tdec->edges[edge].tail >= 0);
    return TU_OKAY;
  }

  TU_TDEC_NODE tail, head;
  TU_CALL( createNode(tu, tdec, &tail) );
  TU_CALL( createNode(tu, tdec, &head) );

  do
  {
    assert(tdec->edges[edge].tail < 0);
    assert(tdec->edges[edge].head < 0);

    tdec->edges[edge].tail = tail;
    tdec->edges[edge].head = head;
    edge = tdec->edges[edge].next;
  }
  while (edge != tdec->members[member].firstEdge);

  return TU_OKAY;
}

static
TU_ERROR splitBond(
  TU* tu,                     /**< \ref TU environment. */
  TU_TDEC* tdec,              /**< t-decomposition. */
  TU_TDEC_MEMBER bond,        /**< A bond member. */
  TU_TDEC_EDGE edge1,         /**< First edge to be moved into the child bond. */
  TU_TDEC_EDGE edge2,         /**< Second edge to be moved into the child bond. */
  TU_TDEC_MEMBER* pChildBond  /**< Pointer to storing the newly created child bond. */
)
{
  assert(tu);
  assert(tdec);
  assert(bond >= 0);
  assert(bond < tdec->memMembers);
  assert(edge1 >= 0);
  assert(edge1 < tdec->memEdges);
  assert(edge2 >= 0);
  assert(edge2 < tdec->memEdges);

  TU_TDEC_MEMBER childBond;
  TU_CALL( createMember(tu, tdec, TDEC_MEMBER_TYPE_BOND, &childBond) );
  TU_TDEC_EDGE markerOfParentBond, markerOfChildBond;
  TU_CALL( createMarkerEdge(tu, tdec, &markerOfParentBond, bond, -1, -1, true) );
  TU_CALL( addEdgeToMembersEdgeList(tu, tdec, markerOfParentBond, bond) );
  TU_CALL( createMarkerEdge(tu, tdec, &markerOfChildBond, childBond, -1, -1, false) );
  TU_CALL( addEdgeToMembersEdgeList(tu, tdec, markerOfChildBond, childBond) );
  tdec->numMarkers += 2;
  tdec->members[childBond].markerOfParent = markerOfParentBond;
  tdec->members[childBond].markerToParent = markerOfChildBond;
  tdec->members[childBond].parentMember = bond;
  tdec->edges[markerOfParentBond].childMember = childBond;
  tdec->edges[markerOfChildBond].childMember = -1;

  TU_CALL( removeEdgeFromMembersEdgeList(tu, tdec, edge1, bond) );
  tdec->edges[edge1].member = childBond;
  TU_CALL( addEdgeToMembersEdgeList(tu, tdec, edge1, childBond) );
  tdec->members[findMember(tdec, tdec->edges[edge1].childMember)].parentMember = childBond;

  TU_CALL( removeEdgeFromMembersEdgeList(tu, tdec, edge2, bond) );
  tdec->edges[edge2].member = childBond;
  TU_CALL( addEdgeToMembersEdgeList(tu, tdec, edge2, childBond) );
  tdec->members[findMember(tdec, tdec->edges[edge2].childMember)].parentMember = childBond;

  if (pChildBond)
    *pChildBond = childBond;

  return TU_OKAY;
}

static
TU_ERROR addColumnProcessBond(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of this member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedMember);
  TU_TDEC_MEMBER member = findMember(tdec, reducedMember->member);
  assert(tdec->members[member].type == TDEC_MEMBER_TYPE_BOND);

  int numOneEnd, numTwoEnds;
  TU_TDEC_EDGE childMarkerEdges[2];
  TU_CALL( countChildrenTypes(tu, tdec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  TUdbgMsg(6 + 2*depth, "addColumnProcessBond for%s member %d (reduced %ld), #one-ends = %d, #two-ends = %d of type %d.\n",
    depth == 0 ? " root" : "", member, (reducedMember - newcolumn->reducedMembers), numOneEnd, numTwoEnds, reducedMember->type);

  if (depth == 0)
  {
    if (numOneEnd == 0 && numTwoEnds == 0)
    {
      assert(reducedMember->firstPathEdge);
      TU_CALL( addTerminal(tu, tdec, reducedComponent, member, -1) );
      TU_CALL( addTerminal(tu, tdec, reducedComponent, member, -1) );

      return TU_OKAY;
    }
    else if (numOneEnd == 1)
    {
      assert(reducedComponent->numTerminals == 1);

      /* We don't merge since the child shall be the terminal member (rule (R3) in the paper)! */
      assert(reducedMember->firstPathEdge);
      TU_TDEC_MEMBER childMember = tdec->edges[childMarkerEdges[0]].childMember;
      TU_TDEC_EDGE childsParentMarker = tdec->members[childMember].markerToParent;
      TU_CALL( addTerminal(tu, tdec, reducedComponent, childMember, findEdgeTail(tdec, childsParentMarker)) );
      tdec->members[childMember].type = TDEC_MEMBER_TYPE_PRIME;
    }
    else if (numOneEnd == 2)
    {
      assert(reducedComponent->numTerminals == 2);

      /* If the bond contains more than 3 edges, we split the two child edges off, which creates a bond with a parent
       * marker and the two children.
       * Test: Graphic.RootBondTwoOneEnds */

      if (tdec->members[member].numEdges > 3)
      {
        TU_CALL( splitBond(tu, tdec, member, childMarkerEdges[0], childMarkerEdges[1], &member) );
        reducedMember->member = member;
      }

      debugDot(tu, tdec, newcolumn);

      assert(tdec->members[member].numEdges == 3);
      TU_CALL( createBondNodes(tu, tdec, member) );
      TU_CALL( mergeMemberIntoParent(tu, tdec, tdec->edges[childMarkerEdges[0]].childMember, true) );
      TU_CALL( mergeMemberIntoParent(tu, tdec, tdec->edges[childMarkerEdges[1]].childMember,
        reducedMember->firstPathEdge == NULL && reducedMember->type != TYPE_4_CONNECTS_TWO_PATHS) );

      debugDot(tu, tdec, newcolumn);

      return TU_OKAY;
    }
    else
    {
      assert(numTwoEnds == 1);
      assert(reducedMember->firstPathEdge);

      /* We don't merge since everything interesting happens in the child. */
      /* Test: Graphic.RootBondOneTwoEnd.  */
    }
  }
  else
  {
    /* Non-root bond. This cannot be a leaf since then it would contain a path edge, i.e., it closes a cycle. */

    assert(numOneEnd == 1);
    assert(reducedComponent->numTerminals >= 1);

    TU_TDEC_NODE tail, head;
    TU_CALL( createNode(tu, tdec, &tail) );
    TU_CALL( createNode(tu, tdec, &head) );
    TUdbgMsg(8 + 2*depth, "Bond's tail node is %d and head node is %d.\n", tail, head);
    TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
    do
    {
      TU_CALL( setEdgeNodes(tu, tdec, edge, tail, head) );
      edge = tdec->edges[edge].next;
    }
    while (edge != tdec->members[member].firstEdge);

    TU_CALL( mergeMemberIntoParent(tu, tdec, tdec->edges[childMarkerEdges[0]].childMember,
      !reducedMember->firstPathEdge) );

    debugDot(tu, tdec, newcolumn);
  }

  return TU_OKAY;
}

inline static
void flipEdge(
  TU_TDEC* tdec,    /**< t-decomposition. */
  TU_TDEC_EDGE edge /**< edge. */
)
{
  assert(tdec);
  assert(edge >= 0);
  assert(edge < tdec->memEdges);

  SWAP_INTS(tdec->edges[edge].tail, tdec->edges[edge].head);
}

static
TU_ERROR addColumnProcessPrime(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of this member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedMember);
  TU_TDEC_MEMBER member = findMember(tdec, reducedMember->member);
  assert(tdec->members[member].type == TDEC_MEMBER_TYPE_PRIME);

  int numOneEnd, numTwoEnds;
  TU_TDEC_EDGE childMarkerEdges[2];
  TU_CALL( countChildrenTypes(tu, tdec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  TUdbgMsg(6 + 2*depth, "addColumnProcessPrime for%s member %d (reduced %ld), #one-ends = %d, #two-ends = %d of type %d.\n",
    depth == 0 ? " root" : "", reducedMember->member, (reducedMember - newcolumn->reducedMembers), numOneEnd,
    numTwoEnds, reducedMember->type);

  TU_TDEC_NODE parentMarkerNodes[2] = {
    tdec->members[member].markerToParent < 0 ? -1 : findEdgeTail(tdec, tdec->members[member].markerToParent),
    tdec->members[member].markerToParent < 0 ? -1 : findEdgeHead(tdec, tdec->members[member].markerToParent)
  };
  TU_TDEC_NODE childMarkerNodes[4] = {
    childMarkerEdges[0] < 0 ? -1 : findEdgeTail(tdec, childMarkerEdges[0]),
    childMarkerEdges[0] < 0 ? -1 : findEdgeHead(tdec, childMarkerEdges[0]),
    childMarkerEdges[1] < 0 ? -1 : findEdgeTail(tdec, childMarkerEdges[1]),
    childMarkerEdges[1] < 0 ? -1 : findEdgeHead(tdec, childMarkerEdges[1])
  };

  TU_TDEC_NODE* endNodes = reducedMember->primeEndNodes;
  int numEndNodes = endNodes[0] < 0 ? 0 : (endNodes[2] < 0 ? 2 : 4 );

  if (depth == 0)
  {
    /* Root prime. */

    if (tdec->members[member].markerToParent >= 0 && newcolumn->edgesInPath[tdec->members[member].markerToParent])
    {
      /* The parent marker is a path edge, so we modify the endNodes array. */
      newcolumn->nodesDegree[parentMarkerNodes[0]]++;
      newcolumn->nodesDegree[parentMarkerNodes[1]]++;
      if (numEndNodes == 0)
      {
        endNodes[0] = parentMarkerNodes[0];
        endNodes[1] = parentMarkerNodes[1];
        numEndNodes = 1;
      }
      else if (numEndNodes == 2)
      {
        if (endNodes[0] == parentMarkerNodes[0])
          endNodes[0] = parentMarkerNodes[1];
        else if (endNodes[0] == parentMarkerNodes[1])
          endNodes[0] = parentMarkerNodes[0];
      }
      else
      {
        endNodes[0] = endNodes[3];
        endNodes[2] = -1;
        endNodes[3] = -1;
        numEndNodes = 2;
      }
    }

    assert(numEndNodes <= 2);

    if (numOneEnd == 0 && numTwoEnds == 0)
    {
      assert(numEndNodes >= 2);

      TU_CALL( addTerminal(tu, tdec, reducedComponent, member, endNodes[0]) );
      TU_CALL( addTerminal(tu, tdec, reducedComponent, member, endNodes[1]) );
    }
    else if (numOneEnd == 1)
    {
      TU_CALL( addTerminal(tu, tdec, reducedComponent, member,
        (endNodes[0] == childMarkerNodes[0] || endNodes[0] == childMarkerNodes[1]) ? endNodes[1] : endNodes[0] ) );

      TU_TDEC_MEMBER childMember = findMember(tdec, tdec->edges[childMarkerEdges[0]].childMember);
      bool headToHead = endNodes[0] == childMarkerNodes[1] || endNodes[1] == childMarkerNodes[1];
      assert(headToHead || endNodes[0] == childMarkerNodes[0] || endNodes[1] == childMarkerNodes[0]);

      TU_CALL( mergeMemberIntoParent(tu, tdec, childMember, headToHead) );
    }
    else
    {
      assert(numOneEnd == 2);
      assert(reducedComponent->numTerminals == 2);

      TU_TDEC_MEMBER childMember0 = findMember(tdec, tdec->edges[childMarkerEdges[0]].childMember);
      TU_TDEC_MEMBER childMember1 = findMember(tdec, tdec->edges[childMarkerEdges[1]].childMember);
      bool headToHead[2];
      for (int i = 0; i < 2; ++i)
      {
        for (int j = 0; j < 4; ++j)
        {
          if (endNodes[i] == childMarkerNodes[j])
            headToHead[j/2] = j % 2 == 1;
        }
      }

      TU_CALL( mergeMemberIntoParent(tu, tdec, childMember0, headToHead[0]) );
      TU_CALL( mergeMemberIntoParent(tu, tdec, childMember1, headToHead[1]) );
    }
  }
  else
  {
    /* Non-root prime. */

    if (numOneEnd == 0 && numTwoEnds == 0)
    {
      assert(reducedMember->firstPathEdge);
      assert(endNodes[0] >= 0);

      TU_CALL( addTerminal(tu, tdec, reducedComponent, member, endNodes[1]) );
      if (parentMarkerNodes[0] == endNodes[0])
        flipEdge(tdec, tdec->members[member].markerToParent);
    }
    else
    {
      assert(numOneEnd == 1);
      
      if (numEndNodes >= 2)
      {
        /* Flip parent if necessary. */
        if (endNodes[0] == parentMarkerNodes[0])
          flipEdge(tdec, tdec->members[member].markerToParent);

        /* Merge child. */
        TU_CALL( mergeMemberIntoParent(tu, tdec, findMember(tdec, tdec->edges[childMarkerEdges[0]].childMember),
          endNodes[1] == childMarkerNodes[1]) );
      }
      else
      {
        /* Parent marker and child marker must be next to each other. */
        if (parentMarkerNodes[0] == childMarkerNodes[0] || parentMarkerNodes[0] == childMarkerNodes[1])
          flipEdge(tdec, tdec->members[member].markerToParent);

        TU_CALL( mergeMemberIntoParent(tu, tdec, tdec->edges[childMarkerEdges[0]].childMember,
          parentMarkerNodes[0] == childMarkerNodes[1] || parentMarkerNodes[1] == childMarkerNodes[1]) );

        debugDot(tu, tdec, newcolumn);
      }
    }
  }

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
  TU_TDEC_EDGE edge,            /**< Edge to be replaced by a bond containing that edge. */
  TU_TDEC_MEMBER* pNewBond      /**< Pointer for storing the new bond. */
)
{
  assert(tu);
  assert(tdec);

  TUdbgMsg(8, "Creating bond for edge %d.\n", edge);

  TU_TDEC_MEMBER parentMember = findEdgeMember(tdec, edge);
  TU_TDEC_MEMBER newBond = -1;
  TU_CALL( createMember(tu, tdec, TDEC_MEMBER_TYPE_BOND, &newBond) );
  tdec->members[newBond].parentMember = parentMember;

  TU_TDEC_EDGE markerOfParent;
  TU_CALL( createMarkerEdge(tu, tdec, &markerOfParent, parentMember, tdec->edges[edge].tail,
    tdec->edges[edge].head, true) );
  tdec->edges[markerOfParent].childMember = newBond;
  tdec->edges[markerOfParent].next = tdec->edges[edge].next;
  tdec->edges[markerOfParent].prev = tdec->edges[edge].prev;
  assert(tdec->edges[markerOfParent].next != markerOfParent);
  tdec->edges[tdec->edges[markerOfParent].next].prev = markerOfParent;
  tdec->edges[tdec->edges[markerOfParent].prev].next = markerOfParent;
  if (tdec->members[parentMember].firstEdge == edge)
    tdec->members[parentMember].firstEdge = markerOfParent;
  tdec->members[newBond].markerOfParent = markerOfParent;

  TU_TDEC_EDGE markerToParent;
  TU_CALL( createMarkerEdge(tu, tdec, &markerToParent, newBond, -1, -1, false) );
  TU_CALL( addEdgeToMembersEdgeList(tu, tdec, markerToParent, newBond) );
  tdec->members[newBond].markerToParent = markerToParent;
  tdec->numMarkers++;

  tdec->edges[edge].member = newBond;
  TU_CALL( addEdgeToMembersEdgeList(tu, tdec, edge, newBond) );

  if (pNewBond)
    *pNewBond = newBond;

  return TU_OKAY;
}

/**
 * \brief Splits polygon into two, connecting them via a bond.
 *
 * Takes all edges of the polygon \p member for which \p edgesPredicate is the same as
 * \p predicateValue.
 */

static
TU_ERROR splitPolygon(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_MEMBER member,              /**< Polygon member to be squeezed. */
  bool* edgesPredicate,               /**< Map from edges to predicate. */
  bool predicateValue,                /**< Value of predicate. */
  TU_TDEC_EDGE* pRepresentativeEdge,  /**< Pointer for storing the child marker edge that links to new bond. */
  TU_TDEC_MEMBER* pNewBond,           /**< Pointer for storing the new bond. */
  TU_TDEC_MEMBER* pNewPolygon         /**< Pointer for storing the new polygon. */
)
{
  assert(tu);
  assert(tdec);
  assert(member >= 0);
  assert(member < tdec->memMembers);
  assert(edgesPredicate);
  assert(tdec->members[member].type == TDEC_MEMBER_TYPE_POLYGON);

#if defined(TU_DEBUG)
  TUdbgMsg(8, "Checking polygon %d for splitting... ", member);
#endif /* TU_DEBUG */

  TU_TDEC_EDGE edge = tdec->members[member].firstEdge;
  TU_TDEC_EDGE someSatisfyingEdge = -1;
  int numSatisfyingEdges = 0;
  do
  {
    bool value = edgesPredicate[edge];
    if ((value && predicateValue) || (!value && !predicateValue))
    {
      someSatisfyingEdge = edge;
      ++numSatisfyingEdges;
    }
    edge = tdec->edges[edge].next;
  } while (edge != tdec->members[member].firstEdge);

  TUdbgMsg(0, "%d of %d edges satisfy the predicate.\n", numSatisfyingEdges, tdec->members[member].numEdges);

  if (numSatisfyingEdges == 0)
  {
    if (pRepresentativeEdge)
      *pRepresentativeEdge = -1;
    if (pNewBond)
      *pNewBond = -1;
    if (pNewPolygon)
      *pNewPolygon = -1;
  }
  else if (numSatisfyingEdges == 1)
  {
    if (pNewPolygon)
      *pNewPolygon = -1;
    if (pRepresentativeEdge)
      *pRepresentativeEdge = someSatisfyingEdge;
    if (pNewBond)
      *pNewBond = -1;
  }
  else
  {
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
    edge = firstEdge;
    bool encounteredStayingEdge = false;
    do
    {
#if defined(TU_DEBUG_SPLITTING)
      TUdbgMsg(8, "Edge %d <%d>", edge, tdec->edges[edge].name);
      if (tdec->edges[edge].childMember >= 0)
        TUdbgMsg(0, " (with child %d)", tdec->edges[edge].childMember);
      if (edge == tdec->members[member].markerToParent)
        TUdbgMsg(0, " (with parent %d)", tdec->members[member].parentMember);
      TUdbgMsg(0, " (prev = %d, next = %d)", tdec->edges[edge].prev, tdec->edges[edge].next);
#endif /* TU_DEBUG_SPLITTING*/

      /* Evaluate predicate. */
      bool value = edgesPredicate[edge];
      if ((value && !predicateValue) || (!value && predicateValue))
      {
#if defined(TU_DEBUG_SPLITTING)
        TUdbgMsg(" does not satisfy the predicate.\n");
#endif /* TU_DEBUG_SPLITTING */
        edge = tdec->edges[edge].next;
        encounteredStayingEdge = true;
        continue;
      }

#if defined(TU_DEBUG_SPLITTING)
      TUdbgMsg(0, " satisfies the predicate.\n");
#endif /* TU_DEBUG_SPLITTING */

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

#if defined(TU_DEBUG_SPLITTING)
    TUdbgMsg(8, "Updated old polygon with these edges:");
    edge = firstEdge;
    do
    {
      TUdbgMsg(0, " %d <%d>", edge, tdec->edges[edge].name);
      edge = tdec->edges[edge].next;
    }
    while (edge != firstEdge);
    TUdbgMsg(0, ".\n");
    TUdbgMsg(8, "New polygon has these edges:");
    edge = polygonParentMarker;
    do
    {
      TUdbgMsg(0, " %d <%d>", edge, tdec->edges[edge].name);
      edge = tdec->edges[edge].next;
    }
    while (edge != polygonParentMarker);
    TUdbgMsg(0, ".\n");
#endif /* TU_DEBUG_SPLITTING */

    TUdbgMsg(8, "Connecting bond is member %d and squeezed polygon is member %d.\n", bond, polygon);

    if (pRepresentativeEdge)
      *pRepresentativeEdge = memberChildMarker;
    if (pNewBond)
      *pNewBond = bond;
    if (pNewPolygon)
      *pNewPolygon = polygon;
  }
  return TU_OKAY;
}

static
TU_ERROR addColumnProcessPolygon(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of this member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedMember);

  TU_TDEC_MEMBER member = findMember(tdec, reducedMember->member);

  int numOneEnd, numTwoEnds;
  TU_TDEC_EDGE childMarkerEdges[2];
  TU_CALL( countChildrenTypes(tu, tdec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  TUdbgMsg(6 + 2*depth, "addColumnProcessPolygon for%s member %d (reduced %ld), #one-ends = %d, #two-ends = %d of type %d.\n",
    depth == 0 ? " root" : "", member, (reducedMember - newcolumn->reducedMembers), numOneEnd, numTwoEnds, reducedMember->type);

  if (depth == 0)
  {
    /* If we have a child containing both ends then we should have moved the reduced root there. */
    assert(numTwoEnds == 0);

    if (numOneEnd == 0)
    {
      /* Root polygon containing both ends. */
      assert(reducedMember->firstPathEdge);
      TU_TDEC_EDGE representativeEdge;
      if (newcolumn->edgesInPath[tdec->members[member].markerToParent])
      {
        TUdbgMsg(8 + 2*depth, "Polygon contains both terminal nodes and the parent marker edge is a path edge.\n");

        /* Squeeze off all non-path edges by moving them to a new polygon and creating a bond to connect it to the
         * remaining polygon. */
        TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, false, &representativeEdge, NULL, NULL) );
      }
      else if (reducedMember->type == TYPE_1_CLOSES_CYCLE)
      {
        TUdbgMsg(8 + 2*depth, "Polygon contains both terminal nodes which are the parent marker edge nodes.\n");

        TU_TDEC_MEMBER parentMember = tdec->members[member].parentMember;
        TU_TDEC_EDGE markerOfParent = tdec->members[member].markerOfParent;
        assert(parentMember >= 0); /* A polygon can only close a cycle if it has a parent. */
        if (tdec->members[parentMember].type == TDEC_MEMBER_TYPE_BOND)
        {
          TU_CALL( addTerminal(tu, tdec, reducedComponent, findEdgeMember(tdec, markerOfParent), -1 ) );
          TU_CALL( addTerminal(tu, tdec, reducedComponent, findEdgeMember(tdec, markerOfParent), -1 ) );
        }
        else
        {
          assert(tdec->members[parentMember].type == TDEC_MEMBER_TYPE_PRIME);
          TU_CALL( addTerminal(tu, tdec, reducedComponent, findEdgeMember(tdec, markerOfParent), findEdgeTail(tdec, markerOfParent) ) );
          TU_CALL( addTerminal(tu, tdec, reducedComponent, findEdgeMember(tdec, markerOfParent), findEdgeHead(tdec, markerOfParent) ) );
        }

        return TU_OKAY;
      }
      else
      {
        TUdbgMsg(8 + 2*depth, "Polygon contains both terminal nodes and a non-path parent marker edge.\n");

        /* Squeeze off all path edges by moving them to a new polygon and creating a bond to connect it to the
         * remaining polygon. */
        TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, true, &representativeEdge, NULL, NULL) );
      }

      debugDot(tu, tdec, newcolumn);

      TU_TDEC_EDGE childMember = tdec->edges[representativeEdge].childMember;
      TU_TDEC_NODE tail = -1;
      TU_TDEC_NODE head = -1;
      if (childMember < 0)
      {
        TU_CALL( createEdgeBond(tu, tdec, representativeEdge, &childMember) );
        debugDot(tu, tdec, newcolumn);
      }
      else if (tdec->members[childMember].type == TDEC_MEMBER_TYPE_PRIME)
      {
        tail = findEdgeTail(tdec, tdec->members[childMember].markerToParent);
        head = findEdgeHead(tdec, tdec->members[childMember].markerToParent);
      }

      assert(reducedComponent->numTerminals == 0);
      TU_CALL( addTerminal(tu, tdec, reducedComponent, childMember, tail) );
      TU_CALL( addTerminal(tu, tdec, reducedComponent, childMember, head) );
      debugDot(tu, tdec, newcolumn);
    }
    else if (numOneEnd == 1)
    {
      /* Root polygon containing one end. */

      if (newcolumn->edgesInPath[tdec->members[member].markerToParent])
      {
        /* There is more than 1 path edge, so we split off all non-path edges and work in the new polygon. */
        if (reducedMember->firstPathEdge)
        {
          TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, false, NULL, NULL, &member) );
          reducedMember->member = member;
          newcolumn->edgesInPath[tdec->members[member].markerToParent] = true;
        }
        newcolumn->edgesInPath[childMarkerEdges[0]] = true;
        TU_TDEC_EDGE nonPathEdge;
        TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
        newcolumn->edgesInPath[childMarkerEdges[0]] = false;

        TU_TDEC_NODE a, b, c;
        TU_CALL( createNode(tu, tdec, &a) );
        TU_CALL( createNode(tu, tdec, &b) );
        if (tdec->members[member].numEdges == 3)
        {
          TU_CALL( createNode(tu, tdec, &c) );
          TU_CALL( setEdgeNodes(tu, tdec, nonPathEdge, a, c) );
        }
        else
          c = a;
        TU_CALL( setEdgeNodes(tu, tdec, tdec->members[member].markerToParent, a, b) );
        TU_CALL( setEdgeNodes(tu, tdec, childMarkerEdges[0], c, b) );
        TU_CALL( addTerminal(tu, tdec, reducedComponent, member, a) );
        TU_CALL( mergeMemberIntoParent(tu, tdec, tdec->edges[childMarkerEdges[0]].childMember, true) );
        tdec->members[member].type = TDEC_MEMBER_TYPE_PRIME;
      }
      else
      {
        /* Parent marker edge is not a path edge. */

        assert(reducedMember->firstPathEdge);

        /* If there is more than 1 path edge, we squeeze off by moving them to a new polygon and creating a bond to
         * connect it to the remaining polygon. */
        TU_TDEC_EDGE pathEdge;
        TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );
        if (pathEdge >= 0)
          newcolumn->edgesInPath[pathEdge] = true;
        debugDot(tu, tdec, newcolumn);

        /* Unless the polygon consists of only the parent marker, the child marker (containing a path end) and a
         * representative edge, we squeeze off the representative edge and the child marker. */
        if (tdec->members[member].numEdges > 3)
        {
          newcolumn->edgesInPath[childMarkerEdges[0]] = true;
          TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, true, NULL, NULL, &member) );
          newcolumn->edgesInPath[childMarkerEdges[0]] = false;
          debugDot(tu, tdec, newcolumn);
        }

        assert(tdec->members[member].numEdges == 3);

        TU_TDEC_NODE a, b, c;
        TU_CALL( createNode(tu, tdec, &a) );
        TU_CALL( createNode(tu, tdec, &b) );
        TU_CALL( createNode(tu, tdec, &c) );
        TU_CALL( setEdgeNodes(tu, tdec, tdec->members[member].markerToParent, b, c) );
        TU_CALL( setEdgeNodes(tu, tdec, pathEdge, a, b) );
        TU_CALL( setEdgeNodes(tu, tdec, childMarkerEdges[0], c, a) );
        TU_CALL( addTerminal(tu, tdec, reducedComponent, member, b) );
        TU_CALL( mergeMemberIntoParent(tu, tdec, tdec->edges[childMarkerEdges[0]].childMember, true) );

        assert("TODO: TESTME");
      }
      debugDot(tu, tdec, newcolumn);
    }
    else
    {
      assert(numOneEnd == 2);

      /* If there is more than 1 path edge, we split off by moving them to a new polygon and creating a bond to
       * connect it to the remaining polygon. */
      TU_TDEC_EDGE pathEdge = -1;
      TU_TDEC_EDGE nonPathEdge = -1;
      if (reducedMember->type != TYPE_4_CONNECTS_TWO_PATHS)
      {
        /* Parent marker is a non-path edge. We split off path edges if more than one. */
        TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );

        /* If there is another non-path edge, we split off the polygon consisting of
         * the two relevant child marker edges and possibly the path edge. We then replace the current polygon by the
         * new one. */

        if (tdec->members[member].numEdges > (pathEdge >= 0 ? 4 : 3))
        {
          if (pathEdge >= 0)
            newcolumn->edgesInPath[pathEdge] = true;
          newcolumn->edgesInPath[childMarkerEdges[0]] = true;
          newcolumn->edgesInPath[childMarkerEdges[1]] = true;
          TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, true, NULL, NULL, &member) );
          reducedMember->member = member;
        }

        nonPathEdge = tdec->members[member].markerToParent;
      }
      else 
      {
        assert(newcolumn->edgesInPath[tdec->members[member].markerToParent]);

        if (reducedMember->firstPathEdge)
        {
          /* Parent marker is one of several path edges. */
          TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, false, NULL, NULL, &member) );
          reducedMember->member = member;
        }
        pathEdge = tdec->members[member].markerToParent;

        if (tdec->members[member].numEdges > 3)
        {
          newcolumn->edgesInPath[childMarkerEdges[0]] = true;
          newcolumn->edgesInPath[childMarkerEdges[1]] = true;
          TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
          newcolumn->edgesInPath[childMarkerEdges[0]] = false;
          newcolumn->edgesInPath[childMarkerEdges[1]] = false;
        }
        else
          nonPathEdge = -1;
      }

      TUdbgMsg(8 + 2*depth,
        "After splitting off, the (potential) path edge is %d and the (potential) non-path edge is %d.\n",
        pathEdge, nonPathEdge);
      debugDot(tu, tdec, newcolumn);
      
      assert(pathEdge >= 0 || nonPathEdge >= 0);

      /* a <----- b ---- c -----> d -------- a
       *   child0   path   child1   non-path
       *
       * or
       *
       * a <----- b ---- c -----> d=a
       *   child0   path   child1
       *
       * or
       *
       * a <----- b=c -----> d -------- a
       *   child0     child1   non-path */
      TU_TDEC_NODE a, b, c, d;
      TU_CALL( createNode(tu, tdec, &a) );
      TU_CALL( createNode(tu, tdec, &b) );
      if (pathEdge >= 0)
        TU_CALL( createNode(tu, tdec, &c) );
      else
        c = b;
      if (nonPathEdge >= 0)
        TU_CALL( createNode(tu, tdec, &d) );
      else
        d = a;
      TU_CALL( setEdgeNodes(tu, tdec, childMarkerEdges[0], a, b) );
      TUdbgMsg(8 + 2*depth, "First child edge %d = {%d, %d}\n", childMarkerEdges[0], a, b);
      TU_CALL( setEdgeNodes(tu, tdec, childMarkerEdges[1], d, c) );
      TUdbgMsg(8 + 2*depth, "Second child edge %d = {%d, %d}\n", childMarkerEdges[1], d, c);
      if (pathEdge >= 0)
      {
        TU_CALL( setEdgeNodes(tu, tdec, pathEdge, b, c) );
        TUdbgMsg(8 + 2*depth, "Path edge %d = {%d, %d}\n", pathEdge, b, c);
      }
      if (nonPathEdge >= 0)
      {
        TU_CALL( setEdgeNodes(tu, tdec, nonPathEdge, d, a) );
        TUdbgMsg(8 + 2*depth, "Non-path edge %d = {%d, %d}\n", nonPathEdge, d, a);
      }

      TU_CALL( mergeMemberIntoParent(tu, tdec, tdec->edges[childMarkerEdges[0]].childMember, true) );
      TU_CALL( mergeMemberIntoParent(tu, tdec, tdec->edges[childMarkerEdges[1]].childMember, true) );
      tdec->members[member].type = TDEC_MEMBER_TYPE_PRIME;
    }
  }
  else
  {
    /* Non-root polygon. */
    assert(reducedMember->type != TYPE_1_CLOSES_CYCLE); /* addColumnProcess should never consider such a member. */
    assert(reducedMember->type != TYPE_2_SHORTCUT); /* For polygons this can never happen. */
    assert(reducedMember->type != TYPE_4_CONNECTS_TWO_PATHS); /* This should only happen at the root. */
    assert(reducedMember->type != TYPE_5_ROOT); /* We are not a root. */
    assert(reducedMember->type == TYPE_3_EXTENSION); /* Only remaining case. */

    if (numOneEnd == 0)
    {
      assert(numOneEnd == 0);
      assert(numTwoEnds == 0);
      assert(reducedComponent->numTerminals < 2);
      assert(reducedMember->firstPathEdge);

      TUdbgMsg(6 + 2*depth, "Non-root polygon %d without one-end children.\n", member);

      /* Squeeze off all path edges by moving them to a new polygon and creating a bond to connect
       * it to the remaining polygon. */
      TU_TDEC_EDGE pathEdge = -1;
      TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );
      assert(pathEdge >= 0);
      newcolumn->edgesInPath[pathEdge] = true;

      /* If necessary, we squeeze off the non-path edges as well. */
      assert(tdec->members[member].numEdges >= 3);
      TU_TDEC_EDGE nonPathEdge;
      if (tdec->members[member].numEdges == 3)
      {
        nonPathEdge = tdec->edges[tdec->members[member].markerToParent].next;
        if (nonPathEdge == pathEdge)
          nonPathEdge = tdec->edges[nonPathEdge].next;
      }
      else
      {
        /* We temporarily mark the parent edge to belong to the path. */
        TU_TDEC_EDGE markerToParent = tdec->members[member].markerToParent;
        newcolumn->edgesInPath[markerToParent] = true;
        TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
        newcolumn->edgesInPath[markerToParent] = false;
        debugDot(tu, tdec, newcolumn);
      }
      assert(tdec->members[member].numEdges == 3);

      /* We now create the nodes of the triangle so that the path leaves it via the parent marker edge's head node. */

      TU_TDEC_NODE a, b, c;
      TU_CALL( createNode(tu, tdec, &a) );
      TU_CALL( createNode(tu, tdec, &b) );
      TU_CALL( createNode(tu, tdec, &c) );
      TU_CALL( setEdgeNodes(tu, tdec, tdec->members[reducedMember->member].markerToParent, a, b) );
      TU_CALL( setEdgeNodes(tu, tdec, pathEdge, b, c) );
      TU_CALL( setEdgeNodes(tu, tdec, nonPathEdge, c, a) );
      TU_CALL( addTerminal(tu, tdec, reducedComponent, reducedMember->member, c) );

      return TU_OKAY;
    }
    else
    {
      assert(numOneEnd == 1);

      /* Squeeze off all path edges by moving them to a new polygon and creating a bond to connect
       * it to the remaining polygon. */
      TU_TDEC_EDGE pathEdge = -1;
      TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );

      /* If necessary, we squeeze off the non-path edges as well. */
      assert(tdec->members[member].numEdges >= 3);
      int numNonPathEdges = tdec->members[member].numEdges - 2 - (pathEdge >= 0 ? 1 : 0);
      TUdbgMsg(8 + 2*depth, "Number of non-path edges is %d.\n", numNonPathEdges);
      TU_TDEC_EDGE nonPathEdge;
      if (numNonPathEdges == 0)
        nonPathEdge = -1;
      else if (numNonPathEdges == 1)
      {
        nonPathEdge = tdec->edges[tdec->members[member].markerToParent].next;
        while (nonPathEdge == childMarkerEdges[0] || nonPathEdge == pathEdge)
          nonPathEdge = tdec->edges[nonPathEdge].next;
      }
      else
      {
        newcolumn->edgesInPath[tdec->members[member].markerToParent] = true;
        newcolumn->edgesInPath[childMarkerEdges[0]] = true;
        TU_CALL( splitPolygon(tu, tdec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
        newcolumn->edgesInPath[tdec->members[member].markerToParent] = false;
        newcolumn->edgesInPath[childMarkerEdges[0]] = false;
        debugDot(tu, tdec, newcolumn);
      }

      TUdbgMsg(8 + 2*depth, "After (potential) splitting: path edge is %d and non-path edge is %d.\n",
        pathEdge, nonPathEdge);
      assert(tdec->members[member].numEdges <= 4);

      /* We now create the nodes of the triangle so that the path leaves it via the parent marker edge's head node. */

      TU_TDEC_NODE a, b, c, d;
      TU_CALL( createNode(tu, tdec, &a) );
      TU_CALL( createNode(tu, tdec, &b) );
      TU_CALL( createNode(tu, tdec, &c) );
      TU_CALL( setEdgeNodes(tu, tdec, tdec->members[reducedMember->member].markerToParent, a, b) );
      if (tdec->members[member].numEdges == 4)
      {
        TU_CALL( createNode(tu, tdec, &d) );
        TU_CALL( setEdgeNodes(tu, tdec, pathEdge, b, c) );
        TU_CALL( setEdgeNodes(tu, tdec, childMarkerEdges[0], d, c) );
        TU_CALL( setEdgeNodes(tu, tdec, nonPathEdge, d, a) );
      }
      else if (nonPathEdge == -1)
      {
        TU_CALL( setEdgeNodes(tu, tdec, pathEdge, b, c) );
        TU_CALL( setEdgeNodes(tu, tdec, childMarkerEdges[0], a, c) );
      }
      else
      {
        assert(pathEdge == -1);
        TU_CALL( setEdgeNodes(tu, tdec, childMarkerEdges[0], c, b) );
        TU_CALL( setEdgeNodes(tu, tdec, nonPathEdge, a, c) );
      }
      debugDot(tu, tdec, newcolumn);

      TU_CALL( mergeMemberIntoParent(tu, tdec, tdec->edges[childMarkerEdges[0]].childMember, true) );
      debugDot(tu, tdec, newcolumn);

      return TU_OKAY;
    }
  }

  return TU_OKAY;
}

/**
 * \brief Process reduced t-decomposition before the actual modification.
 *
 * Processes the reduced members in depth-first search manner and does the following:
 * - Polygons are squeezed.
 * - Terminal nodes and (reduced) members are detected.
 * - Marker edges along unique path between terminal nodes are merged.
 */

static
TU_ERROR addColumnProcessComponent(
  TU* tu,                             /**< \ref TU environment. */
  TU_TDEC* tdec,                      /**< t-decomposition. */
  TU_TDEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(tdec);
  assert(newcolumn);
  assert(reducedComponent);

  TUdbgMsg(6 + 2*depth, "addColumnProcessComponent(member %d = reduced member %ld)\n", reducedMember->member,
    (reducedMember - &newcolumn->reducedMembers[0]));

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

  /* If we are non-root type 1, then we don't need to do anything. */
  if (reducedMember->type == TYPE_1_CLOSES_CYCLE && depth > 0)
  {
    return TU_OKAY;
  }

  /* Handle children recursively. */
  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    ReducedMember* child = reducedMember->children[c];
    if (child->type != TYPE_1_CLOSES_CYCLE)
    {
      TU_CALL( addColumnProcessComponent(tu, tdec, newcolumn, reducedComponent, reducedMember->children[c], depth+1) );
    }
    else
    {
      TUdbgMsg(8 + 2*depth, "Member %d is implicitly replaced by a path edge.\n", findMember(tdec, child->member));
    }
  }

  /* Different behavior for bonds, polygons and prime components. */
  if (tdec->members[reducedMember->member].type == TDEC_MEMBER_TYPE_BOND)
    TU_CALL( addColumnProcessBond(tu, tdec, newcolumn, reducedComponent, reducedMember, depth) );
  else if (tdec->members[reducedMember->member].type == TDEC_MEMBER_TYPE_POLYGON)
    TU_CALL( addColumnProcessPolygon(tu, tdec, newcolumn, reducedComponent, reducedMember, depth) );
  else
    TU_CALL( addColumnProcessPrime(tu, tdec, newcolumn, reducedComponent, reducedMember, depth) );

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

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

  TUdbgMsg(8, "Flipping parenting of member %d with old parent %d and new parent %d.\n", member, oldParent, newParent);

  tdec->members[member].markerToParent = newMarkerToParent;
  tdec->members[member].markerOfParent = markerOfNewParent;
  tdec->members[member].parentMember = newParent;
  tdec->edges[markerOfNewParent].childMember = member;
  tdec->edges[newMarkerToParent].childMember = -1;

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

  TUdbgMsg(4, "Making member %d the new root of its component.\n", newRoot);

  if (tdec->members[newRoot].parentMember >= 0)
  {
    TU_CALL( doReorderComponent(tu, tdec, findMemberParent(tdec, newRoot), newRoot,
      tdec->members[newRoot].markerOfParent, tdec->members[newRoot].markerToParent) );
  }

  return TU_OKAY;
}

static
TU_ERROR doMergeLeafBonds(
  TU* tu,       /**< \ref TU environment. */
  TU_TDEC* tdec /**< t-decomposition. */
)
{
  assert(tu);
  assert(tdec);

  for (TU_TDEC_MEMBER member = 0; member < tdec->numMembers; ++member)
  {
    if (isRepresentativeMember(tdec, member) && tdec->members[member].type == TDEC_MEMBER_TYPE_BOND
      && tdec->members[member].numEdges == 2 && tdec->members[member].parentMember >= 0)
    {
      TU_TDEC_MEMBER parentMember = tdec->members[member].parentMember;
      TU_TDEC_EDGE parentEdge = tdec->members[member].markerOfParent;
      TU_TDEC_EDGE childEdge = tdec->members[member].markerToParent;
      TU_TDEC_EDGE otherBondEdge = tdec->edges[childEdge].next;
      TUdbgMsg(4, "Merging bond %d into its parent %d.\n", member, parentMember);

      /* We just use the nodes of the parent's child marker (even if -1). */
      tdec->edges[otherBondEdge].head = tdec->edges[parentEdge].head;
      tdec->edges[otherBondEdge].tail = tdec->edges[parentEdge].tail;

      /* Identify members. */
      tdec->members[member].representativeMember = parentMember;

      /* We replace the parent's child marker edge by the members non-parent marker edge and remove the two marker edges. */
      if (tdec->members[parentMember].firstEdge == parentEdge)
        tdec->members[parentMember].firstEdge = tdec->edges[parentEdge].next;
      tdec->edges[tdec->edges[parentEdge].next].prev = otherBondEdge;
      tdec->edges[tdec->edges[parentEdge].prev].next = otherBondEdge;
      tdec->edges[otherBondEdge].prev = tdec->edges[parentEdge].prev;
      tdec->edges[otherBondEdge].next = tdec->edges[parentEdge].next;
      tdec->edges[otherBondEdge].member = findMember(tdec, parentMember);
      tdec->numEdges -= 2;
      tdec->edges[parentEdge].next = tdec->firstFreeEdge;
      tdec->edges[childEdge].next = parentEdge;
      tdec->firstFreeEdge = childEdge;
    }
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

  TUdbgMsg(0, "\n  Adding a column with %d 1's.\n", numEntries);

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

  /* Create reduced components for new edges. */
  TU_CALL( completeReducedDecomposition(tu, tdec, newcolumn, entryRows, numEntries) );

  /* Process all reduced components individually. */

  TU_TDEC_EDGE* componentNewEdges = NULL;
  TU_CALL( TUallocStackArray(tu, &componentNewEdges, newcolumn->numReducedComponents) );

  int maxDepthComponent = -1;
  TUdbgMsg(4, "Processing %d reduced components.\n", newcolumn->numReducedComponents);
  for (int i = 0; i < newcolumn->numReducedComponents; ++i)
  {
    ReducedComponent* reducedComponent = &newcolumn->reducedComponents[i];

    TUdbgMsg(4, "Moving root of reduced component %d.\n", i);

    TU_CALL( moveReducedRoot(tu, tdec, newcolumn, reducedComponent) );
    debugDot(tu, tdec, newcolumn);

    TUdbgMsg(4, "Processing reduced component %d of depth %d.\n", i, reducedComponent->rootDepth);

    if (maxDepthComponent < 0 || reducedComponent->rootDepth
      > newcolumn->reducedComponents[maxDepthComponent].rootDepth)
    {
      maxDepthComponent = i;
    }

    TU_CALL( addColumnProcessComponent(tu, tdec, newcolumn, reducedComponent, reducedComponent->root, 0) );

    assert(reducedComponent->numTerminals == 2);
    TUdbgMsg(6, "Terminal members are %d and %d.\n", reducedComponent->terminalMember[0],
      reducedComponent->terminalMember[1]);
    TUdbgMsg(6, "Terminal nodes are %d and %d.\n", reducedComponent->terminalNode[0],
      reducedComponent->terminalNode[1]);
    assert(findMember(tdec, reducedComponent->terminalMember[0]) == findMember(tdec, reducedComponent->terminalMember[1]));

    /* Create new edge for this component. If there is one component, this is a column edge, and
     * otherwise it is a marker edge that will be linked to a new polygon consisting of all these
     * marker edges and the column edge. */
    TU_TDEC_EDGE newEdge;
    TU_CALL( createEdge(tu, tdec, -1, &newEdge) );
    componentNewEdges[i] = newEdge;
    tdec->edges[newEdge].childMember = -1;
    tdec->edges[newEdge].member = findMember(tdec, reducedComponent->terminalMember[0]);
    tdec->edges[newEdge].head = reducedComponent->terminalNode[0];
    tdec->edges[newEdge].tail = reducedComponent->terminalNode[1];
    tdec->edges[newEdge].name = INT_MAX/2;
    TU_CALL( addEdgeToMembersEdgeList(tu, tdec, newEdge, tdec->edges[newEdge].member) );
  }

  for (int c = 0; c < newcolumn->numReducedComponents; ++c)
  {
    if (c == maxDepthComponent)
      TUdbgMsg(6, "Reduced component %d has maximum depth and will remain a root.\n", c);
    else
      TU_CALL( reorderComponent(tu, tdec, findMember(tdec, tdec->edges[componentNewEdges[c]].member)) );
  }

  if (newcolumn->numReducedComponents == 0)
  {
    TU_TDEC_MEMBER loopMember;
    TU_CALL( createMember(tu, tdec, TDEC_MEMBER_TYPE_LOOP, &loopMember) );
    TU_TDEC_MEMBER loopEdge;
    TU_CALL( createEdge(tu, tdec, loopMember, &loopEdge) );
    TU_CALL( addEdgeToMembersEdgeList(tu, tdec, loopEdge, loopMember) );
    tdec->edges[loopEdge].name = -1 - column;
    tdec->edges[loopEdge].childMember = -1;
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
    tdec->edges[columnEdge].name = -1 - column;
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

  debugDot(tu, tdec, newcolumn);

  TU_CALL( TUfreeStackArray(tu, &componentNewEdges) );

  newcolumn->numReducedMembers = 0;
  newcolumn->numReducedComponents = 0;

  TUconsistencyAssert( TUtdecConsistency(tu, tdec) );

  return TU_OKAY;
}

TU_ERROR testGraphicnessTDecomposition(TU* tu, TU_CHRMAT* matrix, TU_CHRMAT* transpose,
  bool* pisGraphic, TU_GRAPH* graph, TU_GRAPH_EDGE* basis, TU_GRAPH_EDGE* cobasis,
  TU_SUBMAT** psubmatrix, int mergeLeafBonds)
{
  assert(tu);
  assert(matrix);
  assert(transpose);
  assert(!graph || (TUgraphNumNodes(graph) == 0 && TUgraphNumEdges(graph) == 0));
  assert(!psubmatrix || !*psubmatrix);
  assert(!basis || graph);
  assert(!cobasis || graph);
  assert(mergeLeafBonds >= 0);
  assert(mergeLeafBonds <= 2);

#if defined(TU_DEBUG)
  TUdbgMsg(0, "testGraphicnessTDecomposition called for a 1-connected %dx%d matrix.\n", matrix->numRows,
    matrix->numColumns);
  TUchrmatPrintDense(stdout, (TU_CHRMAT*) matrix, '0', true);
#endif /* TU_DEBUG */

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

      TUdbgMsg(0, "Constructed graph with %d nodes and %d edges.\n", TUgraphNumNodes(graph), TUgraphNumEdges(graph));
    }
    return TU_OKAY;
  }

  TU_TDEC* tdec = NULL;
  TU_CALL( TUtdecCreate(tu, &tdec, 0, 0, 0, 0, 0) ); /* TODO: avoid reallocations. */

  /* Process each column. */
  TU_TDEC_NEWCOLUMN* newcolumn = NULL;
  TUtdecnewcolumnCreate(tu, &newcolumn);
  *pisGraphic = true;
  for (int column = 0; column < matrix->numColumns && *pisGraphic; ++column)
  {
    TU_CALL( TUtdecAddColumnCheck(tu, tdec, newcolumn,
      &transpose->entryColumns[transpose->rowStarts[column]],
      transpose->rowStarts[column+1] - transpose->rowStarts[column]) );

    debugDot(tu, tdec, newcolumn);

    if (newcolumn->remainsGraphic)
    {
      TU_CALL( TUtdecAddColumnApply(tu, tdec, newcolumn, column, &transpose->entryColumns[transpose->rowStarts[column]],
        transpose->rowStarts[column+1] - transpose->rowStarts[column]) );

      if (mergeLeafBonds == 2)
        TU_CALL( doMergeLeafBonds(tu, tdec) );

      debugDot(tu, tdec, newcolumn);
    }
    else
    {
      *pisGraphic = false;
    }
  }

  if (*pisGraphic && mergeLeafBonds > 0)
  {
    TU_CALL( doMergeLeafBonds(tu, tdec) );
    debugDot(tu, tdec, newcolumn);
  }

  TU_CALL( TUtdecnewcolumnFree(tu, &newcolumn) );

  if (*pisGraphic && graph)
  {
    /* Add members and edges for empty rows. */
    if (tdec->numRows < matrix->numRows)
    {
      /* Reallocate if necessary. */
      if (tdec->memRows < matrix->numRows)
      {
        TUreallocBlockArray(tu, &tdec->rowEdges, matrix->numRows);
        tdec->memRows = matrix->numRows;
      }

      /* Add single-edge bonds for each missing row. */
      for (int r = tdec->numRows; r < matrix->numRows; ++r)
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

        TUdbgMsg(8, "New empty row %d is edge %d of member %d.\n", r, edge, member);

        tdec->rowEdges[r].edge = edge;
      }

      tdec->numRows = matrix->numRows;
    }

    TU_CALL( TUtdecToGraph(tu, tdec, graph, true, basis, cobasis, NULL) );
  }

  TU_CALL( TUtdecFree(tu, &tdec) );

  return TU_OKAY;
}
