#define TU_DEBUG /* Uncomment to debug graphic. */
// #define TU_DEBUG_DOT /* Uncomment to write dot files of t-decompositions. */

#include <tu/graphic.h>

#include "env_internal.h"
#include "matrix_internal.h"
#include "one_sum.h"
#include "heap.h"

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#define SWAP_INTS(a, b) \
  do \
  { \
    int tmp = a; \
    a = b; \
    b = tmp; \
  } \
  while (false)


typedef enum
{
  DEC_MEMBER_TYPE_INVALID = 0,
  DEC_MEMBER_TYPE_PARALLEL = 1, /**< \brief Indicate a parallel member, also known as a bond. */
  DEC_MEMBER_TYPE_SERIES = 2,   /**< \brief Indicate a series member, also known as a polygon. */
  DEC_MEMBER_TYPE_RIGID = 3,    /**< \brief Indicate a parallel member, also known as prime. */
  DEC_MEMBER_TYPE_LOOP = 4      /**< \brief Indicate a loop. */
} DEC_MEMBER_TYPE;

static inline
const char* memberTypeString(
  DEC_MEMBER_TYPE type
)
{
  switch(type)
  {
  case DEC_MEMBER_TYPE_PARALLEL:
    return "parallel";
  case DEC_MEMBER_TYPE_SERIES:
    return "series";
  case DEC_MEMBER_TYPE_RIGID:
    return "ridig";
  case DEC_MEMBER_TYPE_LOOP:
    return "loop";
  default:
    return "invalid";
  }
}

typedef int DEC_EDGE;   /**< \brief Type for refering to an decomposition edge. */
typedef int DEC_NODE;   /**< \brief Type for refering to a decomposition node. */
typedef int DEC_MEMBER; /**< \brief Type for refering to a decomposition member. */

typedef struct
{
  DEC_NODE representativeNode;  /**< \brief Next representative of same node towards root, or -1 if root. */
} DEC_NODE_DATA;

typedef struct
{
  int name;               /**< \brief Name of this edge.
                           *
                           * 0, 1, ..., m-1 indicate rows, -1,-2, ..., -n indicate columns,
                           * and for (small) k >= 0, MAX_INT-k and -MAX_INT+k indicate
                           * markers of the parent and to the parent, respectively. */
  DEC_MEMBER member;      /**< \brief Member this edge belongs to or -1 if in free list. */
  DEC_NODE head;          /**< \brief Head node of this edge for a rigid member, -1 otherwise. */
  DEC_NODE tail;          /**< \brief Tail node of this edge for a rigid member, -1 otherwise. */
  DEC_EDGE prev;          /**< \brief Next edge of this member. */
  DEC_EDGE next;          /**< \brief Previous edge of this member. */
  DEC_MEMBER childMember; /**< \brief Child member linked to this edge, or -1. */
} DEC_EDGE_DATA;

typedef struct
{
  DEC_MEMBER_TYPE type;             /**< \brief Type of member. Only valid for representative member. */
  DEC_MEMBER representativeMember;  /**< \brief Representative of member, or -1 if this is a representative member. */
  DEC_MEMBER parentMember;          /**< \brief Parent member of this member or -1 for a root. Only valid for representative member. */
  int numEdges;                     /**< \brief Number of edges. Only valid for representative member. */
  DEC_EDGE markerToParent;          /**< \brief Parent marker edge. Only valid for representative member. */
  DEC_EDGE markerOfParent;          /**< \brief Child marker of parent to which this member is linked. Only valid if root representative. */
  DEC_EDGE firstEdge;               /**< \brief First edge in doubly-linked edge list of this member. */
} DEC_MEMBER_DATA;

typedef struct
{
  DEC_EDGE edge;  /**< \brief Edge corresponding to this row or -1. */
} DEC_ROW_DATA;

typedef struct
{
  DEC_EDGE edge;  /**< \brief Edge corresponding to this column or -1. */
} DEC_COLUMN_DATA;

typedef struct
{
  TU* tu;                       /**< \brief \ref TU environment. */

  int memMembers;               /**< \brief Allocated memory for members. */
  int numMembers;               /**< \brief Number of members. */
  DEC_MEMBER_DATA* members;     /**< \brief Array of members. */

  int memEdges;                 /**< \brief Allocated memory for edges. */
  int numEdges;                 /**< \brief Number of used edges. */
  DEC_EDGE_DATA* edges;         /**< \brief Array of edges. */
  DEC_EDGE firstFreeEdge;       /**< \brief First edge in free list or -1. */

  int memNodes;                 /**< \brief Allocated memory for nodes. */
  int numNodes;                 /**< \brief Number of nodes. */
  DEC_NODE_DATA* nodes;         /**< \brief Array of nodes. */
  DEC_NODE firstFreeNode;       /**< \brief First node in free list or -1. */

  int memRows;                  /**< \brief Allocated memory for \c rowEdges. */
  int numRows;                  /**< \brief Number of rows. */
  DEC_ROW_DATA* rowEdges;       /**< \brief Maps each row to its edge. */

  int memColumns;               /**< \brief Allocated memory for \c columnEdges. */
  int numColumns;               /**< \brief Number of columns. */
  DEC_COLUMN_DATA* columnEdges; /**< \brief Maps each column to its edge. */

  int numMarkerPairs;           /**< \brief Number of marker edge pairs in t-decomposition. */
} DEC;

/**
 * \brief Returns \c true if and only \p member is a representative member.
 */

static inline
bool isRepresentativeMember(
  DEC* dec,         /**< Decomposition. */
  DEC_MEMBER member /**< Member of \p dec. */
)
{
  assert(dec);

  return dec->members[member].representativeMember < 0;
}

/**
 * \brief Returns the representative member of \p member.
 */

static
DEC_MEMBER findMember(
  DEC* dec,         /**< Decomposition. */
  DEC_MEMBER member /**< Member for which the representative shall be returned. */
)
{
  DEC_MEMBER current = member;
  DEC_MEMBER next;
  while ((next = dec->members[current].representativeMember) >= 0)
    current = next;
  DEC_MEMBER root = current;
  current = member;
  while ((next = dec->members[current].representativeMember) >= 0)
  {
    if (next != root)
      dec->members[current].representativeMember = root;
    current = next;
  }
  return root;
}

/**
 * \brief Returns the representative of the parent member of \p member.
 * 
 * Assumes that \p member is a representative member.
 */

static inline
DEC_MEMBER findMemberParent(
  DEC* dec,         /**< Decomposition. */
  DEC_MEMBER member /**< Member whose parent shall be returned. */
)
{
  assert(member >= 0);
  assert(isRepresentativeMember(dec, member));

  DEC_MEMBER someParent = dec->members[member].parentMember;
  if (someParent >= 0)
    return findMember(dec, someParent);
  else
    return -1;
}

/**
 * \brief Returns the representative member of \p edge.
 */

static inline
DEC_MEMBER findEdgeMember(
  DEC* dec,     /**< Decomposition. */
  DEC_EDGE edge /**< Edge. */
)
{
  return findMember(dec, dec->edges[edge].member);
}

/**
 * \brief Checks whether \p dec has consistent edge data.
 *
 * \returns Explanation of inconsistency, or \c NULL if consistent.
 */

static
char* consistencyEdges(
  DEC* dec  /**< Decomposition. */
)
{
  assert(dec);

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    DEC_EDGE edge = dec->members[member].firstEdge;
    int countEdges = 0;
    if (edge >= 0)
    {
      do
      {
        if (edge < 0 || edge >= dec->memEdges)
          return TUconsistencyMessage("edge %d of member %d out of range.", member, edge);
        if (dec->edges[edge].next < 0 || dec->edges[edge].next > dec->memEdges)
          return TUconsistencyMessage("edge %d of member %d has next out of range", member, edge);
        if (dec->edges[dec->edges[edge].next].prev != edge)
          return TUconsistencyMessage("member %d has inconsistent edge list", member);
        if (findEdgeMember(dec, edge) != member)
          return TUconsistencyMessage("edge %d belongs to member %d but is in member %d's edge list.", edge,
            findEdgeMember(dec, edge), member);
        edge = dec->edges[edge].next;
        countEdges++;
      }
      while (edge != dec->members[member].firstEdge);
    }
    if (countEdges != dec->members[member].numEdges)
    {
      return TUconsistencyMessage("member %d has %d edges, but numEdges %d", member, countEdges,
        dec->members[member].numEdges);
    }
  }

  return NULL;
}

/**
 * \brief Checks whether \p dec has consistent member data.
 *
 * \returns Explanation of inconsistency, or \c NULL if consistent.
 */

static
char* consistencyMembers(
  DEC* dec  /**< Decomposition. */
)
{
  assert(dec);

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (dec->members[member].type != DEC_MEMBER_TYPE_PARALLEL
      && dec->members[member].type != DEC_MEMBER_TYPE_RIGID
      && dec->members[member].type != DEC_MEMBER_TYPE_SERIES
      && dec->members[member].type != DEC_MEMBER_TYPE_LOOP)
    {
      return TUconsistencyMessage("member %d has invalid type", member);
    }
  }

  return NULL;
}

/**
 * \brief Checks whether \p dec has consistent node data.
 *
 * \returns Explanation of inconsistency, or \c NULL if consistent.
 */

static
char* consistencyNodes(
  DEC* dec  /**< Decomposition. */
)
{
  assert(dec);

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    bool isRigid = dec->members[member].type == DEC_MEMBER_TYPE_RIGID;
    DEC_EDGE edge = dec->members[member].firstEdge;
    if (edge < 0)
      continue;
    do
    {
      DEC_NODE head = dec->edges[edge].head;
      DEC_NODE tail = dec->edges[edge].tail;
      if (isRigid)
      {
        if (head < 0)
          return TUconsistencyMessage("edge %d of rigid member %d has invalid head node", edge, member);
        if (tail < 0)
          return TUconsistencyMessage("edge %d of rigid member %d has invalid tail node", edge, member);
        if (head >= dec->memNodes)
          return TUconsistencyMessage("edge %d of rigid member %d has head node out of range", edge, member);
        if (tail >= dec->memNodes)
          return TUconsistencyMessage("edge %d of rigid member %d has tail node out of range", edge, member);
      }
      edge = dec->edges[edge].next;
    }
    while (edge != dec->members[member].firstEdge);
  }

  return NULL;
}

/**
 * \brief Checks whether the members of \p dec form a forest.
 *
 * \returns Explanation of inconsistency, or \c NULL if consistent.
 */

static
char* consistencyTree(
  DEC* dec  /**< Decomposition. */
)
{
  assert(dec);

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    int length = 0;
    DEC_MEMBER current;
    for (current = dec->members[member].parentMember; current >= 0; current = dec->members[current].parentMember)
    {
      ++length;
      if (length > dec->numMembers)
        return "infinite member parent loop";
    }
  }

  return NULL;
}

typedef enum
{
  TYPE_CYCLE_CHILD = 1,   /**< Parent marker edge plus path edges form a cycle.
                           *   If graphic, this means we can replace the child by its child marker, adding the latter to
                           *   the list of path edges. */
  TYPE_SINGLE_CHILD = 2,  /**< One node of the parent marker edge is a path end.
                           *   If graphic, this means that one terminal node is in this member or its children. */
  TYPE_DOUBLE_CHILD = 3,  /**< Path edges form two cycles and adding the parent marker edge yields one.
                           *   If graphic, this means that both terminal nodes are in this member or its children. */
  TYPE_ROOT = 4           /**< Root member of reduced decomposition. */
} Type;

/**
 * \brief Additional information specific to a path edge.
 */

typedef struct _PathEdge
{
  DEC_EDGE edge;          /**< \brief The actual edge in the decomposition. */
  struct _PathEdge* next; /**< \brief Next edge of this reduced member, or \c NULL. */
} PathEdge;

/**
 * \brief Additional member information specfic to a given path.
 *
 * @TODO: Maybe add parent reduced member as well, so we don't have to go via the membersToReducedMembers array.
 */

typedef struct _ReducedMember
{
  DEC_MEMBER member;                /**< \brief The member from the decomposition. */
  DEC_MEMBER rootMember;            /**< \brief The root member of this component of the decomposition. */
  int depth;                        /**< \brief Depth of this member in the reduced decomposition. */
  Type type;                        /**< \brief Type of this member. */
  struct _ReducedMember* parent;    /**< \brief Parent in the reduced decomposition. */
  int numChildren;                  /**< \brief Number of children in the reduced decomposition. */
  struct _ReducedMember** children; /**< \brief Children in the reduced decomposition. */
  PathEdge* firstPathEdge;          /**< \brief First edge in linked list of path edges of \p member. */
  DEC_NODE rigidEndNodes[4];        /**< \brief For rigid members, the end nodes of the paths inside the member (or -1). */
} ReducedMember;

/**
 * \brief A component of the reduced decomposition.
 */

typedef struct _ReducedComponent
{
  int rootDepth;                /**< \brief Depth of reduced root member. */
  ReducedMember* root;          /**< \brief Reduced root member. */
  DEC_NODE terminalNode[2];     /**< \brief Terminal nodes of path. */
  DEC_MEMBER terminalMember[2]; /**< \brief Terminal members of path. */
  int numTerminals;
} ReducedComponent;

/**
 * \brief Information for adding a new column.
 */

typedef struct
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
  int memNodesDegree;                       /**< \brief Allocated memory for \c nodesDegree. */

  bool* edgesInPath;                        /**< \brief Map from edges to indicator for being in the path. */
  int memEdgesInPath;                       /**< \brief Allocated memory for \p edgesInPath. */
} DEC_NEWCOLUMN;

/**
 * \brief Checks whether \p dec has consistent parent/child structure of members.
 *
 * \returns Explanation of inconsistency, or \c NULL.
 */

static
char* consistencyParentChild(
  DEC* dec  /**< Decomposition. */
)
{
  assert(dec);

  if (dec->memMembers < dec->numMembers)
    return TUconsistencyMessage("member count and memory are inconsistent");
  if (dec->numMembers < 0)
    return TUconsistencyMessage("negative member count");

  int* countChildren = NULL;
  if (TUallocStackArray(dec->tu, &countChildren, dec->memMembers) != TU_OKAY)
    return TUconsistencyMessage("stack allocation in consistencyParentChild() failed");
  for (int m = 0; m < dec->memMembers; ++m)
    countChildren[m] = 0;

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    if (dec->members[member].parentMember >= dec->memMembers)
    {
      TUfreeStackArray(dec->tu, &countChildren);
      return TUconsistencyMessage("parent member of %d is out of range", member);
    }
    if (dec->members[member].parentMember >= 0)
      countChildren[dec->members[member].parentMember]++;
  }

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    DEC_EDGE edge = dec->members[member].firstEdge;
    if (edge < 0)
      continue;
    do
    {
      if (dec->edges[edge].childMember >= 0)
      {
        countChildren[member]--;

        if (findMember(dec, dec->members[findMember(dec, dec->edges[edge].childMember)].parentMember) != findMember(dec, member))
        {
          TUfreeStackArray(dec->tu, &countChildren);
          return TUconsistencyMessage("member %d has child edge %d for child %d whose parent member is %d",
            member, edge, findMember(dec, dec->edges[edge].childMember),
            findMember(dec, dec->members[findMember(dec, dec->edges[edge].childMember)].parentMember));
        }
        if (dec->members[findMember(dec, dec->edges[edge].childMember)].markerOfParent != edge)
        {
          TUfreeStackArray(dec->tu, &countChildren);
          return TUconsistencyMessage("member %d has child edge %d for child %d whose parent's markerOfParent is %d",
            member, edge, findMember(dec, dec->edges[edge].childMember),
            dec->members[findMember(dec, dec->edges[edge].childMember)].markerOfParent);
        }
        DEC_EDGE markerChild = dec->members[findMember(dec, dec->edges[edge].childMember)].markerToParent;
        if (dec->edges[markerChild].name != -dec->edges[edge].name)
        {
          TUfreeStackArray(dec->tu, &countChildren);
          return TUconsistencyMessage("marker edges %d and %d of members %d (parent) and %d (child) have names %d and %d.",
            edge, markerChild, member, findEdgeMember(dec, markerChild), dec->edges[edge].name,
            dec->edges[markerChild].name);
        }
      }
      edge = dec->edges[edge].next;
    }
    while (edge != dec->members[member].firstEdge);
  }

  if (TUfreeStackArray(dec->tu, &countChildren) != TU_OKAY)
    return "stack deallocation in consistencyParentChild() failed";

  return NULL;
}

/**
 * \brief Checks whether \p dec is consistent.
 *
 * \returns Explanation of inconsistency, or \c NULL.
 */

static
char* decConsistency(
  DEC* dec  /**< Decomposition. */
)
{
  char* message = NULL;
  if ((message = consistencyMembers(dec)))
    return message;
  if ((message = consistencyEdges(dec)))
    return message;
  if ((message = consistencyNodes(dec)))
    return message;
  if ((message = consistencyParentChild(dec)))
    return message;
  if ((message = consistencyTree(dec)))
    return message;

  return NULL;
}

/**
 * \brief Returns the representative node of \p node.
 * 
 * Assumes \p node is indeed a node, i.e., not -1.
 */

static
DEC_NODE findNode(
  DEC* dec,
  DEC_NODE node
)
{
  assert(node >= 0);

  DEC_NODE current = node;
  DEC_NODE next;
  while ((next = dec->nodes[current].representativeNode) >= 0)
    current = next;
  DEC_NODE root = current;
  current = node;
  while ((next = dec->nodes[current].representativeNode) >= 0)
  {
    if (next != root)
      dec->nodes[current].representativeNode = root;
    current = next;
  }
  return root;
}

/**
 * \brief Returns the representative node of the tail of \p edge.
 */

static inline
DEC_NODE findEdgeTail(
  DEC* dec,     /**< Decomposition. */
  DEC_EDGE edge /**< Edge. */
)
{
  assert(dec);
  assert(edge >= 0);
  assert(edge < dec->memEdges);
  assert(dec->edges[edge].tail >= 0);
  assert(dec->edges[edge].tail < dec->memNodes);

  return findNode(dec, dec->edges[edge].tail);
}

/**
 * \brief Returns the representative node of the head of \p edge.
 */

static inline
DEC_NODE findEdgeHead(
  DEC* dec,     /**< Decomposition. */
  DEC_EDGE edge /**< Edge. */
)
{
  assert(dec);
  assert(edge >= 0);
  assert(edge < dec->memEdges);
  assert(dec->edges[edge].head >= 0);
  assert(dec->edges[edge].head < dec->memNodes);

  return findNode(dec, dec->edges[edge].head);
}

/**
 * \brief Creates a node for some rigid member of the decomposition \p dec.
 */

static
TU_ERROR createNode(
  DEC* dec,       /**< Decomposition. */
  DEC_NODE* pnode /**< Pointer for storing new node. */
)
{
  assert(dec);
  assert(pnode);

  DEC_NODE node = dec->firstFreeNode;
  if (node >= 0)
  {
    TUdbgMsg(10, "createNode returns free node %d.\n", node);
    dec->firstFreeNode = dec->nodes[node].representativeNode;
  }
  else
  {
    /* No node in free list, so we enlarge the array. */

    int newSize = 2 * dec->memNodes + 16;
    TU_CALL( TUreallocBlockArray(dec->tu, &dec->nodes, newSize) );
    for (int v = dec->memNodes + 1; v < newSize; ++v)
      dec->nodes[v].representativeNode = v+1;
    dec->nodes[newSize-1].representativeNode = -1;
    dec->firstFreeNode = dec->memNodes + 1;
    node = dec->memNodes;
    dec->memNodes = newSize;
    TUdbgMsg(12, "createNode enlarges node array to %d and returns node %d.\n", newSize, node);
  }
  dec->nodes[node].representativeNode = -1;
  dec->numNodes++;

  *pnode = node;

  return TU_OKAY;
}

/**
 * \brief Adds \p edge to the edge list of its member.
 */

static
TU_ERROR addEdgeToMembersEdgeList(
  DEC* dec,     /**< Decomposition. */
  DEC_EDGE edge /**< Edge to be added. */
)
{
  assert(dec);
  assert(edge >= 0);

  DEC_MEMBER member = findEdgeMember(dec, edge);
  DEC_EDGE first = dec->members[member].firstEdge;
  if (first >= 0)
  {
    assert(dec->members[member].numEdges > 0);
    DEC_EDGE last = dec->edges[first].prev;
    dec->edges[edge].next = first;
    dec->edges[edge].prev = last;
    dec->edges[first].prev = edge;
    dec->edges[last].next =  edge;
  }
  else
  {
    assert(dec->members[member].numEdges == 0);
    dec->edges[edge].next = edge;
    dec->edges[edge].prev = edge;
  }
  dec->members[member].firstEdge = edge;
  dec->members[member].numEdges++;

  return TU_OKAY;
}

/**
 * \brief Removes \p edge from the edge list of its member.
 */

static
TU_ERROR removeEdgeFromMembersEdgeList(
  DEC* dec,     /**< Decomposition. */
  DEC_EDGE edge /**< Edge to be added. */
)
{
  assert(dec);
  assert(edge >= 0);

  DEC_MEMBER member = findEdgeMember(dec, edge);
  if (dec->members[member].numEdges == 1)
    dec->members[member].firstEdge = -1;
  else
  {
    if (dec->members[member].firstEdge == edge)
      dec->members[member].firstEdge = dec->edges[edge].next;

    assert(dec->members[member].firstEdge != edge);

    dec->edges[dec->edges[edge].prev].next = dec->edges[edge].next;
    dec->edges[dec->edges[edge].next].prev = dec->edges[edge].prev;
  }

  dec->members[member].numEdges--;

  return TU_OKAY;
}


/**
 * \brief Replaced \p oldEdge in its member by \p newEdge.
 * 
 * Assumes that \p newEdge has the same member but is not in its edge list.
 */

static
TU_ERROR replaceEdgeInMembersEdgeList(
  DEC* dec,         /**< Decomposition. */
  DEC_EDGE oldEdge, /**< Edge to be removed. */
  DEC_EDGE newEdge  /**< Edge to be added. */
)
{
  assert(dec);
  assert(oldEdge >= 0);
  assert(newEdge >= 0);

  DEC_MEMBER member = findEdgeMember(dec, oldEdge);
  assert(findEdgeMember(dec, newEdge) == member);

  dec->edges[newEdge].tail = dec->edges[oldEdge].tail;
  dec->edges[newEdge].head = dec->edges[oldEdge].head;
  dec->edges[newEdge].next = dec->edges[oldEdge].next;
  dec->edges[newEdge].prev = dec->edges[oldEdge].prev;
  dec->edges[dec->edges[oldEdge].next].prev = newEdge;
  dec->edges[dec->edges[oldEdge].prev].next = newEdge;
  if (dec->members[member].firstEdge == oldEdge)
    dec->members[member].firstEdge = newEdge;

  return TU_OKAY;
}

/**
 * \brief Creates a new edge in \p member.
 */

static
TU_ERROR createEdge(
  DEC* dec,           /**< Decomposition. */
  DEC_MEMBER member,  /**< Member this edge belongs to. */
  DEC_EDGE* pedge     /**< Pointer for storing the new edge. */
)
{
  assert(dec);
  assert(pedge);
  assert(member < 0 || isRepresentativeMember(dec, member));

  DEC_EDGE edge = dec->firstFreeEdge;
  if (edge >= 0)
  {
    TUdbgMsg(12, "Creating edge %d by using a free edge.\n", edge);
    dec->firstFreeEdge = dec->edges[edge].next;
  }
  else /* No edge in free list, so we enlarge the array. */
  {
    int newSize = 2 * dec->memEdges + 16;
    TU_CALL( TUreallocBlockArray(dec->tu, &dec->edges, newSize) );
    for (int e = dec->memEdges + 1; e < newSize; ++e)
    {
      dec->edges[e].next = e+1;
      dec->edges[e].member = -1;
    }
    dec->edges[newSize-1].next = -1;
    dec->firstFreeEdge = dec->memEdges + 1;
    edge = dec->memEdges;
    dec->memEdges = newSize;
    TUdbgMsg(12, "Creating edge %d and reallocating edge to array %d elements.\n", edge, newSize);
  }

  dec->edges[edge].tail = -1;
  dec->edges[edge].head = -1;
  dec->edges[edge].name = INT_MAX/2;
  dec->edges[edge].member = member;
  dec->numEdges++;

  *pedge = edge;

  return TU_OKAY;
}

/**
 * \brief Creates a pair of marker edges for linking parent \p parentMember to child \p childMember.
 */

static
TU_ERROR createMarkerEdgePair(
  DEC* dec,                     /**< Decomposition. */
  DEC_MEMBER parentMember,      /**< Parent member. */
  DEC_EDGE* pMarkerOfParent,    /**< Pointer for storing the new child marker in the parent. */
  DEC_NODE markerOfParentTail,  /**< Tail node of this *\p pChildMarker. */
  DEC_NODE markerOfParentHead,  /**< Head node of this *\p pChildMarker. */
  DEC_MEMBER childMember,       /**< Child member. */
  DEC_EDGE* pMarkerToParent,    /**< Pointer for storing the new parent marker in the child. */
  DEC_NODE markerToParentTail,  /**< Tail node of this *\p pParentMarker. */
  DEC_NODE markerToParentHead   /**< Head node of this *\p pParentMarker. */
)
{
  assert(dec);
  assert(pMarkerOfParent);
  assert(pMarkerToParent);
  assert(isRepresentativeMember(dec, parentMember));
  assert(isRepresentativeMember(dec, childMember));

  /* Create the child marker edge of the parent member. */

  TU_CALL( createEdge(dec, parentMember, pMarkerOfParent) );
  DEC_EDGE_DATA* data = &dec->edges[*pMarkerOfParent];
  data->tail = markerOfParentTail;
  data->head = markerOfParentHead;
  data->childMember = childMember;
  data->name = -INT_MAX + dec->numMarkerPairs;
  TUdbgMsg(12, "Created child marker edge {%d,%d} <%d> of parent member %d.\n", markerOfParentTail, markerOfParentHead,
    data->name, parentMember);

  /* Create the parent marker edge of the child member. */

  TU_CALL( createEdge(dec, childMember, pMarkerToParent) );
  data = &dec->edges[*pMarkerToParent];
  data->tail = markerToParentTail;
  data->head = markerToParentHead;
  data->childMember = -1;
  dec->members[childMember].parentMember = parentMember;
  dec->members[childMember].markerOfParent = *pMarkerOfParent;
  dec->members[childMember].markerToParent = *pMarkerToParent;
  data->name = INT_MAX - dec->numMarkerPairs;
  TUdbgMsg(12, "Created child marker edge {%d,%d} <%d> of child member %d.\n", markerToParentTail, markerToParentHead,
    data->name, childMember);

  /* Increase counter of used marker pairs. */
  dec->numMarkerPairs++;

  return TU_OKAY;
}

/**
 * \brief Creates a new member.
 */

static
TU_ERROR createMember(
  DEC* dec,             /**< Decomposition. */
  DEC_MEMBER_TYPE type, /**< Type of member. */
  DEC_MEMBER* pmember   /**< Pointer for storing the new member. */
)
{
  assert(dec);

  if (dec->numMembers == dec->memMembers)
  {
    dec->memMembers = 16 + 2 * dec->memMembers;
    TU_CALL( TUreallocBlockArray(dec->tu, &dec->members, dec->memMembers) );
  }

  DEC_MEMBER_DATA* data = &dec->members[dec->numMembers];
  data->markerOfParent = -1;
  data->markerToParent = -1;
  data->firstEdge = -1;
  data->representativeMember = -1;
  data->numEdges = 0;
  data->parentMember = -1;
  data->type = type;
  *pmember = dec->numMembers;
  dec->numMembers++;

  TUdbgMsg(10, "Creating %s member %d.\n", memberTypeString(type), *pmember);

  return TU_OKAY;
}

TU_ERROR decCreate(
  TU* tu,           /**< TU environment. */
  DEC** pdec,  /**< Pointer to new t-decomposition. .*/
  int memEdges,     /**< Initial memory for edges of the t-decomposition. */
  int memNodes,     /**< Initial memory for nodes of the t-decomposition. */
  int memMembers,   /**< Initial memory for members of the t-decomposition. */
  int memRows,      /**< Initial memory for rows. */
  int memColumns    /**< Initial memory for columns. */
)
{
  assert(tu);
  assert(pdec);
  assert(!*pdec);

  TU_CALL( TUallocBlock(tu, pdec) );
  DEC* dec = *pdec;
  dec->tu = tu;
  dec->memMembers = memMembers;
  dec->numMembers = 0;
  dec->members = NULL;
  TU_CALL( TUallocBlockArray(tu, &dec->members, dec->memMembers) );

  if (memNodes < 1)
    memNodes = 1;
  dec->memNodes = memNodes;
  dec->nodes = NULL;
  TU_CALL( TUallocBlockArray(tu, &dec->nodes, memNodes) );
  dec->numNodes = 0;
  for (int v = 0; v < memNodes; ++v)
    dec->nodes[v].representativeNode = v+1;
  dec->nodes[memNodes-1].representativeNode = -1;
  dec->firstFreeNode = 0;

  if (memEdges < 1)
    memEdges = 1;
  dec->memEdges = memEdges;
  dec->edges = NULL;
  TU_CALL( TUallocBlockArray(tu, &dec->edges, memEdges) );
  dec->numEdges = 0;
  dec->numMarkerPairs = 0;

  /* Initialize free list with unused edges. */
  if (memEdges > dec->numEdges)
  {
    for (int e = dec->numEdges; e < memEdges; ++e)
    {
      dec->edges[e].next = e+1;
      dec->edges[e].member = -1;
    }
    dec->edges[memEdges-1].next = -1;
    dec->firstFreeEdge = dec->numEdges;
  }
  else
    dec->firstFreeEdge = -1;

  dec->numRows = 0;
  dec->memRows = memRows;
  dec->rowEdges = NULL;
  TU_CALL( TUallocBlockArray(tu, &dec->rowEdges, dec->memRows) );
  for (int r = 0; r < dec->numRows; ++r)
    dec->rowEdges[r].edge = -1;

  dec->numColumns = 0;
  dec->memColumns = memColumns;
  dec->columnEdges = NULL;
  TU_CALL( TUallocBlockArray(tu, &dec->columnEdges, dec->memColumns) );
  for (int c = 0; c < dec->numColumns; ++c)
    dec->columnEdges[c].edge = -1;

  TUconsistencyAssert( decConsistency(dec) );

  return TU_OKAY;
}

TU_ERROR decFree(
  TU* tu,     /**< \ref TU environment. */
  DEC** pdec /**< Pointer to t-decomposition. .*/
)
{
  assert(pdec);
  assert(*pdec);

  DEC* dec = *pdec;
  TU_CALL( TUfreeBlockArray(tu, &dec->members) );
  TU_CALL( TUfreeBlockArray(tu, &dec->edges) );
  TU_CALL( TUfreeBlockArray(tu, &dec->nodes) );
  TU_CALL( TUfreeBlockArray(tu, &dec->rowEdges) );
  TU_CALL( TUfreeBlockArray(tu, &dec->columnEdges) );
  TU_CALL( TUfreeBlock(tu, pdec) );

  return TU_OKAY;
}

/**
 * \brief Creates a graph represented by given t-decomposition.
 */

TU_ERROR decToGraph(
  TU* tu,                 /**< \ref TU environment. */
  DEC* dec,              /**< t-decomposition. */
  TU_GRAPH* graph,        /**< Graph. */
  bool merge,             /**< Merge and remove corresponding parent and child markers. */
  TU_GRAPH_EDGE* basis,   /**< If not NULL, the edges of a spanning tree are stored here. */
  TU_GRAPH_EDGE* cobasis, /**< If not NULL, the non-basis edges are stored here. */
  int* edgeElements       /**< If not NULL, the elements for each edge are stored here. */
)
{
  assert(tu);
  assert(dec);
  assert(graph);

  TUconsistencyAssert( decConsistency(dec) );

  TUdbgMsg(0, "TUdecToGraph for t-decomposition.\n");

  TU_CALL( TUgraphClear(tu, graph) );

  TU_GRAPH_EDGE* localEdgeElements = NULL;
  if (edgeElements)
    localEdgeElements = edgeElements;
  else if (basis || cobasis)
    TU_CALL( TUallocStackArray(tu, &localEdgeElements, dec->memEdges) );
  TU_GRAPH_NODE* decNodesToGraphNodes = NULL;
  TU_CALL( TUallocStackArray(tu, &decNodesToGraphNodes, dec->memNodes) );
  TU_GRAPH_EDGE* decEdgesToGraphEdges = NULL;
  TU_CALL( TUallocStackArray(tu, &decEdgesToGraphEdges, dec->memEdges) );

  for (int v = 0; v < dec->memNodes; ++v)
  {
    if (dec->nodes[v].representativeNode < 0)
    {
      TU_CALL( TUgraphAddNode(tu, graph, &decNodesToGraphNodes[v]) );
    }
    else
      decNodesToGraphNodes[v] = -1;
  }

  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    DEC_MEMBER_TYPE type = dec->members[member].type;
    TUdbgMsg(2, "Member %d is %s with %d edges.\n", member, memberTypeString(type), dec->members[member].numEdges);

    TU_GRAPH_EDGE graphEdge;
    DEC_EDGE edge = dec->members[member].firstEdge;
    if (type == DEC_MEMBER_TYPE_RIGID)
    {
      do
      {
        DEC_NODE head = findEdgeHead(dec, edge);
        DEC_NODE tail = findEdgeTail(dec, edge);
        TU_CALL( TUgraphAddEdge(tu, graph, decNodesToGraphNodes[head], decNodesToGraphNodes[tail],
          &graphEdge) );
        decEdgesToGraphEdges[edge] = graphEdge;
        if (localEdgeElements)
          localEdgeElements[graphEdge] = dec->edges[edge].name;
        edge = dec->edges[edge].next;
      }
      while (edge != dec->members[member].firstEdge);
    }
    else if (type == DEC_MEMBER_TYPE_PARALLEL)
    {
      TU_GRAPH_NODE graphHead, graphTail;
      TU_CALL( TUgraphAddNode(tu, graph, &graphHead) );
      TU_CALL( TUgraphAddNode(tu, graph, &graphTail) );
      do
      {
        TU_CALL( TUgraphAddEdge(tu, graph, graphHead, graphTail, &graphEdge) );
        decEdgesToGraphEdges[edge] = graphEdge;
        if (localEdgeElements)
          localEdgeElements[graphEdge] = dec->edges[edge].name;
        edge = dec->edges[edge].next;
      }
      while (edge != dec->members[member].firstEdge);
    }
    else if (type == DEC_MEMBER_TYPE_SERIES)
    {
      TU_GRAPH_NODE firstNode, v;
      TU_CALL( TUgraphAddNode(tu, graph, &firstNode) );
      v = firstNode;
      edge = dec->edges[edge].next;
      while (edge != dec->members[member].firstEdge)
      {
        TU_GRAPH_NODE w;
        TU_CALL( TUgraphAddNode(tu, graph, &w) );
        TU_CALL( TUgraphAddEdge(tu, graph, v, w, &graphEdge) );
        decEdgesToGraphEdges[edge] = graphEdge;
        if (localEdgeElements)
          localEdgeElements[graphEdge] = dec->edges[edge].name;

        edge = dec->edges[edge].next;
        v = w;
      }
      TU_CALL( TUgraphAddEdge(tu, graph, v, firstNode, &graphEdge) );
      decEdgesToGraphEdges[edge] = graphEdge;
      if (localEdgeElements)
        localEdgeElements[graphEdge] = dec->edges[edge].name;
    }
    else
    {
      assert(type == DEC_MEMBER_TYPE_LOOP);

      TU_GRAPH_NODE v;
      TU_CALL( TUgraphAddNode(tu, graph, &v) );
      TU_CALL( TUgraphAddEdge(tu, graph, v, v, &graphEdge) );
      decEdgesToGraphEdges[edge] = graphEdge;
      if (localEdgeElements)
        localEdgeElements[graphEdge] = dec->edges[edge].name;
    }
  }

  /* Merge respective parent and child edges. */

  if (merge)
  {
    TUdbgMsg(2, "Before merging, the graph has %d nodes and %d edges.\n", TUgraphNumNodes(graph),
      TUgraphNumEdges(graph));

    for (int m = 0; m < dec->numMembers; ++m)
    {
      if (!isRepresentativeMember(dec, m) || dec->members[m].parentMember < 0)
        continue;

      TU_GRAPH_EDGE parent = decEdgesToGraphEdges[dec->members[m].markerOfParent];
      TU_GRAPH_EDGE child = decEdgesToGraphEdges[dec->members[m].markerToParent];
      TU_GRAPH_NODE parentU = TUgraphEdgeU(graph, parent);
      TU_GRAPH_NODE parentV = TUgraphEdgeV(graph, parent);
      TU_GRAPH_NODE childU = TUgraphEdgeU(graph, child);
      TU_GRAPH_NODE childV = TUgraphEdgeV(graph, child);

      TUdbgMsg(2, "Merging edges %d = {%d,%d} <%d> and %d = {%d,%d} <%d>.\n", parent, parentU, parentV,
        dec->edges[dec->members[m].markerOfParent].name, child, childU, childV,
        dec->edges[dec->members[m].markerToParent].name);

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
    for (int r = 0; r < dec->numRows; ++r)
      basis[r] = INT_MIN;
    for (int c = 0; c < dec->numColumns; ++c)
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
    for (int r = 0; r < dec->numRows; ++r)
      assert(basis[r] >= 0);
    for (int c = 0; c < dec->numColumns; ++c)
      assert(cobasis[c] >= 0);
#endif /* !NDEBUG */
  }

  TU_CALL( TUfreeStackArray(tu, &decEdgesToGraphEdges) );
  TU_CALL( TUfreeStackArray(tu, &decNodesToGraphNodes) );
  if (localEdgeElements != edgeElements)
    TU_CALL( TUfreeStackArray(tu, &localEdgeElements) );

  return TU_OKAY;
}

static
void edgeToDot(
  FILE* stream,           /**< File stream. */
  DEC* dec,          /**< t-decomposition. */
  DEC_MEMBER member,  /**< Member this edge belongs to. */
  DEC_EDGE edge,      /**< Edge. */
  int u,                  /**< First node. */
  int v,                  /**< Second node. */
  bool red                /**< Whether to color it red. */
)
{
  assert(stream);
  assert(member >= 0);
  assert(edge >= 0);

  member = findMember(dec, member);

  char type = (dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL) ?
    'P' : (dec->members[member].type == DEC_MEMBER_TYPE_SERIES ? 'S' : 'R');
  const char* redStyle = red ? ",color=red" : "";
  if (dec->members[member].markerToParent == edge)
  {
    fprintf(stream, "    %c_%d_%d -- %c_p_%d [label=\"%d\",style=dashed%s];\n", type, member, u, type, member, edge, redStyle);
    fprintf(stream, "    %c_p_%d -- %c_%d_%d [label=\"%d\",style=dashed%s];\n", type, member, type, member, v, edge, redStyle);
    fprintf(stream, "    %c_%d_%d [shape=box];\n", type, member, u);
    fprintf(stream, "    %c_%d_%d [shape=box];\n", type, member, v);
    fprintf(stream, "    %c_p_%d [style=dashed];\n", type, member);
  }
  else if (dec->edges[edge].childMember >= 0)
  {
    DEC_MEMBER child = findMember(dec, dec->edges[edge].childMember);
    char childType = (dec->members[child].type == DEC_MEMBER_TYPE_PARALLEL) ?
      'P' : (dec->members[child].type == DEC_MEMBER_TYPE_SERIES ? 'S' : 'R');
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
      edge, dec->edges[edge].name, redStyle);
    fprintf(stream, "    %c_%d_%d [shape=box];\n", type, member, u);
    fprintf(stream, "    %c_%d_%d [shape=box];\n", type, member, v);
  }
}

/**
 * \brief Visualizes a t-decomposition as a graph in \c dot format.
 */

TU_ERROR TUdecToDot(
  TU* tu,                 /**< \ref TU environment. */
  DEC* dec,              /**< t-decomposition. */
  FILE* stream,           /**< Stream to write to. */
  bool* edgesHighlighted  /**< Indicator for edges to be highlighted. */
)
{
  assert(tu);
  assert(dec);
  assert(stream);

  fprintf(stream, "// t-decomposition\n");
  fprintf(stream, "graph dec {\n");
  fprintf(stream, "  compound = true;\n");
  for (DEC_MEMBER member = 0; member < dec->numMembers; ++member)
  {
    if (!isRepresentativeMember(dec, member))
      continue;

    fprintf(stream, "  subgraph member%d {\n", member);
    if (dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL)
    {
      DEC_EDGE edge = dec->members[member].firstEdge;
      do
      {
        edgeToDot(stream, dec, member, edge, 0, 1, edgesHighlighted ? edgesHighlighted[edge] : false);
        edge = dec->edges[edge].next;
      }
      while (edge != dec->members[member].firstEdge);
    }
    else if (dec->members[member].type == DEC_MEMBER_TYPE_RIGID)
    {
      DEC_EDGE edge = dec->members[member].firstEdge;
      do
      {
        DEC_NODE u = findEdgeHead(dec, edge);
        DEC_NODE v = findEdgeTail(dec, edge);
        edgeToDot(stream, dec, member, edge, u, v,
          edgesHighlighted ? edgesHighlighted[edge] : false);
        edge = dec->edges[edge].next;
      }
      while (edge != dec->members[member].firstEdge);
    }
    else if (dec->members[member].type == DEC_MEMBER_TYPE_SERIES)
    {
      DEC_EDGE edge = dec->members[member].firstEdge;
      int i = 0;
      do
      {
        edgeToDot(stream, dec, member, edge, i, (i+1) % dec->members[member].numEdges,
          edgesHighlighted ? edgesHighlighted[edge] : false);
        edge = dec->edges[edge].next;
        i++;
      }
      while (edge != dec->members[member].firstEdge);
    }
    else
    {
      assert(dec->members[member].type == DEC_MEMBER_TYPE_LOOP);

      edgeToDot(stream, dec, member, dec->members[member].firstEdge, 0, 0, false);
    }
    fprintf(stream, "  }\n");
  }
  fprintf(stream, "}\n");

  return TU_OKAY;
}

static inline
void debugDot(
  TU* tu,                       /**< \ref TU environment. */
  DEC* dec,                /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn  /**< new column. */
)
{
#if defined(TU_DEBUG_DOT)
  static int dotFileCounter = 1;
  char name[256];
  snprintf(name, 256, "dec-%03d.dot", dotFileCounter);
  TUdbgMsg(0, "Writing <%s>...", name);
  FILE* dotFile = fopen(name, "w");
  TU_CALL( TUdecToDot(tu, dec, dotFile, newcolumn ? newcolumn->edgesInPath : NULL) );
  fclose(dotFile);
  TUdbgMsg(0, " done.\n");

  dotFileCounter++;
#endif /* TU_DEBUG_DOT */
}

TU_ERROR newcolumnCreate(
  TU* tu,                         /**< TU environment. */
  DEC_NEWCOLUMN** pnewcolumn  /**< new-column structure. */
)
{
  assert(tu);

  TU_CALL( TUallocBlock(tu, pnewcolumn) );
  DEC_NEWCOLUMN* newcolumn = *pnewcolumn;
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
  newcolumn->memNodesDegree = 0;
  
  newcolumn->edgesInPath = NULL;
  newcolumn->memEdgesInPath = 0;

  return TU_OKAY;
}

TU_ERROR newcolumnFree(
  TU* tu,                         /**< TU environment. */
  DEC_NEWCOLUMN** pnewcolumn  /**< new-column structure. */
)
{
  assert(tu);
  assert(*pnewcolumn);

  DEC_NEWCOLUMN* newcolumn = *pnewcolumn;

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
  DEC* dec,                /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn  /**< new column. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);
  assert(newcolumn->numReducedMembers == 0);

  newcolumn->remainsGraphic = true;
  newcolumn->numPathEdges = 0;

  /* memEdges does not suffice since new edges can be created by squeezing off.
   * Each squeezing off introduces 4 new edges, and we might apply this twice for each series member. */
  size_t requiredNumEdgesInPath = dec->memEdges + 8*dec->numMembers + 32;
  if (requiredNumEdgesInPath > newcolumn->memEdgesInPath)
  {
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->edgesInPath, requiredNumEdgesInPath) );
    newcolumn->memEdgesInPath = requiredNumEdgesInPath;
  }

  if (newcolumn->memNodesDegree < dec->memNodes)
  {
    while (newcolumn->memNodesDegree < dec->memNodes)
      newcolumn->memNodesDegree = 16 + 2 * newcolumn->memNodesDegree;
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->nodesDegree, newcolumn->memNodesDegree) );
  }

#if defined(TU_DEBUG_DOT)
  for (int e = 0; e < requiredNumEdgesInPath; ++e)
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
TU_ERROR parallelParentChildCheckMember(
  TU* tu,                     /**< \ref TU environment. */
  DEC* dec,              /**< t-decomposition. */
  bool* visitedMembers,       /**< Map from members to indicator of whether this was already a \p childMember. */
  DEC_MEMBER member,      /**< Parent of \p childMember. */
  DEC_MEMBER childMember  /**< Child member of \p member. */
)
{
  assert(tu);
  assert(dec);
  assert(childMember >= 0);
  assert(member == findMemberParent(dec, childMember));
  assert(visitedMembers);

  /* Stop if childMember was already processed at some point. */
  if (visitedMembers[childMember])
    return TU_OKAY;

  visitedMembers[childMember] = true;
  DEC_MEMBER parentMember = findMemberParent(dec, member);
  TUdbgMsg(10, "Consider child member %d of %d with parent %d.\n", childMember, member, parentMember);
  if (parentMember < 0)
    return TU_OKAY;

  if (dec->members[member].type == DEC_MEMBER_TYPE_RIGID)
  {
    TUdbgMsg(10, "Checking if the child marker of %d for child member %d is not parallel to %d's parent marker.\n",
      member, childMember, member);

    DEC_EDGE childMarkerEdge = dec->members[childMember].markerOfParent;
    DEC_NODE nodes[4] = {
      findEdgeTail(dec, childMarkerEdge),
      findEdgeHead(dec, childMarkerEdge),
      findEdgeTail(dec, dec->members[member].markerToParent),
      findEdgeHead(dec, dec->members[member].markerToParent),
    };

    if ((nodes[0] == nodes[2] && nodes[1] == nodes[3]) || (nodes[0] == nodes[3] && nodes[1] == nodes[2]))
    {
      TUdbgMsg(12, "Child marker edge %d is parallel to parent marker edge %d.\n",
        dec->members[childMember].markerOfParent, dec->members[member].markerToParent);

      if (dec->members[parentMember].type != DEC_MEMBER_TYPE_PARALLEL)
      {
        TUdbgMsg(12, "Parent member is not paralle, so we create a parallel member between edges %d and %d.\n",
          dec->members[member].markerToParent, dec->members[member].markerOfParent);

        DEC_MEMBER newParallel;
        TU_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &newParallel) );
        DEC_EDGE markerOfParent = dec->members[member].markerOfParent;
        DEC_EDGE newMarkerOfParent, newMarkerToParent;
        TU_CALL( createMarkerEdgePair(dec, parentMember, &newMarkerOfParent, dec->edges[markerOfParent].tail,
          dec->edges[markerOfParent].head, newParallel, &newMarkerToParent, -1, -1) );

        TU_CALL( replaceEdgeInMembersEdgeList(dec, markerOfParent, newMarkerOfParent) );
        dec->edges[markerOfParent].childMember = member;
        dec->edges[markerOfParent].member = newParallel;
        dec->edges[markerOfParent].tail = -1;
        dec->edges[markerOfParent].head = -1;
        TU_CALL( addEdgeToMembersEdgeList(dec, markerOfParent) );
        TU_CALL( addEdgeToMembersEdgeList(dec, newMarkerToParent) );
        dec->members[member].parentMember = newParallel;
        parentMember = newParallel;
      }

      assert(dec->members[parentMember].type == DEC_MEMBER_TYPE_PARALLEL);

      TUdbgMsg(12, "Removing child marker edge %d from member %d and adding it to the new parallel member %d.\n",
        childMarkerEdge, member, parentMember);

      TU_CALL( removeEdgeFromMembersEdgeList(dec, childMarkerEdge) );
      dec->edges[childMarkerEdge].member = parentMember;
      TU_CALL( addEdgeToMembersEdgeList(dec, childMarkerEdge) );
      dec->edges[childMarkerEdge].tail = -1;
      dec->edges[childMarkerEdge].head = -1;
      dec->members[childMember].parentMember = parentMember;

      debugDot(tu, dec, NULL);
    }
  }

  TU_CALL( parallelParentChildCheckMember(tu, dec, visitedMembers, parentMember, member) );

  return TU_OKAY;
}

static
TU_ERROR parallelParentChildCheckReducedMembers(
  TU* tu,         /**< \ref TU environment. */
  DEC* dec,  /**< t-decomposition. */
  int* entryRows, /**< Array of rows of new column's enries. */
  int numEntries  /**< Length of \p entryRows. */
)
{
  assert(tu);
  assert(dec);
  assert(entryRows);

  /* Create array of visited members to avoid multiple consideration of a whole path towards the same root. */

  size_t maxNumMembers = 2*dec->numMembers;
  bool* visitedMembers = NULL;
  TU_CALL( TUallocStackArray(tu, &visitedMembers, maxNumMembers) );
  for (int m = 0; m < maxNumMembers; ++m)
    visitedMembers[m] = false;

  /* Start processing at each member containing a row element. */
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    if (row >= dec->numRows)
      continue;
    DEC_EDGE edge = dec->rowEdges[row].edge;
    if (edge < 0)
      continue;
    DEC_MEMBER member = findEdgeMember(dec, edge);
    TUdbgMsg(6, "Entry %d is row %d of %d and corresponds to edge %d of member %d.\n", p, row, dec->numRows, edge,
      member);
    DEC_MEMBER parentMember = findMemberParent(dec, member);
    if (parentMember >= 0)
      TU_CALL( parallelParentChildCheckMember(tu, dec, visitedMembers, parentMember, member) );
  }

  assert(dec->numMembers <= maxNumMembers);
  TU_CALL( TUfreeStackArray(tu, &visitedMembers) );

  return TU_OKAY;
}

/**
 * \brief Creates, if necessary, the reduced member for \p member and calls itself for the parent.
 */

static
TU_ERROR createReducedMembers(
  TU* tu,                             /**< \ref TU environment. */
  DEC* dec,                      /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn,       /**< new column. */
  DEC_MEMBER member,              /**< Member to create reduced member for. */
  ReducedMember** rootDepthMinimizer, /**< Array mapping root members to the depth minimizer. */
  ReducedMember** pReducedMember      /**< Pointer for storing the created reduced member. */
)
{
  assert(tu);
  assert(dec);
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
    reducedMember->rigidEndNodes[0] = -1;
    reducedMember->rigidEndNodes[1] = -1;
    reducedMember->rigidEndNodes[2] = -1;
    reducedMember->rigidEndNodes[3] = -1;

    TUdbgMsg(10, "Reduced member is new.\n");

    DEC_MEMBER parentMember = findMemberParent(dec, member);
    if (parentMember >= 0)
    {
      /* The parent member might have changed. */
      parentMember = findMemberParent(dec, member);

      ReducedMember* parentReducedMember;
      TU_CALL( createReducedMembers(tu, dec, newcolumn, parentMember, rootDepthMinimizer, &parentReducedMember) );
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
  DEC* dec,                /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn, /**< new column. */
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
  size_t maxNumReducedMembers = dec->numMembers + maxRow + 1;
  if (newcolumn->memReducedMembers < maxNumReducedMembers)
  {
    newcolumn->memReducedMembers = maxNumReducedMembers;
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->reducedMembers, newcolumn->memReducedMembers) );
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->membersToReducedMembers,
      newcolumn->memReducedMembers) );
  }

  /* Initialize the mapping from members to reduced members. */
  for (int m = 0; m < dec->numMembers; ++m)
    newcolumn->membersToReducedMembers[m] = NULL;

  ReducedMember** rootDepthMinimizer = NULL;
  /* Factor 2 because of possible new parallels due to ensureChildParentMarkersNotParallel */
  TU_CALL( TUallocStackArray(tu, &rootDepthMinimizer, 2*dec->numMembers) );
  for (int m = 0; m < 2*dec->numMembers; ++m)
    rootDepthMinimizer[m] = NULL;
  newcolumn->numReducedMembers = 0;
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    DEC_EDGE edge = (row < dec->numRows) ? dec->rowEdges[row].edge : -1;
    TUdbgMsg(6, "Entry %d is row %d of %d and corresponds to edge %d.\n", p, row, dec->numRows, edge);
    if (edge >= 0)
    {
      DEC_MEMBER member = findEdgeMember(dec, edge);
      TUdbgMsg(8, "Edge %d exists and belongs to member %d.\n", edge, member);
      ReducedMember* reducedMember;
      TU_CALL( createReducedMembers(tu, dec, newcolumn, member, rootDepthMinimizer, &reducedMember) );

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

    DEC_MEMBER parentMember = findMemberParent(dec, newcolumn->reducedMembers[m].member);
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
  DEC_NEWCOLUMN* newcolumn, /**< new column. */
  DEC_EDGE edge,            /**< The edge it refers to. */
  PathEdge* next,               /**< The next edge in the singly linked list of this reduced member. */
  PathEdge** pNewPathEdge       /**< Pointer for storing the new path edge. */
)
{
  assert(tu);
  assert(newcolumn);

  assert(newcolumn->numPathEdges < newcolumn->memPathEdges);
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
  DEC* dec,                /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn, /**< new column. */
  int* entryRows,               /**< Array of rows of new column's enries. */
  int numEntries                /**< Length of \p entryRows. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);

  /* Check if we need new rows. */

  int newNumRows = dec->numRows-1;
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    DEC_EDGE edge = (row < dec->numRows) ? dec->rowEdges[row].edge : -1;
    if (edge < 0)
    {
      if (row > newNumRows)
        newNumRows = row;
    }
  }
  newNumRows++;

  TUdbgMsg(4, "Completing reduced decomposition: increasing #rows from %d to %d.\n", dec->numRows, newNumRows);

  /* Create single-edge parallel members for all new rows. */

  if (newNumRows > dec->numRows)
  {
    if (newNumRows > dec->memRows)
    {
      dec->memRows = 2*newNumRows;
      TU_CALL( TUreallocBlockArray(tu, &dec->rowEdges, dec->memRows) );
    }

    for (int r = dec->numRows; r < newNumRows; ++r)
    {
      DEC_MEMBER member;
      TU_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &member) );

      DEC_EDGE edge;
      TU_CALL( createEdge(dec, member, &edge) );
      TU_CALL( addEdgeToMembersEdgeList(dec, edge) );
      dec->edges[edge].name = r;
      dec->edges[edge].head = -1;
      dec->edges[edge].tail = -1;
      dec->edges[edge].childMember = -1;

      TUdbgMsg(8, "New row %d is edge %d of member %d.\n", r, edge, member);

      dec->rowEdges[r].edge = edge;
    }
  }

  /* Create reduced members. */
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    if (row >= dec->numRows)
    {
      /* Edge and member for this row were created. We create the reduced member. */
      DEC_EDGE edge = dec->rowEdges[row].edge;
      DEC_MEMBER member = findEdgeMember(dec, edge);

      TUdbgMsg(4, "Creating reduced member for edge %d of member %d.\n", edge, member);

      assert(newcolumn->numReducedMembers < newcolumn->memReducedMembers);
      ReducedMember* reducedMember = &newcolumn->reducedMembers[newcolumn->numReducedMembers];
      newcolumn->numReducedMembers++;
      reducedMember->numChildren = 0;
      reducedMember->member = member;
      reducedMember->depth = 0;
      reducedMember->rootMember = -1;
      reducedMember->type = TYPE_ROOT;

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

  dec->numRows = newNumRows;

  return TU_OKAY;
}

static
TU_ERROR initializeReducedMemberEdgeLists(
  TU* tu,                       /**< \ref TU environment. */
  DEC* dec,                /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn, /**< new column. */
  int* entryRows,               /**< Array of rows of new column's enries. */
  int numEntries                /**< Length of \p entryRows. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);

  TUdbgMsg(4, "Initializing edge lists for members of reduced t-decomposition.\n");

  for (int v = 0; v < dec->memNodes; ++v)
    newcolumn->nodesDegree[v] = 0;
  for (int e = 0; e < dec->memEdges; ++e)
    newcolumn->edgesInPath[e] = false;

  assert(newcolumn->numPathEdges == 0);
  size_t maxNumPathEdges = numEntries + dec->numMembers;
  if (newcolumn->memPathEdges < maxNumPathEdges)
  {
    TU_CALL( TUreallocBlockArray(tu, &newcolumn->pathEdges, maxNumPathEdges) );
    newcolumn->memPathEdges = maxNumPathEdges;
  }

  /* Start with empty lists. */
  for (int i = 0; i < newcolumn->numReducedMembers; ++i)
    newcolumn->reducedMembers[i].firstPathEdge = NULL;

  /* Fill edge lists. */
  for (int p = 0; p < numEntries; ++p)
  {
    int row = entryRows[p];
    DEC_EDGE edge = (row < dec->numRows) ? dec->rowEdges[row].edge : -1;
    if (edge >= 0)
    {
      DEC_MEMBER member = findEdgeMember(dec, edge);
      assert(member >= 0);
      ReducedMember* reducedMember = newcolumn->membersToReducedMembers[member];
      assert(reducedMember);
      PathEdge* pathEdge;
      TU_CALL( createPathEdge(tu, newcolumn, edge, reducedMember->firstPathEdge, &pathEdge) );
      reducedMember->firstPathEdge = pathEdge;
      newcolumn->edgesInPath[edge] = true;
      if (dec->members[member].type == DEC_MEMBER_TYPE_RIGID)
      {
        newcolumn->nodesDegree[findEdgeHead(dec, edge)]++;
        newcolumn->nodesDegree[findEdgeTail(dec, edge)]++;
      }

      TUdbgMsg(6, "Edge %d <%d> belongs to reduced member %ld which is member %d.\n", edge, dec->edges[edge].name,
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
  DEC* dec,                    /**< t-decomposition. */
  ReducedMember* reducedMember,     /**< Reduced member. */
  int* pNumOneEnd,                  /**< Number of children that (recursively) must contain one path end. */
  int* pNumTwoEnds,                 /**< Number of children that (recursively) must contain two path ends. */
  DEC_EDGE childMarkerEdges[2]  /**< Array for storing a child marker edges containing one/two path ends. */
)
{
  assert(tu);
  assert(dec);
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

    if (child->type == TYPE_SINGLE_CHILD)
    {
      if (nextChildMarker < 2)
      {
        childMarkerEdges[nextChildMarker] = dec->members[findMember(dec, child->member)].markerOfParent;
        nextChildMarker++;
      }
      (*pNumOneEnd)++;
    }
    else if (child->type == TYPE_DOUBLE_CHILD)
    {
      if (nextChildMarker < 2)
      {
        childMarkerEdges[nextChildMarker] = dec->members[findMember(dec, child->member)].markerOfParent;
        nextChildMarker++;
      }
      (*pNumTwoEnds)++;
    }
  }

  return TU_OKAY;
}

static
TU_ERROR determineTypeParallel(
  TU* tu,                             /**< \ref TU environment. */
  DEC* dec,                      /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int numOneEnd,                      /**< Number of child markers containing one path end. */
  int numTwoEnds,                     /**< Number of child markers containing two path ends. */
  DEC_EDGE childMarkerEdges[2],   /**< Marker edges of children containing one/two path ends. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);
  assert(reducedComponent);
  assert(reducedMember);
  assert(numOneEnd >= 0);
  assert(numTwoEnds >= 0);
  assert(numOneEnd + 2*numTwoEnds <= 2);
  assert(childMarkerEdges);
  assert(dec->members[findMember(dec, reducedMember->member)].type == DEC_MEMBER_TYPE_PARALLEL);

  if (depth == 0)
  {
    /* A parallel root always works. */
    reducedMember->type = TYPE_ROOT;
    return TU_OKAY;
  }

  /* No children, but a path edge. */
  if (2*numTwoEnds + numOneEnd == 0 && reducedMember->firstPathEdge)
    reducedMember->type = TYPE_CYCLE_CHILD;
  else if (numOneEnd == 1)
    reducedMember->type = TYPE_SINGLE_CHILD;
  else if (numOneEnd + 2 * numTwoEnds == 2)
  {
    if (reducedMember->firstPathEdge)
    {
      /* Tested in TypingInnerParallelWithEdge */
      newcolumn->remainsGraphic = false;
    }
    else
      reducedMember->type = TYPE_DOUBLE_CHILD;
  }
  else
  {
    /* Since no child contains path edges, the parallel must be a leaf of the reduced decomposition and contains one. */
    assert(reducedMember->firstPathEdge);
    reducedMember->type = TYPE_CYCLE_CHILD;
  }

  return TU_OKAY;
}

static
TU_ERROR determineTypeSeries(
  TU* tu,                             /**< \ref TU environment. */
  DEC* dec,                      /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int numOneEnd,                      /**< Number of child markers containing one path end. */
  int numTwoEnds,                     /**< Number of child markers containing two path ends. */
  DEC_EDGE childMarkerEdges[2],   /**< Marker edges of children containing one/two path ends. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);
  assert(reducedComponent);
  assert(reducedMember);
  assert(numOneEnd >= 0);
  assert(numTwoEnds >= 0);
  assert(numOneEnd + 2*numTwoEnds <= 2);
  assert(childMarkerEdges);
  assert(dec->members[findMember(dec, reducedMember->member)].type == DEC_MEMBER_TYPE_SERIES);

  DEC_MEMBER member = findMember(dec, reducedMember->member);

  int countPathEdges = 0;
  for (PathEdge* edge = reducedMember->firstPathEdge; edge != NULL; edge = edge->next)
    ++countPathEdges;
  int numEdges = dec->members[member].numEdges;

  TUdbgMsg(8+2*depth,
    "Determining type of series with %d edges, %d path edges, %d 1-end children and %d 2-end children.\n", numEdges,
    countPathEdges, numOneEnd, numTwoEnds);

  if (depth == 0)
  {
    /* We assume that we are not the root of the whole decomposition. */
    assert(dec->members[member].parentMember >= 0);

    /* Tested in TypingRootSeriesDoubleChild */
    newcolumn->remainsGraphic = (numTwoEnds == 0);
    reducedMember->type = (countPathEdges == numEdges - 1) ? TYPE_CYCLE_CHILD : TYPE_ROOT;
    return TU_OKAY;
  }

  if (countPathEdges == numEdges - 1)
  {
    reducedMember->type = TYPE_CYCLE_CHILD;
  }
  else if (countPathEdges + numTwoEnds == numEdges - 1)
  {
    assert(numTwoEnds == 1);
    reducedMember->type = TYPE_DOUBLE_CHILD;
  }
  else if (numTwoEnds == 1)
  {
    /* Tested in TypingInnerSeries2Child. */
    newcolumn->remainsGraphic = false;
  }
  else if (numOneEnd == 1)
  {
    reducedMember->type = TYPE_SINGLE_CHILD;
  }
  else if (numOneEnd == 2)
  {
    reducedMember->type = TYPE_DOUBLE_CHILD;
  }
  else
  {
    assert(numOneEnd == 0);
    assert(numTwoEnds == 0);
    reducedMember->type = TYPE_SINGLE_CHILD;
  }

  return TU_OKAY;
}

static
TU_ERROR determineTypeRigid(
  TU* tu,                             /**< \ref TU environment. */
  DEC* dec,                      /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int numOneEnd,                      /**< Number of child markers containing one path end. */
  int numTwoEnds,                     /**< Number of child markers containing two path ends. */
  DEC_EDGE childMarkerEdges[2],   /**< Marker edges of children containing one/two path ends. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);
  assert(reducedComponent);
  assert(reducedMember);
  assert(numOneEnd >= 0);
  assert(numTwoEnds >= 0);
  assert(numOneEnd + 2*numTwoEnds <= 2);
  assert(childMarkerEdges);
  assert(dec->members[findMember(dec, reducedMember->member)].type == DEC_MEMBER_TYPE_RIGID);

  DEC_MEMBER member = findMember(dec, reducedMember->member);

  /* Collect nodes of parent marker and of child markers containing one/two path end nodes. */
  DEC_NODE parentMarkerNodes[2] = {
    depth == 0 ? -1 : findEdgeTail(dec, dec->members[member].markerToParent),
    depth == 0 ? -1 : findEdgeHead(dec, dec->members[member].markerToParent)
  };
  DEC_NODE childMarkerNodes[4] = {
    childMarkerEdges[0] < 0 ? -1 : findEdgeTail(dec, childMarkerEdges[0]),
    childMarkerEdges[0] < 0 ? -1 : findEdgeHead(dec, childMarkerEdges[0]),
    childMarkerEdges[1] < 0 ? -1 : findEdgeTail(dec, childMarkerEdges[1]),
    childMarkerEdges[1] < 0 ? -1 : findEdgeHead(dec, childMarkerEdges[1])
  };

  DEC_NODE* pathEndNodes = reducedMember->rigidEndNodes;
  int numPathEndNodes = 0;

  /* Check the node degrees (with respect to path edges) in this component. */
  for (PathEdge* reducedEdge = reducedMember->firstPathEdge; reducedEdge; reducedEdge = reducedEdge->next)
  {
    DEC_NODE nodes[2] = { findEdgeHead(dec, reducedEdge->edge), findEdgeTail(dec, reducedEdge->edge) };
    for (int i = 0; i < 2; ++i)
    {
      DEC_NODE v = nodes[i];
      if (newcolumn->nodesDegree[v] >= 3)
      {
        TUdbgMsg(6 + 2*depth, "Node %d of rigid member %d has path-degree at least 3.\n", v, reducedMember->member);

        /* Tested in TypingAnyRigidDegree3. */
        newcolumn->remainsGraphic = false;
        return TU_OKAY;
      }

      if (newcolumn->nodesDegree[v] == 1)
      {
        if (numPathEndNodes == 4)
        {
          TUdbgMsg(6 + 2*depth, "Rigid member %d has at least five path end nodes: %d, %d, %d, %d and %d.\n",
            reducedMember->member, pathEndNodes[0], pathEndNodes[1], pathEndNodes[2], pathEndNodes[3], v);

          /* Tested in TypingAnyRigidManyPaths.*/
          newcolumn->remainsGraphic = false;
          return TU_OKAY;
        }

        pathEndNodes[numPathEndNodes] = v;
        ++numPathEndNodes;
      }
    }
  }

  /* Exchange end nodes array using these rules:
   * If there are 4 edges, then one path connects 0 with 1 and the other 2 with 3.
   * For each path, parent marker nodes should come first (i.e., index 0 and 2).
   */

  if (numPathEndNodes == 4)
  {
    DEC_EDGE* nodeEdges = NULL;
    TUallocStackArray(tu, &nodeEdges, 2*dec->memNodes);

    /* Initialize relevant entries to -1. */
    for (PathEdge* reducedEdge = reducedMember->firstPathEdge; reducedEdge; reducedEdge = reducedEdge->next)
    {
      DEC_NODE nodes[2] = { findEdgeHead(dec, reducedEdge->edge), findEdgeTail(dec, reducedEdge->edge) };
      for (int i = 0; i < 2; ++i)
      {
        DEC_NODE v = nodes[i];
        nodeEdges[2*v] = -1;
        nodeEdges[2*v + 1] = -1;
      }
    }

    /* Store incident edges for every node. */
    for (PathEdge* reducedEdge = reducedMember->firstPathEdge; reducedEdge; reducedEdge = reducedEdge->next)
    {
      DEC_NODE nodes[2] = { findEdgeHead(dec, reducedEdge->edge), findEdgeTail(dec, reducedEdge->edge) };
      for (int i = 0; i < 2; ++i)
      {
        DEC_NODE v = nodes[i];
        nodeEdges[2*v + (nodeEdges[2*v] == -1 ? 0 : 1)] = reducedEdge->edge;
      }
    }

    /* Start at end node 0 and see where we end. */
    DEC_EDGE previousEdge = -1;
    DEC_NODE currentNode = pathEndNodes[0];
    while (true)
    {
      DEC_EDGE edge = nodeEdges[2*currentNode];
      if (edge == previousEdge)
        edge = nodeEdges[2*currentNode+1];
      if (edge == -1)
        break;
      previousEdge = edge;
      DEC_NODE v = findEdgeHead(dec, edge);
      currentNode = (v != currentNode) ? v : findEdgeTail(dec, edge);
    }
    TUfreeStackArray(tu, &nodeEdges);

    /* Exchange such that we end nodes 0 and 1 are end nodes of the same path. */
    if (currentNode == pathEndNodes[2])
    {
      pathEndNodes[2] = pathEndNodes[1];
      pathEndNodes[1] = currentNode;
    }
    else if (currentNode == pathEndNodes[3])
    {
      pathEndNodes[3] = pathEndNodes[1];
      pathEndNodes[1] = currentNode;
    }

    /* Exchange such that end node 2 is at the parent marker. */
    if (pathEndNodes[2] != parentMarkerNodes[0] && pathEndNodes[2] != parentMarkerNodes[1])
      SWAP_INTS(pathEndNodes[2], pathEndNodes[3]);
  }

  /* Exchange such that end node 0 is at the parent marker. */
  if (numPathEndNodes >= 2 && pathEndNodes[0] != parentMarkerNodes[0] && pathEndNodes[0] != parentMarkerNodes[1])
  {
    SWAP_INTS(pathEndNodes[0], pathEndNodes[1]);
  }

  TUdbgMsg(6 + 2*depth, "Rigid %smember %d has %d path end nodes:", depth == 0 ? "root " : "", member, numPathEndNodes);
  if (numPathEndNodes >= 2)
    TUdbgMsg(0, " a path from %d to %d.", pathEndNodes[0], pathEndNodes[1]);
  if (numPathEndNodes == 4)
    TUdbgMsg(0, " a path from %d to %d.", pathEndNodes[2], pathEndNodes[3]);
  TUdbgMsg(0, "\n");

  if (depth == 0)
  {
    if (numPathEndNodes == 0)
    {
      TUdbgMsg(6 + 2*depth, "Root without paths. Child markers are {%d,%d} and {%d,%d}.\n", childMarkerNodes[0], childMarkerNodes[1],
        childMarkerNodes[2], childMarkerNodes[3]);
      /* No path edges, so there should be two adjacent child marker edges. */
      if (numOneEnd == 2 && (childMarkerNodes[0] == childMarkerNodes[2] || childMarkerNodes[0] == childMarkerNodes[3]
          || childMarkerNodes[1] == childMarkerNodes[2] || childMarkerNodes[1] == childMarkerNodes[3]))
      {
        reducedMember->type = TYPE_ROOT;
      }
      else
      {
        /* Tested in TypingRootRigidNoPathsDisjointSingleChildren. */
        newcolumn->remainsGraphic = false;
      }
    }
    else if (numPathEndNodes == 2)
    {
      if (numOneEnd == 1)
      {
        bool pathAdjacentToChildMarker = false;
        for (int i = 0; i < 2; ++i)
        {
          for (int j = 0; j < 2; ++j)
          {
            if (pathEndNodes[i] == childMarkerNodes[j])
            {
              pathAdjacentToChildMarker = true;
              break;
            }
          }
        }
        if (pathAdjacentToChildMarker)
          reducedMember->type = TYPE_ROOT;
        else
        {
          /* Tested in TypingRootRigidOnePathSingleChild. */
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
            if (reducedMember->rigidEndNodes[i] == childMarkerNodes[j])
              matched[j/2] = true;
          }
        }

        if (matched[0] && matched[1])
          reducedMember->type = TYPE_ROOT;
        else
        {
          /* Tested in TypingRootRigidOnePathTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
        }
      }
      else if (numTwoEnds == 0)
      {
        assert(numOneEnd == 0);
        reducedMember->type = TYPE_ROOT;
      }
      else
      {
        assert(numOneEnd == 0);
        assert(numTwoEnds == 1);

        if ((childMarkerNodes[0] == pathEndNodes[0] && childMarkerNodes[1] == pathEndNodes[1])
          || (childMarkerNodes[0] == pathEndNodes[1] && childMarkerNodes[1] == pathEndNodes[0]))
        {
          reducedMember->type = TYPE_ROOT;
        }
        else
        {
          /* Tested in TypingRootRigidOnePathDoubleChild. */
          newcolumn->remainsGraphic = false;
        }
      }
    }
    else
    {
      assert(numPathEndNodes == 4);

      /* Tested in TypingRootRigidTwoPaths. */
      newcolumn->remainsGraphic = false;
    }
  }
  else
  {
    /* Non-root rigid member. */

    int parentMarkerDegrees[2] = {
      newcolumn->nodesDegree[parentMarkerNodes[0]],
      newcolumn->nodesDegree[parentMarkerNodes[1]]
    };

    if (numPathEndNodes == 0)
    {
      /* We have no path edges, so there must be at least one child containing one/two path ends. */
      assert(numOneEnd + numTwoEnds > 0);
      /* We should not have a child marker edge parallel to the parent marker edge! */
      assert(!(parentMarkerNodes[0] == childMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[1])
        && !(parentMarkerNodes[0] == childMarkerNodes[1] && parentMarkerNodes[1] == childMarkerNodes[0]));

      if (numOneEnd == 0)
      {
        /* Even if parent and child marker (type 4) are adjacent, this is non-graphic. */

        /* Tested in TypingInnerRigidNoPathNoSingleChild. */
        newcolumn->remainsGraphic = false;
      }
      else if (numOneEnd == 1)
      {
        if (childMarkerNodes[0] == parentMarkerNodes[0] || childMarkerNodes[0] == parentMarkerNodes[1]
          || childMarkerNodes[1] == parentMarkerNodes[0] || childMarkerNodes[1] == parentMarkerNodes[1])
        {
          reducedMember->type = TYPE_SINGLE_CHILD;
        }
        else
        {
          /* Tested in TypingInnerRigidNoPathOneSingleChild. */
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
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          /* Tested in TypingInnerRigidNoPathTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
        }
      }
    }
    else if (numPathEndNodes == 2)
    {
      /* Exchange such that end node 0 is at the parent marker (or none of the end nodes is). */
      if (pathEndNodes[0] != parentMarkerNodes[0] && pathEndNodes[0] != parentMarkerNodes[1])
        SWAP_INTS(pathEndNodes[0], pathEndNodes[0]);

      if (numOneEnd == 1)
      {
        TUdbgMsg(6 + 2*depth, "%d-%d-path, parent {%d,%d} and child {%d,%d}\n", pathEndNodes[0], pathEndNodes[1],
          parentMarkerNodes[0], parentMarkerNodes[1], childMarkerNodes[0], childMarkerNodes[1]);

        if (parentMarkerNodes[0] != pathEndNodes[0])
        {
          SWAP_INTS(parentMarkerNodes[0], parentMarkerNodes[1]);
          SWAP_INTS(parentMarkerDegrees[0], parentMarkerDegrees[1]);
        }
        if (parentMarkerNodes[0] != pathEndNodes[0])
        {
          /* Tested in TypingInnerRigidOnePathOneSingleChild. */
          newcolumn->remainsGraphic = false;
          return TU_OKAY;
        }

        if (parentMarkerNodes[1] == pathEndNodes[1])
        {
          /* Path closes a cycle with parent marker edge. */
          if (childMarkerNodes[0] == parentMarkerNodes[0] || childMarkerNodes[0] == parentMarkerNodes[1]
            || childMarkerNodes[1] == parentMarkerNodes[0] || childMarkerNodes[1] == parentMarkerNodes[1])
          {
            reducedMember->type = TYPE_SINGLE_CHILD;
          }
          else
          {
            /* Tested in TypingInnerRigidOnePathOneSingleChild. */
            newcolumn->remainsGraphic = false;
          }
        }
        else
        {
          /* Path end is not incident to parent marker edge. */
          if (childMarkerNodes[0] == pathEndNodes[1] || childMarkerNodes[1] == pathEndNodes[1])
          {
            reducedMember->type = TYPE_SINGLE_CHILD;
          }
          else if (childMarkerNodes[0] == parentMarkerNodes[1] || childMarkerNodes[1] == parentMarkerNodes[1])
          {
            reducedMember->type = TYPE_DOUBLE_CHILD;
          }
          else
          {
            /* Tested in TypingInnerRigidOnePathOneSingleChild. */
            newcolumn->remainsGraphic = false;
          }
        }
      }
      else if (numOneEnd == 2)
      {
        TUdbgMsg(6 + 2*depth, "%d-%d-path, parent marker {%d,%d} and child markers {%d,%d} and {%d,%d}\n",
          pathEndNodes[0], pathEndNodes[1], parentMarkerNodes[0], parentMarkerNodes[1], childMarkerNodes[0],
          childMarkerNodes[1], childMarkerNodes[2], childMarkerNodes[3]);

        /* We know that pathEndNodes[0] is incident to the parent marker. */
        DEC_NODE otherParentNode;
        if (pathEndNodes[0] == parentMarkerNodes[0])
          otherParentNode = parentMarkerNodes[1];
        else if (pathEndNodes[0] == parentMarkerNodes[1])
          otherParentNode = parentMarkerNodes[0];
        else
        {
          /* Tested in TypingInnerRigidOnePathTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
          return TU_OKAY;
        }

        /* Two 1-ends and a path that closes a cycle with the parent marker is only allowed for root members. */
        if (pathEndNodes[1] == otherParentNode)
        {
          /* Tested in TypingInnerRigidOnePathTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
          return TU_OKAY;
        }

        bool childMatched[2] = { false, false };
        bool pathEndMatched = false;
        bool otherParentMatched = false;
        for (int i = 0; i < 4; ++i)
        {
          if (childMarkerNodes[i] == pathEndNodes[1])
          {
            childMatched[i/2] = true;
            pathEndMatched = true;
          }
          if (childMarkerNodes[i] == otherParentNode)
          {
            childMatched[i/2] = true;
            otherParentMatched = true;
          }
        }

        if (childMatched[0] && childMatched[1] && pathEndMatched && otherParentMatched)
        {
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          /* Tested in TypingInnerRigidOnePathTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
        }
      }
      else if (numTwoEnds == 0)
      {
        if (parentMarkerDegrees[0] % 2 == 0 && parentMarkerDegrees[1] == 1)
        {
          reducedMember->type = TYPE_SINGLE_CHILD;
        }
        else if (parentMarkerDegrees[0] == 1 && parentMarkerDegrees[1] % 2 == 0)
        {
          reducedMember->type = TYPE_SINGLE_CHILD;
        }
        else if (parentMarkerDegrees[0] == 1 && parentMarkerDegrees[1] == 1)
        {
          reducedMember->type = TYPE_CYCLE_CHILD;
        }
        else
        {
          /* Both have degree 0 or 2. */
          /* Tested in TypingInnerRigidOnePathNoChildren. */
          newcolumn->remainsGraphic = false;
        }
      }
      else
      {
        assert(numTwoEnds == 1);

        if ((pathEndNodes[0] == parentMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[0]
          && childMarkerNodes[1] == pathEndNodes[1])
          || (pathEndNodes[0] == parentMarkerNodes[0] && parentMarkerNodes[1] == childMarkerNodes[1]
          && childMarkerNodes[0] == pathEndNodes[1])
          || (pathEndNodes[0] == parentMarkerNodes[1] && parentMarkerNodes[0] == childMarkerNodes[0]
          && childMarkerNodes[1] == pathEndNodes[1])
          || (pathEndNodes[0] == parentMarkerNodes[1] && parentMarkerNodes[0] == childMarkerNodes[1]
          && childMarkerNodes[0] == pathEndNodes[1]))
        {
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          /* Tested in TypingInnerRigidOnePathDoubleChild. */
          newcolumn->remainsGraphic = false;
        }
      }
    }
    else if (numPathEndNodes == 4)
    {
      if (reducedMember->rigidEndNodes[0] != parentMarkerNodes[0] && reducedMember->rigidEndNodes[0] != parentMarkerNodes[1])
      {
        TUdbgMsg(6 + 2*depth, "First path does not start at parent marker edge.\n");
        /* Tested in TypingInnerRigidTwoPathsNonadjacentParent. */
        newcolumn->remainsGraphic = false;
        return TU_OKAY;
      }
      if (reducedMember->rigidEndNodes[2] != parentMarkerNodes[0] && reducedMember->rigidEndNodes[2] != parentMarkerNodes[1])
      {
        TUdbgMsg(6 + 2*depth, "Second path does not start at parent marker edge.\n");
        /* Tested in TypingInnerRigidTwoPathsNonadjacentParent. */
        newcolumn->remainsGraphic = false;
        return TU_OKAY;
      }

      if (numOneEnd == 1)
      {
        bool pathConnects[2] = {
          reducedMember->rigidEndNodes[1] == childMarkerNodes[0]
            || reducedMember->rigidEndNodes[1] == childMarkerNodes[1],
          reducedMember->rigidEndNodes[3] == childMarkerNodes[0]
            || reducedMember->rigidEndNodes[3] == childMarkerNodes[1]
        };

        if (pathConnects[0] && pathConnects[1])
        {
          TUdbgMsg(6 + 2*depth, "The both paths end at the child marker edge.\n");
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else if (pathConnects[0])
        {
          TUdbgMsg(6 + 2*depth, "Only the first path ends at the child marker edge.\n");
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else if (pathConnects[1])
        {
          TUdbgMsg(6 + 2*depth, "Only the second path ends at the child marker edge.\n");
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          TUdbgMsg(6 + 2*depth, "No path ends at the child marker edge.\n");
          /* Tested in TypingInnerRigidTwoPathOneSingleChild. */
          newcolumn->remainsGraphic = false;
        }
      }
      else if (numOneEnd == 2)
      {
        bool pathConnected[2] = { false, false };
        bool childConnected[2] = { false, false };
        for (int i = 0; i < 2; ++i)
        {
          for (int j = 0; j < 4; ++j)
          {
            if (reducedMember->rigidEndNodes[1 + 2*i] == childMarkerNodes[j])
            {
              pathConnected[i] = true;
              childConnected[j/2] = true;
            }
          }
        }
        if (pathConnected[0] && pathConnected[1] && childConnected[0] && childConnected[1])
        {
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          TUdbgMsg(6 + 2*depth, "No pairing of paths to nodes of child marker edges possible.\n");
          /* Tested in TypingInnerRigidTwoPathsTwoSingleChildren. */
          newcolumn->remainsGraphic = false;
        }
      }
      else if (numTwoEnds == 0)
      {
        reducedMember->type = TYPE_DOUBLE_CHILD;
      }
      else
      {
        /* One child marker with both ends. We already checked that path end nodes 0 and 2 are parent marker nodes. */
        assert(numTwoEnds == 1);

        if ((pathEndNodes[1] == childMarkerNodes[0] && pathEndNodes[3] == childMarkerNodes[1])
          || (pathEndNodes[1] == childMarkerNodes[1] && pathEndNodes[3] == childMarkerNodes[0]))
        {
          reducedMember->type = TYPE_DOUBLE_CHILD;
        }
        else
        {
          TUdbgMsg(6 + 2*depth, "Paths do not end at child marker.");
          /* Tested in TypingInnerRigidTwoPathsDoubleChild. */
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
  DEC* dec,                      /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);

  TUdbgMsg(4 + 2*depth, "determineTypes(member %d = reduced member %ld)\n", reducedMember->member,
    reducedMember - &newcolumn->reducedMembers[0]);

  /* First handle children recursively. */
  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    TU_CALL( determineTypes(tu, dec, newcolumn, reducedComponent, reducedMember->children[c], depth + 1) );

    if (newcolumn->remainsGraphic)
    {
      TUdbgMsg(6 + 2*depth, "Child member %d of %d has type %d\n", reducedMember->children[c]->member, reducedMember->member,
        reducedMember->children[c]->type);
    }
    else
    {
      TUdbgMsg(6 + 2*depth, "Child prohibits graphicness.\n");
      return TU_OKAY;
    }
  }

  int numOneEnd, numTwoEnds;
  DEC_EDGE childMarkerEdges[2];
  TU_CALL( countChildrenTypes(tu, dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

#if defined(TU_DEBUG)
  int countReducedEdges = 0;
  for (PathEdge* e = reducedMember->firstPathEdge; e; e = e->next)
    ++countReducedEdges;
  TUdbgMsg(6 + 2*depth, "Member %d has %d children with one end and %d with two ends and %d path edges.\n",
    reducedMember->member, numOneEnd, numTwoEnds, countReducedEdges);
#endif /* TU_DEBUG */

  if (2*numTwoEnds + numOneEnd > 2)
  {
    /* Tested in Graphic.TypingManyChildrenTerminals */
    newcolumn->remainsGraphic = false;
    return TU_OKAY;
  }

  bool isRoot = reducedMember == reducedComponent->root;
  DEC_MEMBER member = findMember(dec, reducedMember->member);
  if (dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL)
  {
    TU_CALL( determineTypeParallel(tu, dec, newcolumn, reducedComponent, reducedMember, numOneEnd, numTwoEnds,
      childMarkerEdges, depth) );
  }
  else if (dec->members[member].type == DEC_MEMBER_TYPE_SERIES)
  {
    TU_CALL( determineTypeSeries(tu, dec, newcolumn, reducedComponent, reducedMember, numOneEnd, numTwoEnds,
      childMarkerEdges, depth) );
  }
  else
  {
    assert(dec->members[member].type == DEC_MEMBER_TYPE_RIGID);
    TU_CALL( determineTypeRigid(tu, dec, newcolumn, reducedComponent, reducedMember, numOneEnd, numTwoEnds,
      childMarkerEdges, depth) );
  }

  TUdbgMsg(6 + 2*depth, "Determined type %d (%s).\n", reducedMember->type, newcolumn->remainsGraphic ? "graphic" : "non-graphic");

  /* Parent marker edge closes cycle, so we propagate information to parent. */

  if (newcolumn->remainsGraphic && !isRoot && reducedMember->type == TYPE_CYCLE_CHILD)
  {
    DEC_MEMBER parentMember = findMemberParent(dec, reducedMember->member);
    ReducedMember* reducedParent = newcolumn->membersToReducedMembers[parentMember];
    DEC_EDGE markerOfParent = dec->members[member].markerOfParent;

    TUdbgMsg(6 + 2*depth, "Marker edge closes cycle.\n");
    TUdbgMsg(6 + 2*depth, "Parent member %d is reduced member %ld.\n", parentMember,
      (reducedParent - newcolumn->reducedMembers));

    /* Add marker edge of parent to reduced parent's path edges. */
    PathEdge* reducedEdge = NULL;
    TU_CALL( createPathEdge(tu, newcolumn, markerOfParent, reducedParent->firstPathEdge, &reducedEdge) );
    reducedParent->firstPathEdge = reducedEdge;

    /* Indicate that marker edge of parent belongs to path. */
    newcolumn->edgesInPath[markerOfParent] = true;

    /* Increase node degrees of nodes in a rigid parent. */
    if (dec->members[reducedParent->member].type == DEC_MEMBER_TYPE_RIGID)
    {
      newcolumn->nodesDegree[findEdgeHead(dec, markerOfParent)]++;
      newcolumn->nodesDegree[findEdgeTail(dec, markerOfParent)]++;
    }

    TUdbgMsg(6 + 2*depth, "Added marker edge of parent to list of path edges.\n");
  }

  return TU_OKAY;
}

TU_ERROR addColumnCheck(
  TU* tu,                       /**< TU environment. */
  DEC* dec,                    /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn, /**< new-column structure. */
  int* entryRows,               /**< Array of rows with 1-entry in this column. */
  int numEntries                /**< Number of 1-entries in this column. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);
  assert(entryRows);

  TUdbgMsg(0, "\n  Checking whether we can add a column with %d 1's.\n", numEntries);

  TUconsistencyAssert( decConsistency(dec) );

  /* Check for the (not yet computed) reduced decomposition whether there is a pair of parallel child/parent marker
   * edges in a non-parallel. */
  TU_CALL( parallelParentChildCheckReducedMembers(tu, dec, entryRows, numEntries) );

  TU_CALL( initializeNewColumn(tu, dec, newcolumn) );
  TU_CALL( computeReducedDecomposition(tu, dec, newcolumn, entryRows, numEntries) );
  TU_CALL( initializeReducedMemberEdgeLists(tu, dec, newcolumn, entryRows, numEntries) );

  debugDot(tu, dec, newcolumn);

  for (int i = 0; i < newcolumn->numReducedComponents; ++i)
  {
    TU_CALL( determineTypes(tu, dec, newcolumn, &newcolumn->reducedComponents[i],
      newcolumn->reducedComponents[i].root, 0) );

    TUdbgMsg(6, "After inspecting reduced component %d, graphic = %d.\n", i, newcolumn->remainsGraphic);
  }

  if (newcolumn->remainsGraphic)
    TUdbgMsg(4, "Adding the column would maintain graphicness.\n");

  TUconsistencyAssert( decConsistency(dec) );

  return TU_OKAY;
}

static
TU_ERROR addTerminal(
  TU* tu,                             /**< \ref TU environment. */
  DEC* dec,                      /**< t-decomposition. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  DEC_MEMBER member,              /**< New terminal member. */
  DEC_NODE node                   /**< New terminal node. */
)
{
  assert(reducedComponent);
  assert(member >= 0);
  assert(isRepresentativeMember(dec, member));

  /* For parallels we don't need to provide a node. */
  assert(node >= 0 || dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL);
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
  DEC* dec,                      /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent  /**< Reduced component. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);
  assert(reducedComponent);

  ReducedMember* reducedMember = reducedComponent->root;
  DEC_MEMBER member = reducedMember->member;
  assert(isRepresentativeMember(dec, member));

  int numOneEnd, numTwoEnds;
  DEC_EDGE childMarkerEdges[2];
  TU_CALL( countChildrenTypes(tu, dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  bool cycleWithUniqueEndChild;
  if (dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL)
  {
    cycleWithUniqueEndChild = (numTwoEnds == 1 || numOneEnd == 1);
  }
  else if (dec->members[member].type == DEC_MEMBER_TYPE_RIGID)
  {
    if (numTwoEnds == 1 || numOneEnd == 1)
    {
      /* For root rigid members, we have to check whether the path edges form a cycle with a two-end child. */
      DEC_NODE childMarkerNodes[2] = {
        findEdgeTail(dec, childMarkerEdges[0]),
        findEdgeHead(dec, childMarkerEdges[0])
      };
      cycleWithUniqueEndChild = (reducedMember->rigidEndNodes[2] == -1
        && ((reducedMember->rigidEndNodes[0] == childMarkerNodes[0]
        && reducedMember->rigidEndNodes[1] == childMarkerNodes[1])
        || (reducedMember->rigidEndNodes[0] == childMarkerNodes[1]
        && reducedMember->rigidEndNodes[1] == childMarkerNodes[0])));
    }
    else
      cycleWithUniqueEndChild = false;
  }
  else
  {
    /* For root series, the parent marker is not a path edge, so we cannot close a cycle with a child marker edge. */
    cycleWithUniqueEndChild = false;
  }

  if (!cycleWithUniqueEndChild)
  {
    TUdbgMsg(6, "Reduced root member does not close a cycle with a 1- or 2-end child marker edge.\n");
    return TU_OKAY;
  }

  while (cycleWithUniqueEndChild)
  {
    TUdbgMsg(6, "Reduced %s member %d closes a cycle with a 1- or 2-end child marker edge.\n",
      memberTypeString(dec->members[member].type), member);

    DEC_MEMBER childMember = findMember(dec, dec->edges[childMarkerEdges[0]].childMember);
    newcolumn->edgesInPath[ dec->members[childMember].markerToParent ] = true;

    /* Find the unique child member in order to process that. */
    for (int c = 0; c < reducedMember->numChildren; ++c)
    {
      Type type = reducedMember->children[c]->type;
      if (type == TYPE_SINGLE_CHILD || type == TYPE_DOUBLE_CHILD)
      {
        reducedMember = reducedMember->children[c];
        break;
      }
    }

    member = reducedMember->member;
    assert(isRepresentativeMember(dec, member));
    TU_CALL( countChildrenTypes(tu, dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

    if (dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL)
    {
      assert(!reducedMember->firstPathEdge);
      cycleWithUniqueEndChild = (numOneEnd == 1 || numTwoEnds == 1);
    }
    else if (dec->members[member].type == DEC_MEMBER_TYPE_RIGID)
    {
      if (numTwoEnds == 1 || numOneEnd == 1)
      {
        /* For non-root rigid members we have to check whether the path edges together with the parent marker edge form
         * a cycle with a two-end child. */
        DEC_NODE parentMarkerNodes[2] = {
          findEdgeTail(dec, dec->members[member].markerToParent),
          findEdgeHead(dec, dec->members[member].markerToParent)
        };
        DEC_NODE childMarkerNodes[2] = {
          findEdgeTail(dec, childMarkerEdges[0]),
          findEdgeHead(dec, childMarkerEdges[0])
        };

        int numEndNodes = reducedMember->rigidEndNodes[0] < 0 ? 0 : (reducedMember->rigidEndNodes[2] < 0 ? 2 : 4);
        if (numEndNodes == 0)
        {
          /* Without path edges the child marker would have to be parallel to the parent marker, which is detected
           * during typing. */
          cycleWithUniqueEndChild = false;
          break;
        }

        /* Determine the end nodes of the path including the parent marker edge. */
        DEC_NODE endNodes[2];
        if (numEndNodes == 4)
        {
          endNodes[0] = reducedMember->rigidEndNodes[1];
          endNodes[1] = reducedMember->rigidEndNodes[3];
        }
        else if (reducedMember->rigidEndNodes[0] == parentMarkerNodes[0])
        {
          endNodes[0] = parentMarkerNodes[1];
          endNodes[1] = reducedMember->rigidEndNodes[1];
        }
        else
        {
          assert(reducedMember->rigidEndNodes[0] == parentMarkerNodes[1]);
          endNodes[0] = parentMarkerNodes[0];
          endNodes[1] = reducedMember->rigidEndNodes[1];
        }

        cycleWithUniqueEndChild = (endNodes[0] == childMarkerNodes[0] && endNodes[1] == childMarkerNodes[1])
          || (endNodes[0] == childMarkerNodes[1] && endNodes[1] == childMarkerNodes[0]);
      }
      else
        cycleWithUniqueEndChild = false;
    }
    else
    {
      assert(dec->members[member].type == DEC_MEMBER_TYPE_SERIES);
      if (numOneEnd == 1 || numTwoEnds == 1)
      {
        int countPathEdges = 1;
        for (PathEdge* pathEdge = reducedMember->firstPathEdge; pathEdge != NULL; pathEdge = pathEdge->next)
          ++countPathEdges;
        cycleWithUniqueEndChild = countPathEdges == dec->members[member].numEdges - 1;
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
  DEC* dec,      /**< t-decomposition. */
  DEC_EDGE edge,  /**< Reduced component. */
  DEC_NODE tail,  /**< New tail node. */
  DEC_NODE head   /**< New head node. */
)
{
  assert(tu);
  assert(dec);
  assert(edge >= 0);
  assert(edge < dec->memEdges);

  dec->edges[edge].tail = tail;
  dec->edges[edge].head = head;

  return TU_OKAY;
}

static
TU_ERROR mergeMemberIntoParent(
  TU* tu,                 /**< \ref TU environment. */
  DEC* dec,          /**< t-decomposition. */
  DEC_MEMBER member,  /**< Reduced component. */
  bool headToHead         /**< Whether the heads of the edges shall be joined. */
)
{
  assert(tu);
  assert(dec);
  assert(member >= 0);

  member = findMember(dec, member);
  DEC_MEMBER parentMember = findMemberParent(dec, member);
  assert(parentMember >= 0);

  TUdbgMsg(10, "Merging child member %d into its parent member %d.\n", member, parentMember);

#if defined(TU_DEBUG)
  DEC_EDGE edge = dec->members[member].firstEdge;
  do
  {
    if (dec->edges[edge].head < 0 || dec->edges[edge].tail < 0)
      TUdbgMsg(10, "Edge %d of merge member %d does not have nodes.\n", edge, member);
    assert(dec->edges[edge].tail >= 0);
    assert(dec->edges[edge].head >= 0);
    edge = dec->edges[edge].next;
  }
  while (edge != dec->members[member].firstEdge);
#endif /* TU_DEBUG */

  DEC_EDGE parentEdge = dec->members[member].markerOfParent;
  assert(parentEdge >= 0);
  DEC_EDGE childEdge = dec->members[member].markerToParent;
  assert(childEdge >= 0);

  DEC_NODE parentEdgeNodes[2] = { findEdgeTail(dec, parentEdge), findEdgeHead(dec, parentEdge) };
  DEC_NODE childEdgeNodes[2] = { findEdgeTail(dec, childEdge), findEdgeHead(dec, childEdge) };
  TUdbgMsg(10, "Merging member %d into %d, identifying %d = {%d,%d} with %d = {%d,%d}.\n", member, parentMember,
    childEdge, childEdgeNodes[0], childEdgeNodes[1], parentEdge, parentEdgeNodes[0], parentEdgeNodes[1]);

  /* Identify nodes. */

  dec->nodes[childEdgeNodes[0]].representativeNode = parentEdgeNodes[headToHead ? 0 : 1];
  dec->nodes[childEdgeNodes[1]].representativeNode = parentEdgeNodes[headToHead ? 1 : 0];

  /* Identify members. */

  dec->members[member].representativeMember = parentMember;

  /* We add the member's edges to the parent's edge list and thereby remove the two marker edges. */
  if (dec->members[parentMember].firstEdge == parentEdge)
    dec->members[parentMember].firstEdge = dec->edges[parentEdge].next;

  dec->edges[dec->edges[parentEdge].next].prev = dec->edges[childEdge].prev;
  dec->edges[dec->edges[parentEdge].prev].next = dec->edges[childEdge].next;
  dec->edges[dec->edges[childEdge].next].prev = dec->edges[parentEdge].prev;
  dec->edges[dec->edges[childEdge].prev].next = dec->edges[parentEdge].next;
  dec->members[parentMember].numEdges += dec->members[member].numEdges - 2;
  dec->numEdges -= 2;
  dec->edges[parentEdge].next = dec->firstFreeEdge;
  dec->edges[childEdge].next = parentEdge;
  dec->firstFreeEdge = childEdge;
  dec->members[parentMember].type = DEC_MEMBER_TYPE_RIGID;

  return TU_OKAY;
}

static
TU_ERROR createParallelNodes(
  TU* tu,           /**< \ref TU environment. */
  DEC* dec,        /**< t-decomposition. */
  DEC_MEMBER member /**< A parallel member. */
)
{
  assert(tu);
  assert(dec);
  assert(member >= 0);
  member = findMember(dec, member);
  assert(dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL);

  DEC_EDGE edge = dec->members[member].firstEdge;
  assert(dec->edges[edge].head < 0);
  assert(dec->edges[edge].tail < 0);
  if (dec->edges[edge].head >= 0)
  {
    assert(dec->edges[edge].tail >= 0);
    return TU_OKAY;
  }

  DEC_NODE tail, head;
  TU_CALL( createNode(dec, &tail) );
  TU_CALL( createNode(dec, &head) );

  do
  {
    assert(dec->edges[edge].tail < 0);
    assert(dec->edges[edge].head < 0);

    dec->edges[edge].tail = tail;
    dec->edges[edge].head = head;
    edge = dec->edges[edge].next;
  }
  while (edge != dec->members[member].firstEdge);

  return TU_OKAY;
}

static
TU_ERROR splitParallel(
  TU* tu,                     /**< \ref TU environment. */
  DEC* dec,                  /**< t-decomposition. */
  DEC_MEMBER parallel,        /**< A parallel member. */
  DEC_EDGE edge1,             /**< First edge to be moved into the child parallel. */
  DEC_EDGE edge2,             /**< Second edge to be moved into the child parallel. */
  DEC_MEMBER* pChildParallel  /**< Pointer to storing the newly created child parallel. */
)
{
  assert(tu);
  assert(dec);
  assert(parallel >= 0);
  assert(parallel < dec->memMembers);
  assert(edge1 >= 0);
  assert(edge1 < dec->memEdges);
  assert(edge2 >= 0);
  assert(edge2 < dec->memEdges);

  DEC_MEMBER childParallel;
  TU_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &childParallel) );
  DEC_EDGE markerOfParentParallel, markerOfChildParallel;
  TU_CALL( createMarkerEdgePair(dec, parallel, &markerOfParentParallel, -1, -1, childParallel, &markerOfChildParallel, -1, -1) );
  TU_CALL( addEdgeToMembersEdgeList(dec, markerOfParentParallel) );
  TU_CALL( addEdgeToMembersEdgeList(dec, markerOfChildParallel) );

  TU_CALL( removeEdgeFromMembersEdgeList(dec, edge1) );
  dec->edges[edge1].member = childParallel;
  TU_CALL( addEdgeToMembersEdgeList(dec, edge1) );
  dec->members[findMember(dec, dec->edges[edge1].childMember)].parentMember = childParallel;

  TU_CALL( removeEdgeFromMembersEdgeList(dec, edge2) );
  dec->edges[edge2].member = childParallel;
  TU_CALL( addEdgeToMembersEdgeList(dec, edge2) );
  dec->members[findMember(dec, dec->edges[edge2].childMember)].parentMember = childParallel;

  if (pChildParallel)
    *pChildParallel = childParallel;

  return TU_OKAY;
}

static
TU_ERROR addColumnProcessParallel(
  TU* tu,                             /**< \ref TU environment. */
  DEC* dec,                          /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn,           /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of this member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);
  assert(reducedMember);
  DEC_MEMBER member = findMember(dec, reducedMember->member);
  assert(dec->members[member].type == DEC_MEMBER_TYPE_PARALLEL);

  int numOneEnd, numTwoEnds;
  DEC_EDGE childMarkerEdges[2];
  TU_CALL( countChildrenTypes(tu, dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  TUdbgMsg(6 + 2*depth, "addColumnProcessParallel for%s member %d (reduced %ld), #one-ends = %d, #two-ends = %d of type %d.\n",
    depth == 0 ? " root" : "", member, (reducedMember - newcolumn->reducedMembers), numOneEnd, numTwoEnds, reducedMember->type);

  if (depth == 0)
  {
    if (numOneEnd == 0 && numTwoEnds == 0)
    {
      /* Tested in UpdateRootParallelNoChildren. */

      assert(reducedMember->firstPathEdge);
      TU_CALL( addTerminal(tu, dec, reducedComponent, member, -1) );
      TU_CALL( addTerminal(tu, dec, reducedComponent, member, -1) );
      return TU_OKAY;
    }
    else
    {
      /* Tested in UpdateRootParallelTwoSingleChildren. */

      /* If we have one single-child or a double-child then we should have moved the reduced root. */
      assert(numOneEnd == 2);
      assert(reducedComponent->numTerminals == 2);

      /* If the paralle contains more than 3 edges, we split the two child edges off, which creates a parallel with a
       * parent marker and the two children.
       * Test: Graphic.RootParallelTwoOneEnds */

      if (dec->members[member].numEdges > 3)
      {
        /* Tested in UpdateRootParallelNoChildrenSplit. */
        TU_CALL( splitParallel(tu, dec, member, childMarkerEdges[0], childMarkerEdges[1], &member) );
        reducedMember->member = member;
      }

      debugDot(tu, dec, newcolumn);

      assert(dec->members[member].numEdges == 3);
      TU_CALL( createParallelNodes(tu, dec, member) );
      TU_CALL( mergeMemberIntoParent(tu, dec, dec->edges[childMarkerEdges[0]].childMember, true) );
      TU_CALL( mergeMemberIntoParent(tu, dec, dec->edges[childMarkerEdges[1]].childMember,
        reducedMember->firstPathEdge == NULL && reducedMember->type != TYPE_DOUBLE_CHILD) );

      debugDot(tu, dec, newcolumn);

      return TU_OKAY;
    }
  }
  else
  {
    /* Tested in UpdateInnerParallelOneSingleChild. */

    /* Cannot be a leaf since then it would contain a path edge, i.e., it closes a cycle. */
    assert(numOneEnd == 1);
    assert(reducedComponent->numTerminals >= 1);

    DEC_NODE tail, head;
    TU_CALL( createNode(dec, &tail) );
    TU_CALL( createNode(dec, &head) );
    TUdbgMsg(8 + 2*depth, "Parallel's tail node is %d and head node is %d.\n", tail, head);
    DEC_EDGE edge = dec->members[member].firstEdge;
    do
    {
      TU_CALL( setEdgeNodes(tu, dec, edge, tail, head) );
      edge = dec->edges[edge].next;
    }
    while (edge != dec->members[member].firstEdge);

    TU_CALL( mergeMemberIntoParent(tu, dec, dec->edges[childMarkerEdges[0]].childMember,
      !reducedMember->firstPathEdge) );

    debugDot(tu, dec, newcolumn);
  }

  return TU_OKAY;
}

inline static
void flipEdge(
  DEC* dec,    /**< t-decomposition. */
  DEC_EDGE edge /**< edge. */
)
{
  assert(dec);
  assert(edge >= 0);
  assert(edge < dec->memEdges);

  SWAP_INTS(dec->edges[edge].tail, dec->edges[edge].head);
}

static
TU_ERROR addColumnProcessRigid(
  TU* tu,                             /**< \ref TU environment. */
  DEC* dec,                      /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of this member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);
  assert(reducedMember);
  DEC_MEMBER member = findMember(dec, reducedMember->member);
  assert(dec->members[member].type == DEC_MEMBER_TYPE_RIGID);

  int numOneEnd, numTwoEnds;
  DEC_EDGE childMarkerEdges[2];
  TU_CALL( countChildrenTypes(tu, dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  TUdbgMsg(6 + 2*depth, "addColumnProcessRigid for%s member %d (reduced %ld), #one-ends = %d, #two-ends = %d of type %d.\n",
    depth == 0 ? " root" : "", reducedMember->member, (reducedMember - newcolumn->reducedMembers), numOneEnd,
    numTwoEnds, reducedMember->type);

  DEC_NODE parentMarkerNodes[2] = {
    dec->members[member].markerToParent < 0 ? -1 : findEdgeTail(dec, dec->members[member].markerToParent),
    dec->members[member].markerToParent < 0 ? -1 : findEdgeHead(dec, dec->members[member].markerToParent)
  };
  DEC_NODE childMarkerNodes[4] = {
    childMarkerEdges[0] < 0 ? -1 : findEdgeTail(dec, childMarkerEdges[0]),
    childMarkerEdges[0] < 0 ? -1 : findEdgeHead(dec, childMarkerEdges[0]),
    childMarkerEdges[1] < 0 ? -1 : findEdgeTail(dec, childMarkerEdges[1]),
    childMarkerEdges[1] < 0 ? -1 : findEdgeHead(dec, childMarkerEdges[1])
  };

  DEC_NODE* pathEndNodes = reducedMember->rigidEndNodes;
  int numPathEndNodes = pathEndNodes[0] < 0 ? 0 : (pathEndNodes[2] < 0 ? 2 : 4 );

  if (depth == 0)
  {
    /* Root rigid member. */

    if (dec->members[member].markerToParent >= 0 && newcolumn->edgesInPath[dec->members[member].markerToParent])
    {
      /* The parent marker is a path edge, so we modify the endNodes array. */
      newcolumn->nodesDegree[parentMarkerNodes[0]]++;
      newcolumn->nodesDegree[parentMarkerNodes[1]]++;
      if (numPathEndNodes == 0)
      {
        /* Tested in UpdateRootRigidParentOnly. */

        pathEndNodes[0] = parentMarkerNodes[0];
        pathEndNodes[1] = parentMarkerNodes[1];
        numPathEndNodes = 1;
      }
      else if (numPathEndNodes == 2)
      {
        if (pathEndNodes[0] == parentMarkerNodes[0])
          pathEndNodes[0] = parentMarkerNodes[1];
        else if (pathEndNodes[0] == parentMarkerNodes[1])
          pathEndNodes[0] = parentMarkerNodes[0];
      }
      else
      {
        /* Tested in UpdateRootRigidParentJoinsPaths. */

        pathEndNodes[0] = pathEndNodes[3];
        pathEndNodes[2] = -1;
        pathEndNodes[3] = -1;
        numPathEndNodes = 2;
      }
    }

    assert(numPathEndNodes <= 2);

    if (numOneEnd == 0 && numTwoEnds == 0)
    {
      /* Tested in UpdateRootRigidNoChildren. */

      assert(numPathEndNodes >= 2);
      TU_CALL( addTerminal(tu, dec, reducedComponent, member, pathEndNodes[0]) );
      TU_CALL( addTerminal(tu, dec, reducedComponent, member, pathEndNodes[1]) );
    }
    else if (numOneEnd == 1)
    {
      /* Tested in UpdateRootRigidOneSingleChild. */

      TU_CALL( addTerminal(tu, dec, reducedComponent, member,
        (pathEndNodes[0] == childMarkerNodes[0] || pathEndNodes[0] == childMarkerNodes[1]) ? pathEndNodes[1] : pathEndNodes[0] ) );

      DEC_MEMBER childMember = findMember(dec, dec->edges[childMarkerEdges[0]].childMember);
      bool headToHead = pathEndNodes[0] == childMarkerNodes[1] || pathEndNodes[1] == childMarkerNodes[1];
      assert(headToHead || pathEndNodes[0] == childMarkerNodes[0] || pathEndNodes[1] == childMarkerNodes[0]);

      TU_CALL( mergeMemberIntoParent(tu, dec, childMember, headToHead) );
    }
    else
    {
      /* Tested in UpdateRootRigidTwoSingleChildren. */

      assert(numOneEnd == 2);
      assert(reducedComponent->numTerminals == 2);

      DEC_MEMBER childMember[2] = {
        findMember(dec, dec->edges[childMarkerEdges[0]].childMember),
        findMember(dec, dec->edges[childMarkerEdges[1]].childMember)
      };

      /* Count to how many path end nodes each child marker is incident. */
      int numIncidentPathNodes[2] = { 0, 0 };
      for (int c = 0; c < 2; ++c)
      {
        for (int i = 0; i < numPathEndNodes; ++i)
        {
          for (int j = 0; j < 2; ++j)
          {
            if (pathEndNodes[i] == childMarkerNodes[2*c + j])
              numIncidentPathNodes[c]++;
          }
        }
      }
      TUdbgMsg(8 + 2*depth,
        "Child marker %d = {%d,%d} has %d incident path end nodes and child marker %d = {%d,%d} has %d.\n",
        childMarkerEdges[0], childMarkerNodes[0], childMarkerNodes[1], numIncidentPathNodes[0], childMarkerEdges[1],
        childMarkerNodes[2], childMarkerNodes[3], numIncidentPathNodes[1]);

      /* If a child marker is incident to both path ends, then we ensure it is the second one. */
      if (numIncidentPathNodes[0] == 2)
      {
        SWAP_INTS(childMember[0], childMember[1]);
        SWAP_INTS(childMarkerNodes[0], childMarkerNodes[2]);
        SWAP_INTS(childMarkerNodes[1], childMarkerNodes[3]);
      }

      /* Check if the two child markers are parallel. We then create a parallel with the two. */
      if ((childMarkerNodes[0] == childMarkerNodes[2] && childMarkerNodes[1] == childMarkerNodes[3])
        || (childMarkerNodes[0] == childMarkerNodes[3] && childMarkerNodes[1] == childMarkerNodes[2]))
      {
        /* Tested in UpdateRootRigidTwoParallelSingleChildren. */

        TUdbgMsg(8, "Moving child marker edges %d = {%d,%d} and %d = {%d,%d} to a new parallel.\n", childMarkerEdges[0],
          childMarkerNodes[0], childMarkerNodes[1], childMarkerEdges[1], childMarkerNodes[2], childMarkerNodes[3]);

        DEC_MEMBER newParallel = -1;
        TU_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &newParallel) );
        dec->members[newParallel].parentMember = member;
        dec->members[childMember[0]].parentMember = newParallel;
        dec->members[childMember[1]].parentMember = newParallel;

        DEC_EDGE markerOfParent, markerToParent;
        TU_CALL( createMarkerEdgePair(dec, member, &markerOfParent, childMarkerNodes[0], childMarkerNodes[1],
          newParallel, &markerToParent, -1, -1) );
        
        TU_CALL( replaceEdgeInMembersEdgeList(dec, childMarkerEdges[0], markerOfParent) );
        TU_CALL( addEdgeToMembersEdgeList(dec, markerToParent) );
        dec->edges[childMarkerEdges[0]].member = newParallel;
        TU_CALL( addEdgeToMembersEdgeList(dec, childMarkerEdges[0]) );
        TU_CALL( removeEdgeFromMembersEdgeList(dec, childMarkerEdges[1]) );
        dec->edges[childMarkerEdges[1]].member = newParallel;
        TU_CALL( addEdgeToMembersEdgeList(dec, childMarkerEdges[1]) );

        /* The parallel has no nodes, so we have to get rid of these. */
        dec->edges[childMarkerEdges[0]].tail = -1;
        dec->edges[childMarkerEdges[0]].head = -1;
        dec->edges[childMarkerEdges[1]].tail = -1;
        dec->edges[childMarkerEdges[1]].head = -1;

        debugDot(tu, dec, newcolumn);

        /* We have to merge in the parallel. */
        TU_CALL( createParallelNodes(tu, dec, newParallel) );
        TU_CALL( mergeMemberIntoParent(tu, dec, dec->edges[childMarkerEdges[0]].childMember, true) );
        /* If there is no path, then we merge heads to heads. Otherwise, one head is mapped to tail, such that the path
         * is used. */
        TU_CALL( mergeMemberIntoParent(tu, dec, dec->edges[childMarkerEdges[1]].childMember, numPathEndNodes == 0) );

        debugDot(tu, dec, newcolumn);

        return TU_OKAY;
      }

      if (numPathEndNodes == 0)
      {
        /* We fake existence of a path by setting both path end nodes to the node common to both child markers. */
        if (childMarkerNodes[0] == childMarkerNodes[2] || childMarkerNodes[0] == childMarkerNodes[3])
        {
          pathEndNodes[0] = childMarkerNodes[0];
          pathEndNodes[1] = childMarkerNodes[0];
        }
        else
        {
          assert(childMarkerNodes[1] == childMarkerNodes[2] || childMarkerNodes[1] == childMarkerNodes[3]);
          pathEndNodes[0] = childMarkerNodes[1];
          pathEndNodes[1] = childMarkerNodes[1];
        }
      }

      if (pathEndNodes[0] != childMarkerNodes[0] && pathEndNodes[0] != childMarkerNodes[1])
        SWAP_INTS(pathEndNodes[0], pathEndNodes[1]);

      TUdbgMsg(8 + 2*depth, "After swapping, we have a %d-%d-path as well as child markers {%d,%d} and {%d,%d}.\n",
        pathEndNodes[0], pathEndNodes[1], childMarkerNodes[0], childMarkerNodes[1], childMarkerNodes[2], childMarkerNodes[3]);

      assert(pathEndNodes[0] == childMarkerNodes[0] || pathEndNodes[0] == childMarkerNodes[1]);
      TU_CALL( mergeMemberIntoParent(tu, dec, childMember[0], pathEndNodes[0] == childMarkerNodes[1]) );
      debugDot(tu, dec, newcolumn);

      TU_CALL( mergeMemberIntoParent(tu, dec, childMember[1], pathEndNodes[1] == childMarkerNodes[3]) );
      debugDot(tu, dec, newcolumn);
    }
  }
  else
  {
    /* Non-root rigid member. */

    if (numOneEnd == 0 && numTwoEnds == 0)
    {
      /* Tested in UpdateLeafRigid. */

      assert(reducedMember->firstPathEdge);
      assert(pathEndNodes[0] >= 0);

      TU_CALL( addTerminal(tu, dec, reducedComponent, member, pathEndNodes[1]) );
      if (parentMarkerNodes[0] == pathEndNodes[0])
        flipEdge(dec, dec->members[member].markerToParent);
    }
    else
    {
      assert(numOneEnd == 1);

      if (numPathEndNodes >= 2)
      {
        /* Tested in UpdateInnerRigidOnePath. */

        TUdbgMsg(8 + 2*depth, "%d-%d-path with parent marker {%d,%d} and child marker {%d,%d}.\n", pathEndNodes[0],
          pathEndNodes[1], parentMarkerNodes[0], parentMarkerNodes[1], childMarkerNodes[0], childMarkerNodes[1]);

        /* Ensure that the child marker is incident to the path end node 1. */
        if (pathEndNodes[1] != childMarkerNodes[0] && pathEndNodes[1] != childMarkerNodes[1])
          SWAP_INTS(pathEndNodes[0], pathEndNodes[1]);
        assert(pathEndNodes[1] == childMarkerNodes[0] || pathEndNodes[1] == childMarkerNodes[1]);

        /* Flip parent if necessary. */
        if (pathEndNodes[0] == parentMarkerNodes[0])
        {
          flipEdge(dec, dec->members[member].markerToParent);
          SWAP_INTS(parentMarkerNodes[0], parentMarkerNodes[1]);
        }

        assert(pathEndNodes[0] == parentMarkerNodes[1]);

        /* Merge child. */
        TU_CALL( mergeMemberIntoParent(tu, dec, findMember(dec, dec->edges[childMarkerEdges[0]].childMember),
          pathEndNodes[1] == childMarkerNodes[1]) );
        debugDot(tu, dec, newcolumn);
      }
      else
      {
        /* Tested in UpdateInnerRigidNoPath. */

        /* Parent marker and child marker must be next to each other. */
        if (parentMarkerNodes[0] == childMarkerNodes[0] || parentMarkerNodes[0] == childMarkerNodes[1])
          flipEdge(dec, dec->members[member].markerToParent);

        TU_CALL( mergeMemberIntoParent(tu, dec, dec->edges[childMarkerEdges[0]].childMember,
          parentMarkerNodes[0] == childMarkerNodes[1] || parentMarkerNodes[1] == childMarkerNodes[1]) );

        debugDot(tu, dec, newcolumn);
      }
    }
  }

  return TU_OKAY;
}

/**
 * \brief Replaces an edge by a parallel containing it.
 *
 * The given member should have at least two edges.
 */

static
TU_ERROR createEdgeParallel(
  TU* tu,                   /**< \ref TU environment. */
  DEC* dec,                /**< t-decomposition. */
  DEC_EDGE edge,            /**< Edge to be replaced by a parallel containing that edge. */
  DEC_MEMBER* pNewParallel  /**< Pointer for storing the new parallel. */
)
{
  assert(tu);
  assert(dec);

  TUdbgMsg(8, "Creating parallel for edge %d.\n", edge);

  DEC_MEMBER parentMember = findEdgeMember(dec, edge);
  DEC_MEMBER newParallel = -1;
  TU_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &newParallel) );
  dec->members[newParallel].parentMember = parentMember;

  DEC_EDGE markerOfParent, markerToParent;
  TU_CALL( createMarkerEdgePair(dec, parentMember, &markerOfParent, dec->edges[edge].tail, dec->edges[edge].head,
    newParallel, &markerToParent, -1, -1) );
  dec->edges[markerOfParent].next = dec->edges[edge].next;
  dec->edges[markerOfParent].prev = dec->edges[edge].prev;
  assert(dec->edges[markerOfParent].next != markerOfParent);
  dec->edges[dec->edges[markerOfParent].next].prev = markerOfParent;
  dec->edges[dec->edges[markerOfParent].prev].next = markerOfParent;
  if (dec->members[parentMember].firstEdge == edge)
    dec->members[parentMember].firstEdge = markerOfParent;

  TU_CALL( addEdgeToMembersEdgeList(dec, markerToParent) );

  dec->edges[edge].member = newParallel;
  TU_CALL( addEdgeToMembersEdgeList(dec, edge) );

  if (pNewParallel)
    *pNewParallel = newParallel;

  return TU_OKAY;
}

/**
 * \brief Splits series member into two, connecting them via a parallel.
 *
 * Takes all edges of the series \p member for which \p edgesPredicate is the same as
 * \p predicateValue.
 */

static
TU_ERROR splitSeries(
  TU* tu,                             /**< \ref TU environment. */
  DEC* dec,                      /**< t-decomposition. */
  DEC_MEMBER member,              /**< Series member to be squeezed. */
  bool* edgesPredicate,               /**< Map from edges to predicate. */
  bool predicateValue,                /**< Value of predicate. */
  DEC_EDGE* pRepresentativeEdge,  /**< Pointer for storing the child marker edge that links to new parallel. */
  DEC_MEMBER* pNewParallel,           /**< Pointer for storing the new parallel. */
  DEC_MEMBER* pNewSeries         /**< Pointer for storing the new parallel. */
)
{
  assert(tu);
  assert(dec);
  assert(member >= 0);
  assert(member < dec->memMembers);
  assert(edgesPredicate);
  assert(dec->members[member].type == DEC_MEMBER_TYPE_SERIES);

#if defined(TU_DEBUG)
  TUdbgMsg(8, "Checking series member %d for splitting... ", member);
#endif /* TU_DEBUG */

  DEC_EDGE edge = dec->members[member].firstEdge;
  DEC_EDGE someSatisfyingEdge = -1;
  int numSatisfyingEdges = 0;
  do
  {
    bool value = edgesPredicate[edge];
    if ((value && predicateValue) || (!value && !predicateValue))
    {
      someSatisfyingEdge = edge;
      ++numSatisfyingEdges;
    }
    edge = dec->edges[edge].next;
  } while (edge != dec->members[member].firstEdge);

  TUdbgMsg(0, "%d of %d edges satisfy the predicate.\n", numSatisfyingEdges, dec->members[member].numEdges);

  if (numSatisfyingEdges == 0)
  {
    if (pRepresentativeEdge)
      *pRepresentativeEdge = -1;
    if (pNewParallel)
      *pNewParallel = -1;
    if (pNewSeries)
      *pNewSeries = -1;
  }
  else if (numSatisfyingEdges == 1)
  {
    if (pNewSeries)
      *pNewSeries = -1;
    if (pRepresentativeEdge)
      *pRepresentativeEdge = someSatisfyingEdge;
    if (pNewParallel)
      *pNewParallel = -1;
  }
  else
  {
    /* Initialize new series member. */
    DEC_MEMBER series;
    TU_CALL( createMember(dec, DEC_MEMBER_TYPE_SERIES, &series) );

    /* Initialize new parallel. */
    DEC_MEMBER parallel;
    TU_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &parallel) );

    DEC_EDGE seriesParentMarker, parallelChildMarker;
    TU_CALL( createMarkerEdgePair(dec, parallel, &parallelChildMarker, -1, -1, series, &seriesParentMarker, -1, -1) );
    TU_CALL( addEdgeToMembersEdgeList(dec, seriesParentMarker) );
    TU_CALL( addEdgeToMembersEdgeList(dec, parallelChildMarker) );

    /* Go through old series member. */
    DEC_EDGE firstEdge = dec->members[member].firstEdge;
    edge = firstEdge;
    bool encounteredStayingEdge = false;
    do
    {
#if defined(TU_DEBUG_SPLITTING)
      TUdbgMsg(8, "Edge %d <%d>", edge, dec->edges[edge].name);
      if (dec->edges[edge].childMember >= 0)
        TUdbgMsg(0, " (with child %d)", dec->edges[edge].childMember);
      if (edge == dec->members[member].markerToParent)
        TUdbgMsg(0, " (with parent %d)", dec->members[member].parentMember);
      TUdbgMsg(0, " (prev = %d, next = %d)", dec->edges[edge].prev, dec->edges[edge].next);
#endif /* TU_DEBUG_SPLITTING*/

      /* Evaluate predicate. */
      bool value = edgesPredicate[edge];
      if ((value && !predicateValue) || (!value && predicateValue))
      {
#if defined(TU_DEBUG_SPLITTING)
        TUdbgMsg(" does not satisfy the predicate.\n");
#endif /* TU_DEBUG_SPLITTING */
        edge = dec->edges[edge].next;
        encounteredStayingEdge = true;
        continue;
      }

#if defined(TU_DEBUG_SPLITTING)
      TUdbgMsg(0, " satisfies the predicate.\n");
#endif /* TU_DEBUG_SPLITTING */

      assert(edge != dec->members[member].markerToParent);

      /* Remove edge from old edge list. */
      DEC_EDGE oldPrev = dec->edges[edge].prev;
      DEC_EDGE oldNext = dec->edges[edge].next;
      dec->edges[oldPrev].next = oldNext;
      dec->edges[oldNext].prev = oldPrev;
      dec->members[member].numEdges--;

      /* Add edge to new edge list. */
      DEC_EDGE newPrev = dec->edges[seriesParentMarker].prev;
      dec->edges[newPrev].next = edge;
      dec->edges[seriesParentMarker].prev = edge;
      dec->edges[edge].prev = newPrev;
      dec->edges[edge].next = seriesParentMarker;
      dec->edges[edge].member = series;
      if (dec->edges[edge].childMember >= 0)
      {
        assert( dec->members[dec->edges[edge].childMember].parentMember == member);
        dec->members[dec->edges[edge].childMember].parentMember = series;
      }
      dec->members[series].numEdges++;

      /* Did we move the first edge of this member? */
      if (edge == firstEdge)
      {
        dec->members[member].firstEdge = oldNext;
        firstEdge = oldNext;
        edge = oldNext;
        continue;
      }

      edge = oldNext;
    }
    while (edge != firstEdge || !encounteredStayingEdge);

    DEC_EDGE memberChildMarker, parallelParentMarker;
    TU_CALL( createMarkerEdgePair(dec, member, &memberChildMarker, -1, -1, parallel, &parallelParentMarker, -1, -1) );
    TU_CALL( addEdgeToMembersEdgeList(dec, parallelParentMarker) );
    DEC_EDGE oldPrev = dec->edges[firstEdge].prev;
    dec->edges[memberChildMarker].next = firstEdge;
    dec->edges[memberChildMarker].prev = oldPrev;
    dec->edges[oldPrev].next = memberChildMarker;
    dec->edges[firstEdge].prev = memberChildMarker;
    dec->members[member].numEdges++;

#if defined(TU_DEBUG_SPLITTING)
    TUdbgMsg(8, "Updated old series member with these edges:");
    edge = firstEdge;
    do
    {
      TUdbgMsg(0, " %d <%d>", edge, dec->edges[edge].name);
      edge = dec->edges[edge].next;
    }
    while (edge != firstEdge);
    TUdbgMsg(0, ".\n");
    TUdbgMsg(8, "New series member has these edges:");
    edge = seriesParentMarker;
    do
    {
      TUdbgMsg(0, " %d <%d>", edge, dec->edges[edge].name);
      edge = dec->edges[edge].next;
    }
    while (edge != seriesParentMarker);
    TUdbgMsg(0, ".\n");
#endif /* TU_DEBUG_SPLITTING */

    TUdbgMsg(8, "Connecting parallel is member %d and squeezed series member is member %d.\n", parallel, series);

    if (pRepresentativeEdge)
      *pRepresentativeEdge = memberChildMarker;
    if (pNewParallel)
      *pNewParallel = parallel;
    if (pNewSeries)
      *pNewSeries = series;
  }
  return TU_OKAY;
}

static
TU_ERROR addColumnProcessSeries(
  TU* tu,                             /**< \ref TU environment. */
  DEC* dec,                      /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of this member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);
  assert(reducedMember);

  DEC_MEMBER member = findMember(dec, reducedMember->member);

  int numOneEnd, numTwoEnds;
  DEC_EDGE childMarkerEdges[2];
  TU_CALL( countChildrenTypes(tu, dec, reducedMember, &numOneEnd, &numTwoEnds, childMarkerEdges) );

  TUdbgMsg(6 + 2*depth, "addColumnProcessSeries for%s member %d (reduced %ld), #one-ends = %d, #two-ends = %d of type %d.\n",
    depth == 0 ? " root" : "", member, (reducedMember - newcolumn->reducedMembers), numOneEnd, numTwoEnds, reducedMember->type);

  if (depth == 0)
  {
    /* If we have a child containing both ends then we should have moved the reduced root there. */
    assert(numTwoEnds == 0);

    if (numOneEnd == 0)
    {
      /* Root series member containing both ends. */
      assert(reducedMember->firstPathEdge);
      DEC_EDGE representativeEdge;
      if (newcolumn->edgesInPath[dec->members[member].markerToParent])
      {
        /* Tested in UpdateRootSeriesNoChildrenParent. */

        TUdbgMsg(8 + 2*depth, "Series contains both terminal nodes and the parent marker edge is a path edge.\n");

        /* Squeeze off all non-path edges by moving them to a new series member and creating a parallel member to
         * connect it to the remaining series member. */
        TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, false, &representativeEdge, NULL, NULL) );
      }
      else if (reducedMember->type == TYPE_CYCLE_CHILD)
      {
        /* Tested in UpdateRootSeriesNoChildrenHamiltonianPath. */

        TUdbgMsg(8 + 2*depth, "Series contains both terminal nodes which are the parent marker edge nodes.\n");

        DEC_MEMBER parentMember = findMemberParent(dec, member);
        DEC_EDGE markerOfParent = dec->members[member].markerOfParent;
        assert(parentMember >= 0); /* A series member can only close a cycle if it has a parent. */
        if (dec->members[parentMember].type == DEC_MEMBER_TYPE_PARALLEL)
        {
          TU_CALL( addTerminal(tu, dec, reducedComponent, findEdgeMember(dec, markerOfParent), -1 ) );
          TU_CALL( addTerminal(tu, dec, reducedComponent, findEdgeMember(dec, markerOfParent), -1 ) );
        }
        else
        {
          assert(dec->members[parentMember].type == DEC_MEMBER_TYPE_RIGID);
          TU_CALL( addTerminal(tu, dec, reducedComponent, findEdgeMember(dec, markerOfParent), findEdgeTail(dec, markerOfParent) ) );
          TU_CALL( addTerminal(tu, dec, reducedComponent, findEdgeMember(dec, markerOfParent), findEdgeHead(dec, markerOfParent) ) );
        }

        return TU_OKAY;
      }
      else
      {
        /* Tested in UpdateRootSeriesNoChildren. */
        TUdbgMsg(8 + 2*depth, "Series contains both terminal nodes and a non-path parent marker edge.\n");

        /* Squeeze off all path edges by moving them to a new series member and creating a parallel member to connect it
         * to the remaining series member. */
        TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, true, &representativeEdge, NULL, NULL) );
      }

      debugDot(tu, dec, newcolumn);

      DEC_EDGE childMember = dec->edges[representativeEdge].childMember;
      DEC_NODE tail = -1;
      DEC_NODE head = -1;
      if (childMember < 0)
      {
        TU_CALL( createEdgeParallel(tu, dec, representativeEdge, &childMember) );
        debugDot(tu, dec, newcolumn);
      }
      else if (dec->members[childMember].type == DEC_MEMBER_TYPE_RIGID)
      {
        tail = findEdgeTail(dec, dec->members[childMember].markerToParent);
        head = findEdgeHead(dec, dec->members[childMember].markerToParent);
      }

      assert(reducedComponent->numTerminals == 0);
      TU_CALL( addTerminal(tu, dec, reducedComponent, childMember, tail) );
      TU_CALL( addTerminal(tu, dec, reducedComponent, childMember, head) );
      debugDot(tu, dec, newcolumn);
    }
    else if (numOneEnd == 1)
    {
      /* Root series member containing one end. */

      if (newcolumn->edgesInPath[dec->members[member].markerToParent])
      {
        /* Tested in UpdateRootSeriesOneSingleChildParent. */

        /* There is more than 1 path edge, so we split off all non-path edges and work in the new series member. */
        if (reducedMember->firstPathEdge)
        {
          TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, false, NULL, NULL, &member) );
          reducedMember->member = member;
          newcolumn->edgesInPath[dec->members[member].markerToParent] = true;
        }
        newcolumn->edgesInPath[childMarkerEdges[0]] = true;
        DEC_EDGE nonPathEdge;
        TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
        newcolumn->edgesInPath[childMarkerEdges[0]] = false;

        DEC_NODE a, b, c;
        TU_CALL( createNode(dec, &a) );
        TU_CALL( createNode(dec, &b) );
        if (dec->members[member].numEdges == 3)
        {
          TU_CALL( createNode(dec, &c) );
          TU_CALL( setEdgeNodes(tu, dec, nonPathEdge, a, c) );
        }
        else
          c = a;
        TU_CALL( setEdgeNodes(tu, dec, dec->members[member].markerToParent, a, b) );
        TU_CALL( setEdgeNodes(tu, dec, childMarkerEdges[0], c, b) );
        TU_CALL( addTerminal(tu, dec, reducedComponent, member, a) );
        TU_CALL( mergeMemberIntoParent(tu, dec, dec->edges[childMarkerEdges[0]].childMember, true) );
        dec->members[member].type = DEC_MEMBER_TYPE_RIGID;
      }
      else
      {
        /* Tested in UpdateRootSeriesOneSingleChild. */

        /* Parent marker edge is not a path edge. */
        assert(reducedMember->firstPathEdge);

        /* If there is more than 1 path edge, we squeeze off by moving them to a new series member and creating a
         * parallel member to connect it to the remaining series member. */
        DEC_EDGE pathEdge;
        TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );
        if (pathEdge >= 0)
          newcolumn->edgesInPath[pathEdge] = true;
        debugDot(tu, dec, newcolumn);

        /* Unless the series member consists of only the parent marker, the child marker (containing a path end) and a
         * representative edge, we squeeze off the representative edge and the child marker. */
        if (dec->members[member].numEdges > 3)
        {
          newcolumn->edgesInPath[childMarkerEdges[0]] = true;
          TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, true, NULL, NULL, &member) );
          newcolumn->edgesInPath[childMarkerEdges[0]] = false;
          debugDot(tu, dec, newcolumn);
        }

        assert(dec->members[member].numEdges == 3);

        DEC_NODE a, b, c;
        TU_CALL( createNode(dec, &a) );
        TU_CALL( createNode(dec, &b) );
        TU_CALL( createNode(dec, &c) );
        TU_CALL( setEdgeNodes(tu, dec, dec->members[member].markerToParent, b, c) );
        TU_CALL( setEdgeNodes(tu, dec, pathEdge, a, b) );
        TU_CALL( setEdgeNodes(tu, dec, childMarkerEdges[0], c, a) );
        TU_CALL( addTerminal(tu, dec, reducedComponent, member, b) );
        TU_CALL( mergeMemberIntoParent(tu, dec, dec->edges[childMarkerEdges[0]].childMember, true) );
      }
      debugDot(tu, dec, newcolumn);
    }
    else
    {
      assert(numOneEnd == 2);

      /* If there is more than 1 path edge, we split off by moving them to a new series member and creating a parallel
       * member to connect it to the remaining series member. */
      DEC_EDGE pathEdge = -1;
      DEC_EDGE nonPathEdge = -1;
      if (reducedMember->type != TYPE_DOUBLE_CHILD)
      {
        /* Tested in UpdateRootSeriesTwoSingleChildren. */

        /* Parent marker is a non-path edge. We split off path edges if more than one. */
        TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );

        /* If there is another non-path edge, we split off the series member consisting of the two relevant child marker
         * edges and possibly the path edge. We then replace the current series member by the new one. */

        if (dec->members[member].numEdges > (pathEdge >= 0 ? 4 : 3))
        {
          if (pathEdge >= 0)
            newcolumn->edgesInPath[pathEdge] = true;
          newcolumn->edgesInPath[childMarkerEdges[0]] = true;
          newcolumn->edgesInPath[childMarkerEdges[1]] = true;
          TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, true, NULL, NULL, &member) );
          reducedMember->member = member;
        }

        nonPathEdge = dec->members[member].markerToParent;
      }
      else
      {
        /* Tested in UpdateRootSeriesTwoSingleChildrenParent. */

        assert(newcolumn->edgesInPath[dec->members[member].markerToParent]);

        if (reducedMember->firstPathEdge)
        {
          /* Parent marker is one of several path edges. */
          TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, false, NULL, NULL, &member) );
          reducedMember->member = member;
          newcolumn->edgesInPath[dec->members[member].markerToParent] = true;
        }
        pathEdge = dec->members[member].markerToParent;
        debugDot(tu, dec, newcolumn);

        if (dec->members[member].numEdges > 3)
        {
          newcolumn->edgesInPath[childMarkerEdges[0]] = true;
          newcolumn->edgesInPath[childMarkerEdges[1]] = true;
          TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
          newcolumn->edgesInPath[childMarkerEdges[0]] = false;
          newcolumn->edgesInPath[childMarkerEdges[1]] = false;
        }
        else
          nonPathEdge = -1;
      }

      TUdbgMsg(8 + 2*depth,
        "After splitting off, the (potential) path edge is %d and the (potential) non-path edge is %d.\n",
        pathEdge, nonPathEdge);
      debugDot(tu, dec, newcolumn);

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
      DEC_NODE a, b, c, d;
      TU_CALL( createNode(dec, &a) );
      TU_CALL( createNode(dec, &b) );
      if (pathEdge >= 0)
        TU_CALL( createNode(dec, &c) );
      else
        c = b;
      if (nonPathEdge >= 0)
        TU_CALL( createNode(dec, &d) );
      else
        d = a;
      TU_CALL( setEdgeNodes(tu, dec, childMarkerEdges[0], a, b) );
      TUdbgMsg(8 + 2*depth, "First child edge %d = {%d, %d}\n", childMarkerEdges[0], a, b);
      TU_CALL( setEdgeNodes(tu, dec, childMarkerEdges[1], d, c) );
      TUdbgMsg(8 + 2*depth, "Second child edge %d = {%d, %d}\n", childMarkerEdges[1], d, c);
      if (pathEdge >= 0)
      {
        TU_CALL( setEdgeNodes(tu, dec, pathEdge, b, c) );
        TUdbgMsg(8 + 2*depth, "Path edge %d = {%d, %d}\n", pathEdge, b, c);
      }
      if (nonPathEdge >= 0)
      {
        TU_CALL( setEdgeNodes(tu, dec, nonPathEdge, d, a) );
        TUdbgMsg(8 + 2*depth, "Non-path edge %d = {%d, %d}\n", nonPathEdge, d, a);
      }

      TU_CALL( mergeMemberIntoParent(tu, dec, dec->edges[childMarkerEdges[0]].childMember, true) );
      TU_CALL( mergeMemberIntoParent(tu, dec, dec->edges[childMarkerEdges[1]].childMember, true) );
      dec->members[member].type = DEC_MEMBER_TYPE_RIGID;
    }
  }
  else
  {
    /* Non-root series member. */
    assert(reducedMember->type != TYPE_CYCLE_CHILD); /* addColumnProcess should never consider such a member. */
    assert(reducedMember->type != TYPE_DOUBLE_CHILD); /* This should only happen at the root. */
    assert(reducedMember->type != TYPE_ROOT); /* We are not a root. */
    assert(reducedMember->type == TYPE_SINGLE_CHILD); /* Only remaining case. */

    if (numOneEnd == 0)
    {
      /* Tested in UpdateLeafSeries. */

      assert(numOneEnd == 0);
      assert(numTwoEnds == 0);
      assert(reducedComponent->numTerminals < 2);
      assert(reducedMember->firstPathEdge);

      TUdbgMsg(6 + 2*depth, "Non-root series member %d without one-end children.\n", member);

      /* Squeeze off all path edges by moving them to a new series member and creating a parallel member to connect it
       * to the remaining series member. */
      DEC_EDGE pathEdge = -1;
      TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );
      assert(pathEdge >= 0);
      newcolumn->edgesInPath[pathEdge] = true;

      /* If necessary, we squeeze off the non-path edges as well. */
      assert(dec->members[member].numEdges >= 3);
      DEC_EDGE nonPathEdge;
      if (dec->members[member].numEdges == 3)
      {
        nonPathEdge = dec->edges[dec->members[member].markerToParent].next;
        if (nonPathEdge == pathEdge)
          nonPathEdge = dec->edges[nonPathEdge].next;
      }
      else
      {
        /* We temporarily mark the parent edge to belong to the path. */
        DEC_EDGE markerToParent = dec->members[member].markerToParent;
        newcolumn->edgesInPath[markerToParent] = true;
        TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
        newcolumn->edgesInPath[markerToParent] = false;
      }
      assert(dec->members[member].numEdges == 3);
      debugDot(tu, dec, newcolumn);

      /* We now create the nodes of the triangle so that the path leaves it via the parent marker edge's head node. */

      DEC_NODE a, b, c;
      TU_CALL( createNode(dec, &a) );
      TU_CALL( createNode(dec, &b) );
      TU_CALL( createNode(dec, &c) );
      TU_CALL( setEdgeNodes(tu, dec, dec->members[reducedMember->member].markerToParent, a, b) );
      TU_CALL( setEdgeNodes(tu, dec, pathEdge, b, c) );
      TU_CALL( setEdgeNodes(tu, dec, nonPathEdge, c, a) );
      TU_CALL( addTerminal(tu, dec, reducedComponent, reducedMember->member, c) );

      return TU_OKAY;
    }
    else
    {
      /* Tested in UpdateInnerSeriesOneSingleChild. */
      assert(numOneEnd == 1);

      /* Squeeze off all path edges by moving them to a new series member and creating a parallel member to connect it
       * to the remaining series member. */
      DEC_EDGE pathEdge = -1;
      TUdbgMsg(8 + 2*depth, "Splitting of path edges.\n");
      TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, true, &pathEdge, NULL, NULL) );
      if (pathEdge >= 0)
        newcolumn->edgesInPath[pathEdge] = true;
      debugDot(tu, dec, newcolumn);

      /* If necessary, we squeeze off the non-path edges as well. */
      assert(dec->members[member].numEdges >= 3);
      int numNonPathEdges = dec->members[member].numEdges - 2 - (pathEdge >= 0 ? 1 : 0);
      TUdbgMsg(8 + 2*depth, "Number of non-path edges is %d.\n", numNonPathEdges);
      DEC_EDGE nonPathEdge;
      if (numNonPathEdges == 0)
        nonPathEdge = -1;
      else if (numNonPathEdges == 1)
      {
        nonPathEdge = dec->edges[dec->members[member].markerToParent].next;
        while (nonPathEdge == childMarkerEdges[0] || nonPathEdge == pathEdge)
          nonPathEdge = dec->edges[nonPathEdge].next;
      }
      else
      {
        newcolumn->edgesInPath[dec->members[member].markerToParent] = true;
        newcolumn->edgesInPath[childMarkerEdges[0]] = true;
        TU_CALL( splitSeries(tu, dec, member, newcolumn->edgesInPath, false, &nonPathEdge, NULL, NULL) );
        newcolumn->edgesInPath[dec->members[member].markerToParent] = false;
        newcolumn->edgesInPath[childMarkerEdges[0]] = false;
        debugDot(tu, dec, newcolumn);
      }

      TUdbgMsg(8 + 2*depth, "After (potential) splitting: path edge is %d and non-path edge is %d.\n",
        pathEdge, nonPathEdge);
      assert(dec->members[member].numEdges <= 4);

      /* We now create the nodes of the triangle so that the path leaves it via the parent marker edge's head node. */

      DEC_NODE a, b, c, d;
      TU_CALL( createNode(dec, &a) );
      TU_CALL( createNode(dec, &b) );
      TU_CALL( createNode(dec, &c) );
      TU_CALL( setEdgeNodes(tu, dec, dec->members[reducedMember->member].markerToParent, a, b) );
      if (dec->members[member].numEdges == 4)
      {
        TU_CALL( createNode(dec, &d) );
        TU_CALL( setEdgeNodes(tu, dec, pathEdge, b, c) );
        TU_CALL( setEdgeNodes(tu, dec, childMarkerEdges[0], d, c) );
        TU_CALL( setEdgeNodes(tu, dec, nonPathEdge, d, a) );
      }
      else if (nonPathEdge == -1)
      {
        TU_CALL( setEdgeNodes(tu, dec, pathEdge, b, c) );
        TU_CALL( setEdgeNodes(tu, dec, childMarkerEdges[0], a, c) );
      }
      else
      {
        assert(pathEdge == -1);
        TU_CALL( setEdgeNodes(tu, dec, childMarkerEdges[0], c, b) );
        TU_CALL( setEdgeNodes(tu, dec, nonPathEdge, a, c) );
      }
      debugDot(tu, dec, newcolumn);

      TU_CALL( mergeMemberIntoParent(tu, dec, dec->edges[childMarkerEdges[0]].childMember, true) );
      debugDot(tu, dec, newcolumn);

      return TU_OKAY;
    }
  }

  return TU_OKAY;
}

/**
 * \brief Process reduced t-decomposition before the actual modification.
 *
 * Processes the reduced members in depth-first search manner and does the following:
 * - Seriess are squeezed.
 * - Terminal nodes and (reduced) members are detected.
 * - Marker edges along unique path between terminal nodes are merged.
 */

static
TU_ERROR addColumnProcessComponent(
  TU* tu,                             /**< \ref TU environment. */
  DEC* dec,                      /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn,       /**< new-column structure. */
  ReducedComponent* reducedComponent, /**< Reduced component. */
  ReducedMember* reducedMember,       /**< Reduced member. */
  int depth                           /**< Depth of member in reduced t-decomposition. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);
  assert(reducedComponent);

  TUdbgMsg(6 + 2*depth, "addColumnProcessComponent(member %d = reduced member %ld)\n", reducedMember->member,
    (reducedMember - &newcolumn->reducedMembers[0]));

  TUconsistencyAssert( decConsistency(dec) );

  /* If we are non-root type 1, then we don't need to do anything. */
  if (reducedMember->type == TYPE_CYCLE_CHILD && depth > 0)
  {
    return TU_OKAY;
  }

  /* Handle children recursively. */
  for (int c = 0; c < reducedMember->numChildren; ++c)
  {
    ReducedMember* child = reducedMember->children[c];
    if (child->type != TYPE_CYCLE_CHILD)
    {
      TU_CALL( addColumnProcessComponent(tu, dec, newcolumn, reducedComponent, reducedMember->children[c], depth+1) );
    }
    else
    {
      TUdbgMsg(8 + 2*depth, "Member %d is implicitly replaced by a path edge.\n", findMember(dec, child->member));
    }
  }

  /* Different behavior for parallel members, series members and rigid components. */
  if (dec->members[reducedMember->member].type == DEC_MEMBER_TYPE_PARALLEL)
    TU_CALL( addColumnProcessParallel(tu, dec, newcolumn, reducedComponent, reducedMember, depth) );
  else if (dec->members[reducedMember->member].type == DEC_MEMBER_TYPE_SERIES)
    TU_CALL( addColumnProcessSeries(tu, dec, newcolumn, reducedComponent, reducedMember, depth) );
  else
    TU_CALL( addColumnProcessRigid(tu, dec, newcolumn, reducedComponent, reducedMember, depth) );

  TUconsistencyAssert( decConsistency(dec) );

  return TU_OKAY;
}

static
TU_ERROR doReorderComponent(
  TU* tu,                           /**< \ref TU environment. */
  DEC* dec,                    /**< t-decomposition. */
  DEC_MEMBER member,            /**< Member to be processed. */
  DEC_MEMBER newParent,         /**< New parent of \p member. */
  DEC_MEMBER newMarkerToParent, /**< New marker edge linked to new parent of \p member. */
  DEC_MEMBER markerOfNewParent  /**< Counterpart to \p newMarkerToParent. */
)
{
  assert(tu);
  assert(dec);
  assert(member >= 0);
  assert(newParent >= 0);

  DEC_MEMBER oldParent = findMemberParent(dec, member);
  DEC_EDGE oldMarkerToParent = dec->members[member].markerToParent;
  DEC_EDGE oldMarkerOfParent = dec->members[member].markerOfParent;

  TUdbgMsg(8, "Flipping parenting of member %d with old parent %d and new parent %d.\n", member, oldParent, newParent);

  dec->members[member].markerToParent = newMarkerToParent;
  dec->members[member].markerOfParent = markerOfNewParent;
  dec->members[member].parentMember = newParent;
  dec->edges[markerOfNewParent].childMember = member;
  dec->edges[newMarkerToParent].childMember = -1;

  if (oldMarkerToParent >= 0)
    TU_CALL( doReorderComponent(tu, dec, oldParent, member, oldMarkerOfParent, oldMarkerToParent) );

  return TU_OKAY;
}

static
TU_ERROR reorderComponent(
  TU* tu,                 /**< \ref TU environment. */
  DEC* dec,          /**< t-decomposition. */
  DEC_MEMBER newRoot  /**< The new root of the component. */
)
{
  assert(tu);
  assert(dec);
  assert(newRoot >= 0 && newRoot < dec->memMembers);
  assert(isRepresentativeMember(dec, newRoot));

  TUdbgMsg(4, "Making member %d the new root of its component.\n", newRoot);

  if (dec->members[newRoot].parentMember >= 0)
  {
    TU_CALL( doReorderComponent(tu, dec, findMemberParent(dec, newRoot), newRoot,
      dec->members[newRoot].markerOfParent, dec->members[newRoot].markerToParent) );
  }

  return TU_OKAY;
}

TU_ERROR addColumnApply(
  TU* tu,                     /**< \ref TU environment. */
  DEC* dec,                  /**< t-decomposition. */
  DEC_NEWCOLUMN* newcolumn,  /**< new-column structure. */
  int column,                 /**< Index of new column to be added. */
  int* entryRows,             /**< Array of rows with 1-entry in this column. */
  int numEntries              /**< Number of 1-entries in this column. */
)
{
  assert(tu);
  assert(dec);
  assert(newcolumn);
  assert(newcolumn->remainsGraphic);

  TUdbgMsg(0, "\n  Adding a column with %d 1's.\n", numEntries);

  TUconsistencyAssert( decConsistency(dec) );

  /* Create reduced components for new edges. */
  TU_CALL( completeReducedDecomposition(tu, dec, newcolumn, entryRows, numEntries) );

  /* Process all reduced components individually. */

  DEC_EDGE* componentNewEdges = NULL;
  TU_CALL( TUallocStackArray(tu, &componentNewEdges, newcolumn->numReducedComponents) );

  int maxDepthComponent = -1;
  TUdbgMsg(4, "Processing %d reduced components.\n", newcolumn->numReducedComponents);
  for (int i = 0; i < newcolumn->numReducedComponents; ++i)
  {
    ReducedComponent* reducedComponent = &newcolumn->reducedComponents[i];

    TUdbgMsg(4, "Moving root of reduced component %d.\n", i);

    TU_CALL( moveReducedRoot(tu, dec, newcolumn, reducedComponent) );
    debugDot(tu, dec, newcolumn);

    TUdbgMsg(4, "Processing reduced component %d of depth %d.\n", i, reducedComponent->rootDepth);

    if (maxDepthComponent < 0 || reducedComponent->rootDepth
      > newcolumn->reducedComponents[maxDepthComponent].rootDepth)
    {
      maxDepthComponent = i;
    }

    TU_CALL( addColumnProcessComponent(tu, dec, newcolumn, reducedComponent, reducedComponent->root, 0) );

    assert(reducedComponent->numTerminals == 2);
    TUdbgMsg(6, "Terminal members are %d and %d.\n", reducedComponent->terminalMember[0],
      reducedComponent->terminalMember[1]);
    TUdbgMsg(6, "Terminal nodes are %d and %d.\n", reducedComponent->terminalNode[0],
      reducedComponent->terminalNode[1]);
    assert(findMember(dec, reducedComponent->terminalMember[0])
      == findMember(dec, reducedComponent->terminalMember[1]));

    /* Create new edge for this component. If there is one component, this is a column edge, and otherwise it is a
     * marker edge that will be linked to a new series member consisting of all these marker edges and the column edge.
     */
    DEC_EDGE newEdge;
    TU_CALL( createEdge(dec, -1, &newEdge) );
    componentNewEdges[i] = newEdge;
    dec->edges[newEdge].childMember = -1;
    dec->edges[newEdge].member = findMember(dec, reducedComponent->terminalMember[0]);
    dec->edges[newEdge].head = reducedComponent->terminalNode[0];
    dec->edges[newEdge].tail = reducedComponent->terminalNode[1];
    dec->edges[newEdge].name = INT_MAX/2;
    TU_CALL( addEdgeToMembersEdgeList(dec, newEdge) );
  }

  for (int c = 0; c < newcolumn->numReducedComponents; ++c)
  {
    if (c == maxDepthComponent)
      TUdbgMsg(6, "Reduced component %d has maximum depth and will remain a root.\n", c);
    else
      TU_CALL( reorderComponent(tu, dec, findMember(dec, dec->edges[componentNewEdges[c]].member)) );
  }

  if (newcolumn->numReducedComponents == 0)
  {
    DEC_MEMBER loopMember;
    TU_CALL( createMember(dec, DEC_MEMBER_TYPE_LOOP, &loopMember) );
    DEC_MEMBER loopEdge;
    TU_CALL( createEdge(dec, loopMember, &loopEdge) );
    TU_CALL( addEdgeToMembersEdgeList(dec, loopEdge) );
    dec->edges[loopEdge].name = -1 - column;
    dec->edges[loopEdge].childMember = -1;
  }
  else if (newcolumn->numReducedComponents == 1)
  {
    DEC_EDGE columnEdge = componentNewEdges[0];
    dec->edges[columnEdge].name = -1 - column;
    dec->edges[columnEdge].childMember = -1;
  }
  else
  {
    /* We create another edge for the column as well as a series member containing all new edges and this one. */

    DEC_MEMBER series;
    TU_CALL( createMember(dec, DEC_MEMBER_TYPE_SERIES, &series) );

    DEC_EDGE columnEdge;
    TU_CALL( createEdge(dec, series, &columnEdge) );
    dec->edges[columnEdge].childMember = -1;
    dec->edges[columnEdge].head = -1;
    dec->edges[columnEdge].tail = -1;
    dec->edges[columnEdge].name = -1 - column;
    TU_CALL( addEdgeToMembersEdgeList(dec, columnEdge) );

    for (int i = 0; i < newcolumn->numReducedComponents; ++i)
    {
      DEC_EDGE newEdge = componentNewEdges[i];
      DEC_EDGE markerEdge;
      TU_CALL( createEdge(dec, series, &markerEdge) );
      TU_CALL( addEdgeToMembersEdgeList(dec, markerEdge) );
      dec->edges[markerEdge].head = -1;
      dec->edges[markerEdge].tail = -1;

      DEC_MEMBER partnerMember = findEdgeMember(dec, newEdge);

      if (i == maxDepthComponent)
      {
        dec->edges[markerEdge].childMember = -1;
        dec->edges[markerEdge].name = INT_MAX - dec->numMarkerPairs;
        dec->members[series].parentMember = partnerMember;
        dec->members[series].markerToParent = markerEdge;
        dec->members[series].markerOfParent = newEdge;
        dec->edges[newEdge].name = -INT_MAX + dec->numMarkerPairs;
        dec->edges[newEdge].childMember = series;
      }
      else
      {
        dec->edges[markerEdge].childMember = partnerMember;
        dec->members[partnerMember].markerOfParent = markerEdge;
        dec->members[partnerMember].markerToParent = newEdge;
        dec->members[partnerMember].parentMember = series;
        dec->edges[markerEdge].name = INT_MAX - dec->numMarkerPairs;
        dec->edges[newEdge].name = -INT_MAX + dec->numMarkerPairs;
      }

      dec->numMarkerPairs++;
    }
  }

  debugDot(tu, dec, newcolumn);

  TU_CALL( TUfreeStackArray(tu, &componentNewEdges) );

  newcolumn->numReducedMembers = 0;
  newcolumn->numReducedComponents = 0;

  TUconsistencyAssert( decConsistency(dec) );

  return TU_OKAY;
}

TU_ERROR TUtestGraphicness(TU* tu, TU_CHRMAT* transpose, bool* pisGraphic, TU_GRAPH** pgraph, TU_GRAPH_EDGE** pbasis,
  TU_GRAPH_EDGE** pcobasis, TU_SUBMAT** psubmatrix)
{
  assert(tu);
  assert(transpose);
  assert(!psubmatrix || !*psubmatrix);
  assert(!pbasis || pgraph);
  assert(!pcobasis || pgraph);

#if defined(TU_DEBUG)
  TUdbgMsg(0, "TUtestGraphicness called a %dx%d matrix whose transpose is \n", transpose->numColumns, transpose->numRows);
  TUchrmatPrintDense(stdout, (TU_CHRMAT*) transpose, '0', true);
#endif /* TU_DEBUG */

  TU_GRAPH* graph = NULL;
  if (pgraph)
  {
    if (*pgraph)
    {
      TU_CALL( TUgraphClear(tu, *pgraph) );
    }
    else
    {
      TU_CALL( TUgraphCreateEmpty(tu, pgraph, transpose->numColumns + 2 * transpose->numRows,
        transpose->numColumns + 3 * transpose->numRows) );
    }
    graph = *pgraph;
  }

  int* basis = NULL;
  if (pbasis)
  {
    if (!*pbasis)
      TU_CALL( TUallocBlockArray(tu, pbasis, transpose->numColumns) );
    basis = *pbasis;
  }
  int* cobasis = NULL;
  if (pcobasis)
  {
    if (!*pcobasis)
      TU_CALL( TUallocBlockArray(tu, pcobasis, transpose->numRows) );
    cobasis = *pcobasis;
  }

  if (transpose->numNonzeros == 0)
  {
    *pisGraphic = true;
    if (graph)
    {
      /* Construct a path with numRows edges and with numColumns loops at 0. */

      TU_GRAPH_NODE s;
      TU_CALL( TUgraphAddNode(tu, graph, &s) );
      for (int c = 0; c < transpose->numRows; ++c)
      {
        TU_GRAPH_EDGE e;
        TU_CALL( TUgraphAddEdge(tu, graph, s, s, &e) );
        if (cobasis)
          *cobasis++ = e;
      }
      for (int r = 0; r < transpose->numColumns; ++r)
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

  DEC* dec = NULL;
  TU_CALL( decCreate(tu, &dec, 4096, 1024, 256, 256, 256) );

  /* Process each column. */
  DEC_NEWCOLUMN* newcolumn = NULL;
  TU_CALL( newcolumnCreate(tu, &newcolumn) );
  *pisGraphic = true;
  for (int column = 0; column < transpose->numRows && *pisGraphic; ++column)
  {
    TU_CALL( addColumnCheck(tu, dec, newcolumn,
      &transpose->entryColumns[transpose->rowStarts[column]],
      transpose->rowStarts[column+1] - transpose->rowStarts[column]) );

    debugDot(tu, dec, newcolumn);

    if (newcolumn->remainsGraphic)
    {
      TU_CALL( addColumnApply(tu, dec, newcolumn, column, &transpose->entryColumns[transpose->rowStarts[column]],
        transpose->rowStarts[column+1] - transpose->rowStarts[column]) );

      debugDot(tu, dec, newcolumn);
    }
    else
    {
      *pisGraphic = false;
    }
  }

  TU_CALL( newcolumnFree(tu, &newcolumn) );

  if (*pisGraphic && graph)
  {
    /* Add members and edges for empty rows. */
    if (dec->numRows < transpose->numColumns)
    {
      /* Reallocate if necessary. */
      if (dec->memRows < transpose->numColumns)
      {
        TUreallocBlockArray(tu, &dec->rowEdges, transpose->numColumns);
        dec->memRows = transpose->numColumns;
      }

      /* Add single-edge parallel for each missing row. */
      for (int r = dec->numRows; r < transpose->numColumns; ++r)
      {
        DEC_MEMBER member;
        TU_CALL( createMember(dec, DEC_MEMBER_TYPE_PARALLEL, &member) );

        DEC_EDGE edge;
        TU_CALL( createEdge(dec, member, &edge) );
        TU_CALL( addEdgeToMembersEdgeList(dec, edge) );
        dec->edges[edge].name = r;
        dec->edges[edge].head = -1;
        dec->edges[edge].tail = -1;
        dec->edges[edge].childMember = -1;

        TUdbgMsg(8, "New empty row %d is edge %d of member %d.\n", r, edge, member);

        dec->rowEdges[r].edge = edge;
      }

      dec->numRows = transpose->numColumns;
    }

    TU_CALL( decToGraph(tu, dec, graph, true, basis, cobasis, NULL) );
  }

  TU_CALL( decFree(tu, &dec) );

  return TU_OKAY;
}

typedef enum
{
  UNKNOWN = 0,    /**< The node was not considered by the shortest-path, yet. */
  SEEN = 1,       /**< Some path to the node is known. */
  COMPLETED = 2,  /**< The shortest path to the node is known. */
  BASIC = 3,      /**< The rootEdge of that node belongs to the spanning forest. */
} STAGE;

/**
 * \brief Node information for shortest-path computation in \ref TUconvertGraphToBinaryMatrix.
 */

typedef struct
{
  STAGE stage;            /**< At which stage of the algorithm is this node? */
  int predecessor;        /**< Predecessor node in shortest-path branching, or -1 for a root. */
  TU_GRAPH_EDGE rootEdge; /**< The actual edge towards the predecessor, or -1 for a root./ */
} NodeData;

/**
 * \brief Comparator for sorting ints increasingly.
 */

int compareInt(const void* A, const void* B)
{
  int* a = (int*) A;
  int* b = (int*) B;
  return *a - *b;
}

TU_ERROR TUconvertGraphToBinaryMatrix(TU* tu, TU_GRAPH* graph, TU_CHRMAT** pmatrix,
  int numBasisEdges, TU_GRAPH_EDGE* basisEdges, int numCobasisEdges, TU_GRAPH_EDGE* cobasisEdges)
{
  assert(tu);
  assert(graph);
  assert(pmatrix);
  assert(!*pmatrix);
  assert(numBasisEdges == 0 || basisEdges);
  assert(numCobasisEdges == 0 || cobasisEdges);
  TUassertStackConsistency(tu);

  NodeData* nodeData = NULL;
  
  TU_CALL( TUallocStackArray(tu, &nodeData, TUgraphMemNodes(graph)) );
  
  TU_INTHEAP heap;
  TU_CALL( TUintheapInitStack(tu, &heap, TUgraphMemNodes(graph)) );
  int* lengths = NULL;
  TU_CALL( TUallocStackArray(tu, &lengths, TUgraphMemEdges(graph)) );
  for (TU_GRAPH_NODE v = TUgraphNodesFirst(graph); TUgraphNodesValid(graph, v);
    v = TUgraphNodesNext(graph, v))
  {
    nodeData[v].stage = UNKNOWN;
  }
  for (TU_GRAPH_ITER i = TUgraphEdgesFirst(graph); TUgraphEdgesValid(graph, i);
    i = TUgraphEdgesNext(graph, i))
  {
    TU_GRAPH_EDGE e = TUgraphEdgesEdge(graph, i);
    lengths[e] = 1;
  }
  for (int b = 0; b < numBasisEdges; ++b)
  {
    TUdbgMsg(0, "basis element %d is edge %d = {%d,%d}\n", b, basisEdges[b], TUgraphEdgeU(graph, basisEdges[b]),
      TUgraphEdgeV(graph, basisEdges[b]));
    lengths[basisEdges[b]] = 0;
  }

  TUassertStackConsistency(tu);

  /* Start Dijkstra's algorithm at each node. */
  int countComponents = 0;
  for (TU_GRAPH_NODE s = TUgraphNodesFirst(graph); TUgraphNodesValid(graph, s);
    s = TUgraphNodesNext(graph, s))
  {
    if (nodeData[s].stage != UNKNOWN)
      continue;

    TUdbgMsg(2, "Executing Dijkstra at starting node %d.\n", s);
    nodeData[s].predecessor = -1;
    nodeData[s].rootEdge = -1;
    ++countComponents;
    TUintheapInsert(&heap, s, 0);
    while (!TUintheapEmpty(&heap))
    {
      int distance = TUintheapMinimumValue(&heap);
      TU_GRAPH_NODE v = TUintheapExtractMinimum(&heap);
      TUdbgMsg(4, "Processing node %d at distance %d.\n", v, distance);
      nodeData[v].stage = COMPLETED;
      for (TU_GRAPH_ITER i = TUgraphIncFirst(graph, v); TUgraphIncValid(graph, i);
        i = TUgraphIncNext(graph, i))
      {
        assert(TUgraphIncSource(graph, i) == v);
        TU_GRAPH_NODE w = TUgraphIncTarget(graph, i);

        /* Skip if already completed. */
        if (nodeData[w].stage == COMPLETED)
          continue;

        TU_GRAPH_EDGE e = TUgraphIncEdge(graph, i);
        int newDistance = distance + lengths[e];
        if (newDistance < TUintheapGetValueInfinity(&heap, w))
        {
          TUdbgMsg(6, "Updating distance of %d via predecessor %d from %d to %d.\n", w, v,
            TUintheapGetValueInfinity(&heap, w), newDistance);
          nodeData[w].stage = SEEN;
          nodeData[w].predecessor = v;
          nodeData[w].rootEdge = e;
          TUintheapDecreaseInsert(&heap, w, newDistance);
        }
      }
    }
  }

  TUassertStackConsistency(tu);
  TU_CALL( TUfreeStackArray(tu, &lengths) );
  TU_CALL( TUintheapClearStack(tu, &heap) );

  /* Now nodeData[.].predecessor is an arborescence for each connected component. */

  TU_GRAPH_NODE* nodesRows = NULL; /* Non-root node v is mapped to row of edge {v,predecessor(v)}. */
  TU_CALL( TUallocStackArray(tu, &nodesRows, TUgraphMemNodes(graph)) );
  for (TU_GRAPH_NODE v = TUgraphNodesFirst(graph); TUgraphNodesValid(graph, v); v = TUgraphNodesNext(graph, v))
  {
    nodesRows[v] = -1;
  }
  int numRows = 0;
  for (int i = 0; i < numBasisEdges; ++i)
  {
    TU_GRAPH_NODE u = TUgraphEdgeU(graph, basisEdges[i]);
    TU_GRAPH_NODE v = TUgraphEdgeV(graph, basisEdges[i]);
    TUdbgMsg(2, "Basic edge %d = {%d,%d}.\n", basisEdges[i], u, v);
    if (nodeData[u].predecessor == v)
    {
      nodesRows[u] = numRows;
      ++numRows;
      nodeData[u].stage = BASIC;
      TUdbgMsg(2, "Basic edge {%d,%d}: %d is predecessor of %d; node %d is row %d.\n", u, v, v, u, u, nodesRows[u]);
    }
    else if (nodeData[v].predecessor == u)
    {
      nodesRows[v] = numRows;
      ++numRows;
      nodeData[v].stage = BASIC;
      TUdbgMsg(2, "Basic edge {%d,%d}: %d is predecessor of %d; node %d is row %d.\n", u, v, u, v, v, nodesRows[v]);
    }
  }
  if (numRows < TUgraphNumNodes(graph) - countComponents)
  {
    /* Some basis edges are missing. */
    for (TU_GRAPH_NODE v = TUgraphNodesFirst(graph); TUgraphNodesValid(graph, v);
      v = TUgraphNodesNext(graph, v))
    {
      if (nodeData[v].predecessor >= 0 && nodeData[v].stage != BASIC)
      {
        nodesRows[v] = numRows;
        ++numRows;
        nodeData[v].stage = BASIC;
        TUdbgMsg(2, "Predecessor edge {%d,%d} not basic; node %d is row %d.\n", nodeData[v].predecessor, v, v,
          nodesRows[v]);
      }
    }
  }

  TUassertStackConsistency(tu);

  TU_CHRMAT* transposed = NULL;
  TU_CALL( TUchrmatCreate(tu, &transposed, TUgraphNumEdges(graph) - numRows, numRows,
    16 * numRows) );
  int numNonzeros = 0; /* Current number of nonzeros. transpose->numNonzeros is the memory. */
  int numColumns = 0;
  TU_GRAPH_EDGE* edgeColumns = NULL;
  TU_CALL( TUallocStackArray(tu, &edgeColumns, TUgraphMemEdges(graph)) );
  for (TU_GRAPH_ITER i = TUgraphEdgesFirst(graph); TUgraphEdgesValid(graph, i);
    i = TUgraphEdgesNext(graph, i))
  {
    TU_GRAPH_EDGE e = TUgraphEdgesEdge(graph, i);
    TU_GRAPH_NODE u = TUgraphEdgeU(graph, e);
    TU_GRAPH_NODE v = TUgraphEdgeV(graph, e);
    edgeColumns[TUgraphEdgesEdge(graph, i)] =
      (nodeData[u].rootEdge == e || nodeData[v].rootEdge == e) ? -1 : -2;
  }
  TU_GRAPH_NODE* uPath = NULL;
  TU_CALL( TUallocStackArray(tu, &uPath, numRows) );
  TU_GRAPH_NODE* vPath = NULL;
  TU_CALL( TUallocStackArray(tu, &vPath, numRows) );
  TU_GRAPH_ITER iter = TUgraphEdgesFirst(graph);
  int cobasicIndex = 0;
  while (TUgraphEdgesValid(graph, iter))
  {
    TU_GRAPH_EDGE e;
    if (cobasicIndex < numCobasisEdges)
    {
      e = cobasisEdges[cobasicIndex];
      ++cobasicIndex;
    }
    else
    {
      e = TUgraphEdgesEdge(graph, iter);
      iter = TUgraphEdgesNext(graph, iter);
    }

    if (edgeColumns[e] >= -1)
      continue;

    TU_GRAPH_NODE u = TUgraphEdgeU(graph, e);
    TU_GRAPH_NODE v = TUgraphEdgeV(graph, e);

    transposed->rowStarts[numColumns] = numNonzeros;
    edgeColumns[e] = numColumns;

    /* Enlarge space for nonzeros if necessary. */
    if (numNonzeros + numRows > transposed->numNonzeros)
      TU_CALL( TUchrmatChangeNumNonzeros(tu, transposed, 2 * transposed->numNonzeros) );

    /* Compute u-root path. */
    int uPathLength = 0;
    TU_GRAPH_NODE w = u;
    while (nodeData[w].predecessor != -1)
    {
      uPath[uPathLength] = w;
      ++uPathLength;
      w = nodeData[w].predecessor;
    }

    /* Compute v-root path. */
    int vPathLength = 0;
    w = v;
    while (nodeData[w].predecessor != -1)
    {
      vPath[vPathLength] = w;
      ++vPathLength;
      w = nodeData[w].predecessor;
    }

    /* Remove common part of u-root path and v-root path. */
    while (uPathLength > 0 && vPathLength > 0 && uPath[uPathLength-1] == vPath[vPathLength-1])
    {
      --uPathLength;
      --vPathLength;
    }

    for (int j = 0; j < uPathLength; ++j)
    {
      assert(nodesRows[uPath[j]] >= 0);
      transposed->entryColumns[numNonzeros] = nodesRows[uPath[j]];
      transposed->entryValues[numNonzeros] = 1;
      ++numNonzeros;
    }
    for (int j = 0; j < vPathLength; ++j)
    {
      assert(nodesRows[vPath[j]] >= 0);
      transposed->entryColumns[numNonzeros] = nodesRows[vPath[j]];
      transposed->entryValues[numNonzeros] = 1;
      ++numNonzeros;
    }
    qsort(&transposed->entryColumns[transposed->rowStarts[numColumns]], uPathLength + vPathLength,
      sizeof(int), compareInt);

    ++numColumns;
  }

  TUassertStackConsistency(tu);

  TU_CALL( TUfreeStackArray(tu, &vPath) );
  TU_CALL( TUfreeStackArray(tu, &uPath) );
  TU_CALL( TUfreeStackArray(tu, &edgeColumns) );

  transposed->rowStarts[numColumns] = numNonzeros;
  if (numNonzeros == 0 && transposed->numNonzeros > 0)
  {
    TU_CALL( TUfreeBlockArray(tu, &transposed->entryColumns) );
    TU_CALL( TUfreeBlockArray(tu, &transposed->entryValues) );
  }
  transposed->numNonzeros = numNonzeros;

  TU_CALL( TUchrmatTranspose(tu, transposed, pmatrix) );
  TU_CALL( TUchrmatFree(tu, &transposed) );

  /* We now process the nonbasic edges. */

  TUassertStackConsistency(tu);
  TU_CALL( TUfreeStackArray(tu, &nodesRows) );
  TUassertStackConsistency(tu);
  TU_CALL( TUfreeStackArray(tu, &nodeData) );

  TUassertStackConsistency(tu);

  return TU_OKAY;
}
