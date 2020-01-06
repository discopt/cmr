#ifndef TU_GRAPH_H
#define TU_GRAPH_H

#include <tu/env.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  int numNodes;
  int numArcs;
} TU_GRAPH;

void TUfreeGraph(
  TU* tu,           /**< TU environment. */
  TU_GRAPH** pgraph /**< Graph */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_GRAPH_H */
