#include <tu/graph.h>

#include <assert.h>

void TUfreeGraph(TU* tu, TU_GRAPH** pgraph)
{
  assert(tu);
  assert(pgraph);

  TUfreeBlock(tu, pgraph);
  *pgraph = NULL;
}
