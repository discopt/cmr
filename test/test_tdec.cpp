#include <gtest/gtest.h>

#include "common.h"
#include <tu/graphic.h>

TEST(TDec, SingleColumn)
{
  TU* tu = NULL;
  TUcreateEnvironment(&tu);

  TU_CHRMAT* matrix = NULL;
  stringToCharMatrix(tu, &matrix, "6 1 "
    "1 "
    "1 "
    "1 "
    "1 "
    "1 "
    "1 "
  );

  TU_GRAPH* graph = NULL;
  bool isGraphic;
  TUtestGraphicnessChr(tu, matrix, &isGraphic, &graph, NULL, NULL, NULL);
  ASSERT_TRUE(isGraphic);
  ASSERT_TRUE(graph);
  
  TUgraphPrint(stdout, graph);
  
  TUgraphFree(tu, &graph);

  TUchrmatFree(tu, &matrix);
  TUfreeEnvironment(&tu);
}
