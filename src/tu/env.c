#include "env_internal.h"

#include <assert.h>
#include <stdlib.h>

void TUinit(
  TU** tu
  )
{
  assert(tu != NULL);
  *tu = (TU*) malloc(sizeof(TU));
  (*tu)->output = stdout;
  (*tu)->closeOutput = false;
  (*tu)->numThreads = 1;
  (*tu)->verbosity = 1;
}

/**
 * \brief Frees a TU environment
 */

void TUfree(
  TU** tu
  )
{
  assert(tu != NULL);
  if ((*tu)->closeOutput)
    fclose((*tu)->output);
  free(*tu);
}
