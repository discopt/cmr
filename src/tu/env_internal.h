#ifndef TU_ENV_INTERNAL_H
#define TU_ENV_INTERNAL_H

#include <stdio.h>
#include <stdbool.h>

struct _TU
{
  FILE* output; /**< Output stream or \c NULL if silent */
  bool closeOutput; /**< Whether to close the output stream at the end */
  int verbosity; /**< Verbosity level */
  int numThreads; /**< Number of threads to use */
};

#include <tu/env.h>

#endif /* TU_ENV_INTERNAL_H */
