#ifndef TU_ENV_H
#define TU_ENV_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/config.h>
#include <tu/export.h>
#include <stdbool.h>

struct _TU;

/**
 * \brief Environment for computations
 * 
 * Manages memory, threading, output and parameters.
 */

typedef struct _TU TU;

/**
 * \brief Initializes a default TU
 * 
 * It has default parameters and outputs to stdout.
 */

TU_EXPORT
void TUinit(
  TU** tu /** Pointer to TU environment */
);

/**
 * \brief Frees a TU environment
 */

TU_EXPORT
void TUfree(
  TU** tu /** Pointer to TU environment */
);

#ifdef __cplusplus
}
#endif

#endif /* TU_ENV_H */

