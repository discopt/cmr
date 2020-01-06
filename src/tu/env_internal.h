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

/**
 * \brief Allocates a memory block on the stack.
 */

#define TUallocStack(tu, ptr) \
  _TUallocStack(tu, (void**) ptr, sizeof(**ptr))

TU_EXPORT
void _TUallocStack(TU* tu, void** ptr, size_t size);

/**
 * \brief Frees a memory block from the stack, which must be the top memory chunk.
 */

#define TUfreeStack(tu, ptr) \
  _TUfreeStack(tu, (void**) ptr, sizeof(**ptr))

TU_EXPORT
void _TUfreeStack(TU* tu, void** ptr, size_t size);


/**
 * \brief Allocates memory for an array of blocks on the stack.
 */

#define TUallocStackArray(tu, ptr, length) \
  _TUallocStackArray(tu, (void**) ptr, sizeof(**ptr), length)

TU_EXPORT
void _TUallocStackArray(TU* tu, void** ptr, size_t size, size_t length);

/**
 * \brief Frees memory of an array of blocks on the stack. Must be the top memory chunk.
 */

#define TUfreeStackArray(tu, ptr) \
  _TUfreeStackArray(tu, (void**) ptr)

TU_EXPORT
void _TUfreeStackArray(TU* tu, void** ptr);

#endif /* TU_ENV_INTERNAL_H */
