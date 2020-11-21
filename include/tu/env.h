#ifndef TU_ENV_H
#define TU_ENV_H

#ifdef __cplusplus
extern "C" {
#endif

#include <tu/config.h>
#include <tu/export.h>
#include <stdbool.h>
#include <stdlib.h>

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
void TUcreateEnvironment(
  TU** tu /** Pointer to TU environment */
);

/**
 * \brief Frees a TU environment
 */

TU_EXPORT
void TUfreeEnvironment(
  TU** tu /** Pointer to TU environment */
);

/**
 * \brief Allocates a memory block.
 */

#define TUallocBlock(tu, ptr) \
  _TUallocBlock(tu, (void**) ptr, sizeof(**ptr))

TU_EXPORT
void _TUallocBlock(TU* tu, void** ptr, size_t size);

/**
 * \brief Frees a memory block.
 */

#define TUfreeBlock(tu, ptr) \
  _TUfreeBlock(tu, (void**) ptr, sizeof(**ptr))

TU_EXPORT
void _TUfreeBlock(TU* tu, void** ptr, size_t size);


/**
 * \brief Allocates memory for an array of blocks.
 */

#define TUallocBlockArray(tu, ptr, length) \
  _TUallocBlockArray(tu, (void**) ptr, sizeof(**ptr), length)

TU_EXPORT
void _TUallocBlockArray(TU* tu, void** ptr, size_t size, size_t length);

/**
 * \brief Reallocates memory for an array of blocks.
 */

#define TUreallocBlockArray(tu, ptr, length) \
  _TUreallocBlockArray(tu, (void**) ptr, sizeof(**ptr), length)

TU_EXPORT
void _TUreallocBlockArray(TU* tu, void** ptr, size_t size, size_t length);

/**
 * \brief Frees memory of an array of blocks.
 */

#define TUfreeBlockArray(tu, ptr) \
  _TUfreeBlockArray(tu, (void**) ptr)

TU_EXPORT
void _TUfreeBlockArray(TU* tu, void** ptr);

#ifdef __cplusplus
}
#endif

#endif /* TU_ENV_H */

