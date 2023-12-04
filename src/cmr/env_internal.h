#ifndef CMR_ENV_INTERNAL_H
#define CMR_ENV_INTERNAL_H

#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>

#include <cmr/env.h>

#if defined(CMR_DEBUG)

static inline
void CMRdbgMsg(int indent, const char* format, ...)
{
  va_list args;

  for (int i = 0; i < indent; ++i)
    putchar(' ');
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
  fflush(stdout);
}

#else /* !CMR_DEBUG */

static inline
void CMRdbgMsg(int indent, const char* format, ...)
{
  CMR_UNUSED(indent);
  CMR_UNUSED(format);
}
/*#define CMRdbgMsg(...) */

#endif /* CMR_DEBUG */


typedef struct
{
  char* memory; /**< \brief Raw memory. */
  size_t top;   /**< \brief First used byte. */
} CMR_STACK;

struct CMR_ENVIRONMENT
{
  char* errorMessage;   /**< \brief Error message. */

  FILE* output;         /**< \brief Output stream or \c NULL if silent. */
  bool closeOutput;     /**< \brief Whether to close the output stream at the end. */
  int verbosity;        /**< \brief Verbosity level. */
  int numThreads;       /**< \brief Number of threads to use. */

  size_t numStacks;     /**< \brief Number of allocated stacks in stack array. */
  size_t memStacks;     /**< \brief Memory for stack array. */
  size_t currentStack;  /**< \brief Index of last used stack. */
  CMR_STACK* stacks;     /**< \brief Array of stacks. */
};

#include <cmr/env.h>

/**
 * \brief Allocates statck memory for *\p ptr.
 *
 * Stack memory shall be freed with \ref CMRfreeStack in the reverse order of allocation.
 * The size is determined automatically.
 */

#define CMRallocStack(cmr, ptr) \
  _CMRallocStack(cmr, (void**) ptr, sizeof(**ptr))

/**
 * \brief Carries out the allocation for \ref CMRallocStack.
 *
 * \note Use \ref CMRallocStack to allocate stack memory.
 */

CMR_EXPORT
CMR_ERROR _CMRallocStack(
  CMR* cmr,     /**< \ref CMR environment. */
  void** ptr, /**< Pointer where the space shall be allocated. */
  size_t size /**< Space to allocate. */
);

/**
 * \brief Frees a stack memory chunk allocated with \ref CMRallocStack.
 */

#define CMRfreeStack(cmr, ptr) \
  _CMRfreeStack(cmr, (void**) ptr, sizeof(**ptr))

/**
 * \brief Carries out the deallocation for \ref CMRfreeStack.
 *
 * \note Use \ref CMRfreeStack to free stack memory.
 */

CMR_EXPORT
CMR_ERROR _CMRfreeStack(
  CMR* cmr,     /**< \ref CMR environment. */
  void** ptr  /**< Pointer of space to be freed. */
);

/**
 * \brief Allocates memory for an array of blocks on the stack.
 */

#define CMRallocStackArray(cmr, ptr, length) \
  _CMRallocStack(cmr, (void**) ptr, sizeof(**ptr) * (length))

/**
 * \brief Frees memory of an array of blocks on the stack.
 */

#define CMRfreeStackArray(cmr, ptr) \
  _CMRfreeStack(cmr, (void**) ptr)

#if !defined(NDEBUG)

/**
 * \brief Checks stack protection fields for corruption.
 *
 * Useful for debugging memory errors.
 */

void CMRassertStackConsistency(
  CMR* cmr  /**< \ref CMR environment. */
);

#else

static inline
void CMRassertStackConsistency(
  CMR* cmr  /**< \ref CMR environment. */
)
{
  CMR_UNUSED(cmr);
  assert(cmr);
}

#endif /* !NDEBUG */

void CMRraiseErrorMessage(
  CMR* cmr,               /**< \ref CMR environment. */
  const char* format, ... /**< \ref Variadic arguments in printf-style. */
);

size_t CMRgetStackUsage(
  CMR* cmr  /**< \ref CMR environment. */
);

char* CMRconsistencyMessage(const char* format, ...);

#if !defined(NDEBUG)

/**
 * \brief Asserts that \p call reports consistency. Otherwise, the inconsistency explanation is printed and the program
 *        terminates.
 * 
 * The following example code checks \c matrix (of type \ref CMR_CHRMAT*) for consistency using
 * \ref CMRchrmatConsistency.
 * 
 *     CMRconsistencyAssert( CMRchrmatConsistency(matrix) );
 */

#define CMRconsistencyAssert( call ) \
  do \
  { \
    char* __message = call; \
    if (__message) \
    { \
      fflush(stdout); \
      fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, __message); \
      fflush(stderr); \
      free(__message); \
      assert(!"Consistency assertion raised!"); \
    } \
  } \
  while (false);

#else

#define CMRconsistencyAssert( call )

#endif

#endif /* CMR_ENV_INTERNAL_H */
