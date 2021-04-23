#ifndef TU_ENV_INTERNAL_H
#define TU_ENV_INTERNAL_H

#include <stdio.h>
#include <stdbool.h>
#include <stdarg.h>

#if defined(TU_DEBUG)

static
void TUdbgMsg(int indent, const char* format, ...)
{
  va_list args;

  for (int i = 0; i < indent; ++i)
    putchar(' ');
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
  fflush(stdout);
}

#else /* TU_DEBUG */

#define TUdbgMsg(...)

#endif /* TU_DEBUG */


typedef struct
{
  char* memory; /**< \brief Raw memory. */
  size_t top;   /**< \brief First used byte. */
} TU_STACK;

struct TU_ENVIRONMENT
{
  FILE* output;         /**< \brief Output stream or \c NULL if silent. */
  bool closeOutput;     /**< \brief Whether to close the output stream at the end. */
  int verbosity;        /**< \brief Verbosity level. */
  int numThreads;       /**< \brief Number of threads to use. */

  size_t numStacks;     /**< \brief Number of allocated stacks in stack array. */
  size_t memStacks;     /**< \brief Memory for stack array. */
  size_t currentStack;  /**< \brief Index of last used stack. */
  TU_STACK* stacks;     /**< \brief Array of stacks. */
};

#include <tu/env.h>

/**
 * \brief Allocates statck memory for *\p ptr.
 *
 * Stack memory shall be freed with \ref TUfreeStack in the reverse order of allocation.
 * The size is determined automatically.
 */

#define TUallocStack(tu, ptr) \
  _TUallocStack(tu, (void**) ptr, sizeof(**ptr))

/**
 * \brief Carries out the allocation for \ref TUallocStack.
 * 
 * \note Use \ref TUallocStack to allocate stack memory.
 */

TU_EXPORT
TU_ERROR _TUallocStack(
  TU* tu,     /**< \ref TU environment. */
  void** ptr, /**< Pointer where the space shall be allocated. */
  size_t size /**< Space to allocate. */
);

/**
 * \brief Frees a stack memory chunk allocated with \ref TUallocStack.
 */

#define TUfreeStack(tu, ptr) \
  _TUfreeStack(tu, (void**) ptr, sizeof(**ptr))

/**
 * \brief Carries out the deallocation for \ref TUfreeStack.
 * 
 * \note Use \ref TUfreeStack to free stack memory.
 */

TU_EXPORT
TU_ERROR _TUfreeStack(
  TU* tu,     /**< \ref TU environment. */
  void** ptr  /**< Pointer of space to be freed. */
);

/**
 * \brief Allocates memory for an array of blocks on the stack.
 */

#define TUallocStackArray(tu, ptr, length) \
  _TUallocStack(tu, (void**) ptr, sizeof(**ptr) * (length))

/**
 * \brief Frees memory of an array of blocks on the stack.
 */
  
#define TUfreeStackArray(tu, ptr) \
  _TUfreeStack(tu, (void**) ptr)

#if !defined(NDEBUG)

/**
 * \brief Checks stack protection fields for corruption. 
 *
 * Useful for debugging memory errors.
 */

void TUassertStackConsistency(
  TU* tu  /**< \ref TU environment. */
);

#else

static inline
void TUassertStackConsistency(
  TU* tu  /**< \ref TU environment. */
)
{
  assert(tu);
}


#endif /* !NDEBUG */
  
char* TUconsistencyMessage(const char* format, ...);

#if !defined(NDEBUG)

#define TUconsistencyAssert( call ) \
  do \
  { \
    char* __message = call; \
    if (__message) \
    { \
      fflush(stdout); \
      fprintf(stderr, "%s:%d: %s\n", __FILE__, __LINE__, __message); \
      fflush(stderr); \
      free(__message); \
      exit(1); \
    } \
  } \
  while (false);

#else

#define TUconsistencyAssert( call )
  
#endif

#endif /* TU_ENV_INTERNAL_H */
