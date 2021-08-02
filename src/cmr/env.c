// #define CMR_DEBUG /* Uncomment for general debugging. */
// #define DEBUG_STACK /* Uncomment to debug CMRallocStack and CMRfreeStack. */
// #define REPLACE_STACK_BY_MALLOC /* Uncomment to not use a stack at all, which may help to detect memory corruption. */

#include "env_internal.h"

#include <assert.h>
#include <stdlib.h>
#include <limits.h>
#include <stdarg.h>
#include <string.h>

static const size_t FIRST_STACK_SIZE = 4096L; /**< Size of the first stack. */
static const int INITIAL_MEM_STACKS = 16;     /**< Initial number of allocated stacks. */

#if !defined(NDEBUG) && !defined(REPLACE_STACK_BY_MALLOC)
static const int PROTECTION = INT_MIN / 42;   /**< Protection bytes to detect corruption. */
#endif /* !NDEBUG */

CMR_ERROR CMRcreateEnvironment(CMR** ptu)
{
  if (!ptu)
    return CMR_ERROR_INPUT;

  *ptu = (CMR*) malloc(sizeof(CMR));
  CMR* cmr = *ptu;
  if (!cmr)
    return CMR_ERROR_MEMORY;

  cmr->output = stdout;
  cmr->closeOutput = false;
  cmr->numThreads = 1;
  cmr->verbosity = 1;

  /* Initialize stack memory. */
  cmr->stacks = malloc(INITIAL_MEM_STACKS * sizeof(CMR_STACK));
  if (!cmr->stacks)
  {
    free(*ptu);
    *ptu = NULL;
    return CMR_ERROR_MEMORY;
  }
  cmr->stacks[0].memory = malloc(FIRST_STACK_SIZE * sizeof(char));
  if (!cmr->stacks[0].memory)
  {
    free(cmr->stacks);
    free(cmr);
    *ptu = NULL;
    return CMR_ERROR_MEMORY;
  }
  cmr->stacks[0].top = FIRST_STACK_SIZE;
  cmr->memStacks = INITIAL_MEM_STACKS;
  cmr->numStacks = 1;
  cmr->currentStack = 0;

  return CMR_OKAY;
}

CMR_ERROR CMRfreeEnvironment(CMR** ptu)
{
  if (!ptu)
    return CMR_ERROR_INPUT;

  CMR* cmr = *ptu;

  if (cmr->closeOutput)
    fclose(cmr->output);

  for (int s = 0; s < cmr->numStacks; ++s)
    free(cmr->stacks[s].memory);
  free(cmr->stacks);
  free(*ptu);
  *ptu = NULL;

  return CMR_OKAY;
}

CMR_ERROR _CMRallocBlock(CMR* cmr, void** ptr, size_t size)
{
  assert(cmr);
  assert(ptr);
  assert(*ptr == NULL);
  *ptr = malloc(size);

  return *ptr ? CMR_OKAY : CMR_ERROR_MEMORY;
}

CMR_ERROR _CMRfreeBlock(CMR* cmr, void** ptr, size_t size)
{
  assert(cmr);
  assert(ptr);
  assert(*ptr);
  free(*ptr);
  *ptr = NULL;

  return CMR_OKAY;
}

CMR_ERROR _CMRallocBlockArray(CMR* cmr, void** ptr, size_t size, size_t length)
{
  assert(cmr);
  assert(ptr);
  assert(*ptr == NULL);
  *ptr = malloc(size * length);

  return *ptr ? CMR_OKAY : CMR_ERROR_MEMORY;
}


CMR_ERROR _CMRreallocBlockArray(CMR* cmr, void** ptr, size_t size, size_t length)
{
  assert(cmr);
  assert(ptr);
  *ptr = realloc(*ptr, size * length);

  return *ptr ? CMR_OKAY : CMR_ERROR_MEMORY;
}

CMR_ERROR _CMRduplicateBlockArray(CMR* cmr, void** ptr, size_t size, size_t length, void* source)
{
  assert(cmr);
  assert(ptr);
  assert(source);

  CMR_CALL( _CMRallocBlockArray(cmr, ptr, size, length) );
  size_t numBytes = size*length;
  for (size_t i = 0; i < numBytes; ++i)
    ((char*)(*ptr))[i] = ((char*)source)[i];

  return CMR_OKAY;
}

CMR_ERROR _CMRfreeBlockArray(CMR* cmr, void** ptr)
{
  assert(cmr);
  assert(ptr);
  assert(*ptr);
  free(*ptr);
  *ptr = NULL;

  return CMR_OKAY;
}

#if defined(REPLACE_STACK_BY_MALLOC)

CMR_ERROR _CMRallocStack(
  CMR* cmr,
  void** ptr,
  size_t size
)
{
  assert(cmr);
  assert(ptr);

  /* Avoid allocation of zero bytes. */
  if (size < 4)
    size = 4;

  *ptr = malloc(size);

  return *ptr ? CMR_OKAY : CMR_ERROR_MEMORY;
}

CMR_ERROR _CMRfreeStack(
  CMR* cmr,
  void** ptr
)
{
  assert(cmr);
  assert(ptr);
  assert(*ptr);

  free(*ptr);

  return CMR_OKAY;
}

void CMRassertStackConsistency(
  CMR* cmr
)
{

}

#else

#define STACK_SIZE(k) \
  (FIRST_STACK_SIZE << k)

CMR_ERROR _CMRallocStack(
  CMR* cmr,
  void** ptr,
  size_t size
)
{
  assert(cmr);
  assert(ptr);
  assert(*ptr == NULL);

  /* Avoid allocation of zero bytes. */
  if (size < 4)
    size = 4;

  size_t requiredSpace = size + sizeof(void*);
#if !defined(NDEBUG)
  requiredSpace += sizeof(int);
#endif /* !NDEBUG */

#if defined(DEBUG_STACK)
  printf("CMRallocStack() called for %ld bytes; current stack: %ld, numStacks: %ld, memStack: %ld.\n",
    size, cmr->currentStack, cmr->numStacks, cmr->memStacks);
  fflush(stdout);
  printf("Current stack has capacity %ld and %ld free bytes.\n",
    FIRST_STACK_SIZE << cmr->currentStack, cmr->stacks[cmr->currentStack].top);
  fflush(stdout);
#endif /* DEBUG_STACK */

  while (cmr->stacks[cmr->currentStack].top < requiredSpace)
  {
    ++cmr->currentStack;
    if (cmr->currentStack == cmr->numStacks)
    {
      /* If necessary, enlarge the slacks array. */
      if (cmr->numStacks == cmr->memStacks)
      {
        cmr->stacks = realloc(cmr->stacks, 2 * cmr->memStacks * sizeof(CMR_STACK));
        size_t newSize = 2*cmr->memStacks;
        for (int s = cmr->memStacks; s < newSize; ++s)
        {
          cmr->stacks[s].memory = NULL;
          cmr->stacks[s].top = FIRST_STACK_SIZE << s;
        }
        cmr->memStacks = newSize;
      }

      cmr->stacks[cmr->numStacks].top = FIRST_STACK_SIZE << cmr->numStacks;
      cmr->stacks[cmr->numStacks].memory = malloc(cmr->stacks[cmr->numStacks].top * sizeof(char));
      ++cmr->numStacks;
    }

    assert(cmr->stacks[cmr->currentStack].top == (FIRST_STACK_SIZE << cmr->currentStack));
  }

  /* The chunk fits into the last stack. */

  CMR_STACK* pstack = &cmr->stacks[cmr->currentStack];
  pstack->top -= size;
  *ptr = &pstack->memory[pstack->top];
#if !defined(NDEBUG)
  pstack->top -= sizeof(int);
  *((int*) &pstack->memory[pstack->top]) = PROTECTION;
#endif /* !NDEBUG */
  pstack->top -= sizeof(void*);
  *((size_t*) &pstack->memory[pstack->top]) = size;

#if defined(DEBUG_STACK)
  printf("Writing size %ld to %p.\n", size, &pstack->memory[pstack->top]);
#endif /* DEBUG_STACK */

  return CMR_OKAY;
}

CMR_ERROR _CMRfreeStack(CMR* cmr, void** ptr)
{
  assert(cmr);
  assert(ptr);
  assert(*ptr);

  CMR_STACK* stack = &cmr->stacks[cmr->currentStack];
  CMRdbgMsg(0, "CMRfreeStack called for pointer %p. Last stack is %d.\n", *ptr, cmr->currentStack);
  size_t size = *((size_t*) &stack->memory[stack->top]);

#if defined(DEBUG_STACK)
  printf("CMRfreeStack() called for %ld bytes (size stored at %p).\n", size,
    &stack->memory[stack->top]);
  fflush(stdout);
#endif /* DEBUG_STACK */

  assert(size < (FIRST_STACK_SIZE << cmr->numStacks));

#if !defined(NDEBUG)
  if (*((int*) (&stack->memory[stack->top] + sizeof(void*))) != PROTECTION)
  {
    fprintf(stderr, "Memory corruption of stack detected!\n");
    fflush(stderr);
    assert(false);
  }
#endif /* !NDEBUG */


#ifndef NDEBUG
  if (&stack->memory[stack->top + sizeof(int) + sizeof(void*)] != *ptr)
  {
    fprintf(stderr,
      "Wrong order of CMRfreeStack(Array) detected. Top chunk on stack has size %ld!\n",
      size);
    fflush(stderr);
  }
#endif /* !NDEBUG */

  stack->top += size + sizeof(void*);
#if !defined(NDEBUG)
  stack->top += sizeof(int);
#endif /* !NDEBUG */

  while (stack->top == (FIRST_STACK_SIZE << cmr->currentStack) && cmr->currentStack > 0)
  {
    --cmr->currentStack;
    stack = &cmr->stacks[cmr->currentStack];
  }
  *ptr = NULL;

  return CMR_OKAY;
}

#if !defined(NDEBUG)

void CMRassertStackConsistency(
  CMR* cmr
)
{
  assert(cmr);

  for (int s = 0; s <= cmr->currentStack; ++s)
  {
    CMR_STACK* stack = &cmr->stacks[s];

    void* ptr = &stack->memory[stack->top];
    CMRdbgMsg(2, "Stack %d of size %d has memory range [%p,%p). top is %p\n", s, STACK_SIZE(s), stack->memory,
      stack->memory + STACK_SIZE(s), ptr);
    while (ptr < (void*)stack->memory + STACK_SIZE(s))
    {
      CMRdbgMsg(4, "pointer is %p.", ptr);
      size_t size = *((size_t*) ptr);
      CMRdbgMsg(0, " It indicates a chunk of size %d.\n", size);
      ptr += sizeof(size_t*);
      assert(*((int*)ptr) == PROTECTION);
      ptr += size + sizeof(int);
    }
  }
}

#endif /* !NDEBUG */

#endif /* else REPLACE_STACK_BY_MALLOC */

char* CMRconsistencyMessage(const char* format, ...)
{
  assert(format);

  char buffer[256];
  va_list argptr;
  va_start(argptr, format);
  vsprintf(buffer, format, argptr);
  va_end(argptr);

  return strdup(buffer);
}
