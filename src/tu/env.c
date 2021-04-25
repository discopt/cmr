// #define TU_DEBUG /* Uncomment for general debugging. */
// #define DEBUG_STACK /* Uncomment to debug TUallocStack and TUfreeStack. */
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

TU_ERROR TUcreateEnvironment(TU** ptu)
{
  if (!ptu)
    return TU_ERROR_INPUT;

  *ptu = (TU*) malloc(sizeof(TU));
  TU* tu = *ptu;
  if (!tu)
    return TU_ERROR_MEMORY;

  tu->output = stdout;
  tu->closeOutput = false;
  tu->numThreads = 1;
  tu->verbosity = 1;

  /* Initialize stack memory. */
  tu->stacks = malloc(INITIAL_MEM_STACKS * sizeof(TU_STACK));
  if (!tu->stacks)
  {
    free(*ptu);
    *ptu = NULL;
    return TU_ERROR_MEMORY;
  }
  tu->stacks[0].memory = malloc(FIRST_STACK_SIZE * sizeof(char));
  if (!tu->stacks[0].memory)
  {
    free(tu->stacks);
    free(tu);
    *ptu = NULL;
    return TU_ERROR_MEMORY;
  }
  tu->stacks[0].top = FIRST_STACK_SIZE;
  tu->memStacks = INITIAL_MEM_STACKS;
  tu->numStacks = 1;
  tu->currentStack = 0;

  return TU_OKAY;
}

TU_ERROR TUfreeEnvironment(TU** ptu)
{
  if (!ptu)
    return TU_ERROR_INPUT;

  TU* tu = *ptu;
  
  if (tu->closeOutput)
    fclose(tu->output);

  for (int s = 0; s < tu->numStacks; ++s)
    free(tu->stacks[s].memory);
  free(tu->stacks);
  free(*ptu);
  *ptu = NULL;

  return TU_OKAY;
}

TU_ERROR _TUallocBlock(TU* tu, void** ptr, size_t size)
{
  assert(tu);
  assert(ptr);
  assert(*ptr == NULL);
  *ptr = malloc(size);

  return *ptr ? TU_OKAY : TU_ERROR_MEMORY;
}

TU_ERROR _TUfreeBlock(TU* tu, void** ptr, size_t size)
{
  assert(tu);
  assert(ptr);
  assert(*ptr);
  free(*ptr);
  *ptr = NULL;

  return TU_OKAY;
}

TU_ERROR _TUallocBlockArray(TU* tu, void** ptr, size_t size, size_t length)
{
  assert(tu);
  assert(ptr);
  assert(*ptr == NULL);
  *ptr = malloc(size * length);

  return *ptr ? TU_OKAY : TU_ERROR_MEMORY;
}

TU_ERROR _TUreallocBlockArray(TU* tu, void** ptr, size_t size, size_t length)
{
  assert(tu);
  assert(ptr);
  *ptr = realloc(*ptr, size * length);

  return *ptr ? TU_OKAY : TU_ERROR_MEMORY;
}

TU_ERROR _TUfreeBlockArray(TU* tu, void** ptr)
{
  assert(tu);
  assert(ptr);
  assert(*ptr);
  free(*ptr);
  *ptr = NULL;

  return TU_OKAY;
}

#if defined(REPLACE_STACK_BY_MALLOC)

TU_ERROR _TUallocStack(
  TU* tu,
  void** ptr,
  size_t size
)
{
  assert(tu);
  assert(ptr);

  /* Avoid allocation of zero bytes. */
  if (size < 4)
    size = 4;

  *ptr = malloc(size);

  return *ptr ? TU_OKAY : TU_ERROR_MEMORY;
}

TU_ERROR _TUfreeStack(
  TU* tu,
  void** ptr
)
{
  assert(tu);
  assert(ptr);
  assert(*ptr);

  free(*ptr);

  return TU_OKAY;
}

void TUassertStackConsistency(
  TU* tu
)
{

}

#else

#define STACK_SIZE(k) \
  (FIRST_STACK_SIZE << k)

TU_ERROR _TUallocStack(
  TU* tu,
  void** ptr,
  size_t size
)
{
  assert(tu);
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
  printf("TUallocStack() called for %ld bytes; current stack: %ld, numStacks: %ld, memStack: %ld.\n",
    size, tu->currentStack, tu->numStacks, tu->memStacks);
  fflush(stdout);
  printf("Current stack has capacity %ld and %ld free bytes.\n",
    FIRST_STACK_SIZE << tu->currentStack, tu->stacks[tu->currentStack].top);
  fflush(stdout);
#endif /* DEBUG_STACK */

  while (tu->stacks[tu->currentStack].top < requiredSpace)
  {
    ++tu->currentStack;
    if (tu->currentStack == tu->numStacks)
    {
      /* If necessary, enlarge the slacks array. */
      if (tu->numStacks == tu->memStacks)
      {
        tu->stacks = realloc(tu->stacks, 2 * tu->memStacks * sizeof(TU_STACK));
        size_t newSize = 2*tu->memStacks;
        for (int s = tu->memStacks; s < newSize; ++s)
        {
          tu->stacks[s].memory = NULL;
          tu->stacks[s].top = FIRST_STACK_SIZE << s;
        }
        tu->memStacks = newSize;
      }

      tu->stacks[tu->numStacks].top = FIRST_STACK_SIZE << tu->numStacks;
      tu->stacks[tu->numStacks].memory = malloc(tu->stacks[tu->numStacks].top * sizeof(char));
      ++tu->numStacks;
    }

    assert(tu->stacks[tu->currentStack].top == (FIRST_STACK_SIZE << tu->currentStack));
  }

  /* The chunk fits into the last stack. */

  TU_STACK* pstack = &tu->stacks[tu->currentStack];
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

  return TU_OKAY;
}

TU_ERROR _TUfreeStack(TU* tu, void** ptr)
{
  assert(tu);
  assert(ptr);
  assert(*ptr);

  TU_STACK* stack = &tu->stacks[tu->currentStack];
  TUdbgMsg(0, "TUfreeStack called for pointer %p. Last stack is %d.\n", *ptr, tu->currentStack);
  size_t size = *((size_t*) &stack->memory[stack->top]);

#if defined(DEBUG_STACK)
  printf("TUfreeStack() called for %ld bytes (size stored at %p).\n", size,
    &stack->memory[stack->top]);
  fflush(stdout);
#endif /* DEBUG_STACK */

  assert(size < (FIRST_STACK_SIZE << tu->numStacks));

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
      "Wrong order of TUfreeStack(Array) detected. Top chunk on stack has size %ld!\n",
      size);
    fflush(stderr);
  }
#endif /* !NDEBUG */

  stack->top += size + sizeof(void*);
#if !defined(NDEBUG)
  stack->top += sizeof(int);
#endif /* !NDEBUG */

  while (stack->top == (FIRST_STACK_SIZE << tu->currentStack) && tu->currentStack > 0)
  {
    --tu->currentStack;
    stack = &tu->stacks[tu->currentStack];
  }
  *ptr = NULL;

  return TU_OKAY;
}

#if !defined(NDEBUG)

void TUassertStackConsistency(
  TU* tu
)
{
  assert(tu);

  for (int s = 0; s <= tu->currentStack; ++s)
  {
    TU_STACK* stack = &tu->stacks[s];

    void* ptr = &stack->memory[stack->top];
    TUdbgMsg(2, "Stack %d of size %d has memory range [%p,%p). top is %p\n", s, STACK_SIZE(s), stack->memory,
      stack->memory + STACK_SIZE(s), ptr);
    while (ptr < (void*)stack->memory + STACK_SIZE(s))
    {
      TUdbgMsg(4, "pointer is %p.", ptr);
      size_t size = *((size_t*) ptr);
      TUdbgMsg(0, " It indicates a chunk of size %d.\n", size);
      ptr += sizeof(size_t*);
      assert(*((int*)ptr) == PROTECTION);
      ptr += size + sizeof(int);
    }
  }
}

#endif /* !NDEBUG */

#endif /* else REPLACE_STACK_BY_MALLOC */

char* TUconsistencyMessage(const char* format, ...)
{
  assert(format);

  char buffer[256];
  va_list argptr;
  va_start(argptr, format);
  vsprintf(buffer, format, argptr);
  va_end(argptr);

  return strdup(buffer);
}
