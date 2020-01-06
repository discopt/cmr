#include "env_internal.h"

#include <assert.h>
#include <stdlib.h>

void TUcreateEnvironment(TU** tu)
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

void TUfreeEnvironment(TU ** tu)
{
  assert(tu != NULL);
  if ((*tu)->closeOutput)
    fclose((*tu)->output);
  free(*tu);
}

void _TUallocBlock(TU* tu, void** ptr, size_t size)
{
  assert(tu);
  assert(ptr);
  *ptr = malloc(size);
}

void _TUfreeBlock(TU* tu, void** ptr, size_t size)
{
  assert(tu);
  assert(ptr);
  assert(*ptr);
  free(*ptr);
  *ptr = NULL;
}

void _TUallocBlockArray(TU* tu, void** ptr, size_t size, size_t length)
{
  assert(tu);
  assert(ptr);
  *ptr = malloc(size * length);
}

void _TUfreeBlockArray(TU* tu, void** ptr)
{
  assert(tu);
  assert(ptr);
  assert(*ptr);
  free(*ptr);
  *ptr = NULL;
}



void _TUallocStack(TU* tu, void** ptr, size_t size)
{
  assert(tu);
  assert(ptr);
  *ptr = malloc(size);
}

void _TUfreeStack(TU* tu, void** ptr, size_t size)
{
  assert(tu);
  assert(ptr);
  assert(*ptr);
  free(*ptr);
  *ptr = NULL;
}

void _TUallocStackArray(TU* tu, void** ptr, size_t size, size_t length)
{
  assert(tu);
  assert(ptr);
  *ptr = malloc(size * length);
}

void _TUfreeStackArray(TU* tu, void** ptr)
{
  assert(tu);
  assert(ptr);
  assert(*ptr);
  free(*ptr);
  *ptr = NULL;
}
