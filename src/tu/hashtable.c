// #define TU_DEBUG /* Uncomment to debug the hash table. */

#include "hashtable.h"

#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"

typedef struct
{
  size_t keyIndex;          /**< \brief Position of first key byte in \ref TU_HASHTABLE::keyStorage. */
  size_t keyLength;       /**< \brief Length of key array. */
  TU_HASHTABLE_HASH hash; /**< \brief Hash value of key array. */
  const void* value;            /**< \brief Stored value. */
} TableData;

struct _TU_HASHTABLE
{
  size_t size;                /**< \brief Size of the hash table. */
  TableData* table;           /**< \brief Actual hash table. */

  unsigned char* keyStorage;  /**< \brief Storage for keys. */
  long freeKeyIndex;        /**< \brief First unused byte in \ref keyStorage. */
  size_t memKeyStorage;       /**< \brief Length of \ref keyStorage. */

  size_t numElements;         /**< \brief Number of stored key/value pairs. */
};

static inline
size_t hashEntry(TU_HASHTABLE* hashtable, TU_HASHTABLE_HASH hash)
{
  assert(hashtable);

  return hash % hashtable->size;
}

TU_ERROR TUhashtableCreate(TU* tu, TU_HASHTABLE** phashtable, size_t initialSize, size_t initialKeyMemory)
{
  assert(tu);
  assert(phashtable);
  assert(!*phashtable);
  assert(initialSize > 0);
  assert(initialKeyMemory > 0);

  TU_CALL( TUallocBlock(tu, phashtable) );
  TU_HASHTABLE* hashtable = *phashtable;

  hashtable->size = initialSize;
  hashtable->table = NULL;
  TU_CALL( TUallocBlockArray(tu, &hashtable->table, initialSize) );
  for (size_t i = 0; i < initialSize; ++i)
    hashtable->table[i].keyLength = 0;

  hashtable->freeKeyIndex = 0;
  hashtable->memKeyStorage = initialKeyMemory;
  hashtable->keyStorage = NULL;
  TU_CALL( TUallocBlockArray(tu, &hashtable->keyStorage, initialKeyMemory) );

  hashtable->numElements = 0;

  return TU_OKAY;
}

TU_ERROR TUhashtableFree(TU* tu, TU_HASHTABLE** phashtable)
{
  assert(tu);
  assert(phashtable);
  assert(*phashtable);

  TU_HASHTABLE* hashtable = *phashtable;

  TU_CALL( TUfreeBlockArray(tu, &hashtable->table) );
  TU_CALL( TUfreeBlockArray(tu, &hashtable->keyStorage) );
  TU_CALL( TUfreeBlock(tu, phashtable) );
  *phashtable = NULL;

  return TU_OKAY;
}

const void* TUhashtableKey(TU_HASHTABLE* hashtable, TU_HASHTABLE_HASH hash, size_t* pKeyLength)
{
  assert(hashtable);
  assert(hash >= 0);
  assert(pKeyLength);

  TableData* data = &hashtable->table[hashEntry(hashtable, hash)];
  *pKeyLength = data->keyLength;
  return &hashtable->keyStorage[data->keyIndex];
}

const void* TUhashtableValue(TU_HASHTABLE* hashtable, TU_HASHTABLE_HASH hash)
{
  assert(hashtable);
  assert(hash >= 0);

  return hashtable->table[hashEntry(hashtable, hash)].value;
}

bool TUhashtableFind(TU_HASHTABLE* hashtable, const void* keyArray, size_t keyLength, TU_HASHTABLE_ENTRY* pentry,
  TU_HASHTABLE_HASH* phash)
{
  assert(hashtable);
  assert(keyArray);
  assert(keyLength > 0);
  assert(pentry);
  assert(phash);

  TU_HASHTABLE_HASH hash = 5381;
  for (size_t i = 0; i < keyLength; ++i)
    hash = ((hash << 5) + hash) + ((unsigned char*)keyArray)[i];

  *phash = hash;
  TUdbgMsg(0, "TUhashtableFind computed hash %ld\n", hash);

  TU_HASHTABLE_ENTRY entry = hash % hashtable->size;
  while (true)
  {
    TUdbgMsg(2, "Checking entry %d\n", entry);
    /* If bucket is empty, then key does not exist. */
    if (!hashtable->table[entry].keyLength)
    {
      TUdbgMsg(2, "-> found empty bucket.\n");
      *pentry = entry;
      return false;
    }

    /* If bucket is used, we have to compare. */
    bool equal = keyLength == hashtable->table[entry].keyLength;
    size_t index = hashtable->table[entry].keyIndex;
    for (size_t i = 0; equal && i < keyLength; ++i)
      equal = equal && ((unsigned char*)keyArray)[i] == hashtable->keyStorage[index + i];
    if (equal)
    {
      TUdbgMsg(2, "-> found bucket %d with key.\n", entry);
      *pentry = entry;
      return true;
    }

    entry = (entry + 1) % hashtable->size;
  }
}

TU_ERROR TUhashtableInsertEntryHash(TU* tu, TU_HASHTABLE* hashtable, const void* keyArray, size_t keyLength,
  TU_HASHTABLE_ENTRY entry, TU_HASHTABLE_HASH hash, const void* value)
{
  assert(tu);
  assert(hashtable);
  assert(keyArray);
  assert(keyLength > 0);
  assert(hash >= 0);

  TableData* data = &hashtable->table[entry];
  if (data->keyLength)
  {
    /* Key exists already. */
    assert(data->keyIndex < hashtable->freeKeyIndex);
    assert(data->keyLength == keyLength);
    data->value = value;
    return TU_OKAY;
  }

  /* Enlarge key storage if necessary. */
  if (hashtable->freeKeyIndex + keyLength > hashtable->memKeyStorage)
  {
    do
    {
      hashtable->memKeyStorage *= 2;
    }
    while (hashtable->freeKeyIndex + keyLength > hashtable->memKeyStorage);
    TU_CALL( TUreallocBlockArray(tu, &hashtable->keyStorage, hashtable->memKeyStorage) );
  }

  /* Store entry by creating a copy of key in storage. */

  data->hash = hash;
  data->value = value;
  data->keyIndex = hashtable->freeKeyIndex;
  data->keyLength = keyLength;
  hashtable->freeKeyIndex += keyLength;
  const unsigned char* source = keyArray;
  unsigned char* target = &hashtable->keyStorage[data->keyIndex];
  for (size_t i = keyLength; i; --i)
    *target++ = *source++;
  hashtable->numElements++;

  if (hashtable->numElements > hashtable->size / 8)
  {
    TUdbgMsg(0, "Enlarging hash table.\n");
    
    /* We now double the size of the hash table. */

    size_t newSize = 2 * hashtable->size;
    TU_CALL( TUreallocBlockArray(tu, &hashtable->table, newSize) );
    for (size_t i = hashtable->size; i < newSize; ++i)
      hashtable->table[i].keyLength = 0;
    size_t oldSize = hashtable->size;
    hashtable->size = newSize;

    /* We re-insert each element based on its hash. */
    for (size_t i = 0; i < oldSize; ++i)
    {
      if (!hashtable->table[i].keyLength)
        continue;
      
      size_t j = hashEntry(hashtable, hashtable->table[i].hash);
      TUdbgMsg(2, "Hash %ld was at %d before and would like to be at %d.", hashtable->table[i].hash, i, j);
      while (j != i && hashtable->table[j].keyLength)
        j = (j+1) % hashtable->size;
      if (j == i)
      {
        TUdbgMsg(1, "-> next available entry is old one %d.\n", j);
      }
      else
      {
        TUdbgMsg(1, "-> next available entry is %d. Moving it there.\n", j);
        hashtable->table[j].hash = hashtable->table[i].hash;
        hashtable->table[j].keyIndex = hashtable->table[i].keyIndex;
        hashtable->table[j].keyLength = hashtable->table[i].keyLength;
        hashtable->table[j].value = hashtable->table[i].value;
        hashtable->table[i].keyLength = 0;
      }
    }
  }
  
  return TU_OKAY;
}

TU_ERROR TUhashtableInsert(TU* tu, TU_HASHTABLE* hashtable, const void* keyArray, size_t keyLength, const void* value)
{
  assert(tu);
  assert(hashtable);
  assert(keyArray);
  assert(keyLength > 0);

  TU_HASHTABLE_ENTRY entry;
  TU_HASHTABLE_HASH hash;
  TU_CALL( TUhashtableFind(hashtable, keyArray, keyLength, &entry, &hash) );
  TU_CALL( TUhashtableInsertEntryHash(tu, hashtable, keyArray, keyLength, entry, hash, value) );

  return TU_OKAY;
}



