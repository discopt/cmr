// #define CMR_DEBUG /* Uncomment to debug the hash table. */

#include "hashtable.h"

#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>

#include "env_internal.h"

typedef struct
{
  size_t keyIndex;          /**< \brief Position of first key byte in \ref CMR_HASHTABLE::keyStorage. */
  size_t keyLength;       /**< \brief Length of key array. */
  CMR_HASHTABLE_HASH hash; /**< \brief Hash value of key array. */
  const void* value;            /**< \brief Stored value. */
} TableData;

struct _CMR_HASHTABLE
{
  size_t size;                /**< \brief Size of the hash table. */
  TableData* table;           /**< \brief Actual hash table. */

  unsigned char* keyStorage;  /**< \brief Storage for keys. */
  long freeKeyIndex;        /**< \brief First unused byte in \ref keyStorage. */
  size_t memKeyStorage;       /**< \brief Length of \ref keyStorage. */

  size_t numElements;         /**< \brief Number of stored key/value pairs. */
};

static inline
size_t hashEntry(CMR_HASHTABLE* hashtable, CMR_HASHTABLE_HASH hash)
{
  assert(hashtable);

  return hash % hashtable->size;
}

CMR_ERROR CMRhashtableCreate(CMR* cmr, CMR_HASHTABLE** phashtable, size_t initialSize, size_t initialKeyMemory)
{
  assert(cmr);
  assert(phashtable);
  assert(!*phashtable);
  assert(initialSize > 0);
  assert(initialKeyMemory > 0);

  CMR_CALL( CMRallocBlock(cmr, phashtable) );
  CMR_HASHTABLE* hashtable = *phashtable;

  hashtable->size = initialSize;
  hashtable->table = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &hashtable->table, initialSize) );
  for (size_t i = 0; i < initialSize; ++i)
    hashtable->table[i].keyLength = 0;

  hashtable->freeKeyIndex = 0;
  hashtable->memKeyStorage = initialKeyMemory;
  hashtable->keyStorage = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &hashtable->keyStorage, initialKeyMemory) );

  hashtable->numElements = 0;

  return CMR_OKAY;
}

CMR_ERROR CMRhashtableFree(CMR* cmr, CMR_HASHTABLE** phashtable)
{
  assert(cmr);
  assert(phashtable);
  assert(*phashtable);

  CMR_HASHTABLE* hashtable = *phashtable;

  CMR_CALL( CMRfreeBlockArray(cmr, &hashtable->table) );
  CMR_CALL( CMRfreeBlockArray(cmr, &hashtable->keyStorage) );
  CMR_CALL( CMRfreeBlock(cmr, phashtable) );
  *phashtable = NULL;

  return CMR_OKAY;
}

const void* CMRhashtableKey(CMR_HASHTABLE* hashtable, CMR_HASHTABLE_HASH hash, size_t* pKeyLength)
{
  assert(hashtable);
  assert(hash >= 0);
  assert(pKeyLength);

  TableData* data = &hashtable->table[hashEntry(hashtable, hash)];
  *pKeyLength = data->keyLength;
  return &hashtable->keyStorage[data->keyIndex];
}

const void* CMRhashtableValue(CMR_HASHTABLE* hashtable, CMR_HASHTABLE_HASH hash)
{
  assert(hashtable);
  assert(hash >= 0);

  return hashtable->table[hashEntry(hashtable, hash)].value;
}

bool CMRhashtableFind(CMR_HASHTABLE* hashtable, const void* keyArray, size_t keyLength, CMR_HASHTABLE_ENTRY* pentry,
  CMR_HASHTABLE_HASH* phash)
{
  assert(hashtable);
  assert(keyArray);
  assert(keyLength > 0);
  assert(pentry);
  assert(phash);

  CMR_HASHTABLE_HASH hash = 5381;
  for (size_t i = 0; i < keyLength; ++i)
    hash = ((hash << 5) + hash) + ((unsigned char*)keyArray)[i];

  *phash = hash;
  CMRdbgMsg(0, "CMRhashtableFind computed hash %ld\n", hash);

  CMR_HASHTABLE_ENTRY entry = hash % hashtable->size;
  while (true)
  {
    CMRdbgMsg(2, "Checking entry %d\n", entry);
    /* If bucket is empty, then key does not exist. */
    if (!hashtable->table[entry].keyLength)
    {
      CMRdbgMsg(2, "-> found empty bucket.\n");
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
      CMRdbgMsg(2, "-> found bucket %d with key.\n", entry);
      *pentry = entry;
      return true;
    }

    entry = (entry + 1) % hashtable->size;
  }
}

CMR_ERROR CMRhashtableInsertEntryHash(CMR* cmr, CMR_HASHTABLE* hashtable, const void* keyArray, size_t keyLength,
  CMR_HASHTABLE_ENTRY entry, CMR_HASHTABLE_HASH hash, const void* value)
{
  assert(cmr);
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
    return CMR_OKAY;
  }

  /* Enlarge key storage if necessary. */
  if (hashtable->freeKeyIndex + keyLength > hashtable->memKeyStorage)
  {
    do
    {
      hashtable->memKeyStorage *= 2;
    }
    while (hashtable->freeKeyIndex + keyLength > hashtable->memKeyStorage);
    CMR_CALL( CMRreallocBlockArray(cmr, &hashtable->keyStorage, hashtable->memKeyStorage) );
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
    CMRdbgMsg(0, "Enlarging hash table.\n");
    
    /* We now double the size of the hash table. */

    size_t newSize = 2 * hashtable->size;
    CMR_CALL( CMRreallocBlockArray(cmr, &hashtable->table, newSize) );
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
      CMRdbgMsg(2, "Hash %ld was at %d before and would like to be at %d.", hashtable->table[i].hash, i, j);
      while (j != i && hashtable->table[j].keyLength)
        j = (j+1) % hashtable->size;
      if (j == i)
      {
        CMRdbgMsg(1, "-> next available entry is old one %d.\n", j);
      }
      else
      {
        CMRdbgMsg(1, "-> next available entry is %d. Moving it there.\n", j);
        hashtable->table[j].hash = hashtable->table[i].hash;
        hashtable->table[j].keyIndex = hashtable->table[i].keyIndex;
        hashtable->table[j].keyLength = hashtable->table[i].keyLength;
        hashtable->table[j].value = hashtable->table[i].value;
        hashtable->table[i].keyLength = 0;
      }
    }
  }
  
  return CMR_OKAY;
}

CMR_ERROR CMRhashtableInsert(CMR* cmr, CMR_HASHTABLE* hashtable, const void* keyArray, size_t keyLength, const void* value)
{
  assert(cmr);
  assert(hashtable);
  assert(keyArray);
  assert(keyLength > 0);

  CMR_HASHTABLE_ENTRY entry;
  CMR_HASHTABLE_HASH hash;
  CMR_CALL( CMRhashtableFind(hashtable, keyArray, keyLength, &entry, &hash) );
  CMR_CALL( CMRhashtableInsertEntryHash(cmr, hashtable, keyArray, keyLength, entry, hash, value) );

  return CMR_OKAY;
}



