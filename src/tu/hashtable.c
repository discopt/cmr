// #define TU_DEBUG /* Uncomment to debug the hash table. */

#include "hashtable.h"

#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include "env_internal.h"

typedef struct
{
  size_t keyIndex;          /**< \brief Position of first key byte in \ref TU_HASHTABLE::keyStorage. */
  size_t keyLength;       /**< \brief Length of key array. */
  TU_LINEARHASHTABLE_HASH hash; /**< \brief Hash value of key array. */
  const void* value;            /**< \brief Stored value. */
} LinearhashtableArrayBucket;

struct _TU_LINEARHASHTABLE_ARRAY
{
  size_t numBuckets;                        /**< \brief Size of the hash table. */
  LinearhashtableArrayBucket* buckets;  /**< \brief Actual hash table. */

  unsigned char* keyStorage;          /**< \brief Storage for keys. */
  long freeKeyIndex;                  /**< \brief First unused byte in \ref keyStorage. */
  size_t memKeyStorage;               /**< \brief Length of \ref keyStorage. */

  size_t numElements;                 /**< \brief Number of stored key/value pairs. */
};

static inline
size_t linearhashtableArrayHashToBucket(TU_LINEARHASHTABLE_ARRAY* hashtable, TU_LINEARHASHTABLE_HASH hash)
{
  assert(hashtable);

  return hash % hashtable->numBuckets;
}

TU_ERROR TUlinearhashtableArrayCreate(TU* tu, TU_LINEARHASHTABLE_ARRAY** phashtable, size_t initialSize, size_t initialKeyMemory)
{
  assert(tu);
  assert(phashtable);
  assert(!*phashtable);
  assert(initialSize > 0);
  assert(initialKeyMemory > 0);

  TU_CALL( TUallocBlock(tu, phashtable) );
  TU_LINEARHASHTABLE_ARRAY* hashtable = *phashtable;

  hashtable->numBuckets = initialSize;
  hashtable->buckets = NULL;
  TU_CALL( TUallocBlockArray(tu, &hashtable->buckets, initialSize) );
  for (size_t i = 0; i < initialSize; ++i)
    hashtable->buckets[i].keyLength = 0;

  hashtable->freeKeyIndex = 0;
  hashtable->memKeyStorage = initialKeyMemory;
  hashtable->keyStorage = NULL;
  TU_CALL( TUallocBlockArray(tu, &hashtable->keyStorage, initialKeyMemory) );

  hashtable->numElements = 0;

  return TU_OKAY;
}

TU_ERROR TUlinearhashtableArrayFree(TU* tu, TU_LINEARHASHTABLE_ARRAY** phashtable)
{
  assert(tu);
  assert(phashtable);
  assert(*phashtable);

  TU_LINEARHASHTABLE_ARRAY* hashtable = *phashtable;

  TU_CALL( TUfreeBlockArray(tu, &hashtable->buckets) );
  TU_CALL( TUfreeBlockArray(tu, &hashtable->keyStorage) );
  TU_CALL( TUfreeBlock(tu, phashtable) );
  *phashtable = NULL;

  return TU_OKAY;
}

const void* TUlinearhashtableArrayKey(TU_LINEARHASHTABLE_ARRAY* hashtable, TU_LINEARHASHTABLE_HASH hash, size_t* pKeyLength)
{
  assert(hashtable);
  assert(hash >= 0);
  assert(pKeyLength);

  LinearhashtableArrayBucket* data = &hashtable->buckets[linearhashtableArrayHashToBucket(hashtable, hash)];
  *pKeyLength = data->keyLength;
  return &hashtable->keyStorage[data->keyIndex];
}

const void* TUlinearhashtableArrayValue(TU_LINEARHASHTABLE_ARRAY* hashtable, TU_LINEARHASHTABLE_HASH hash)
{
  assert(hashtable);
  assert(hash >= 0);

  return hashtable->buckets[linearhashtableArrayHashToBucket(hashtable, hash)].value;
}

bool TUlinearhashtableArrayFind(TU_LINEARHASHTABLE_ARRAY* hashtable, const void* keyArray, size_t keyLength,
  TU_LINEARHASHTABLE_BUCKET* pbucket, TU_LINEARHASHTABLE_HASH* phash)
{
  assert(hashtable);
  assert(keyArray);
  assert(keyLength > 0);
  assert(pbucket);
  assert(phash);

  TU_LINEARHASHTABLE_HASH hash = 5381;
  for (size_t i = 0; i < keyLength; ++i)
    hash = ((hash << 5) + hash) + ((unsigned char*)keyArray)[i];

  *phash = hash;
  TUdbgMsg(0, "TUhashtableFind computed hash %ld\n", hash);

  TU_LINEARHASHTABLE_BUCKET bucket = hash % hashtable->numBuckets;
  while (true)
  {
    TUdbgMsg(2, "Checking bucket %d\n", bucket);
    /* If bucket is empty, then key does not exist. */
    if (!hashtable->buckets[bucket].keyLength)
    {
      TUdbgMsg(2, "-> found empty bucket.\n");
      *pbucket = bucket;
      return false;
    }

    /* If bucket is used, we have to compare. */
    bool equal = keyLength == hashtable->buckets[bucket].keyLength;
    size_t index = hashtable->buckets[bucket].keyIndex;
    for (size_t i = 0; equal && i < keyLength; ++i)
      equal = equal && ((unsigned char*)keyArray)[i] == hashtable->keyStorage[index + i];
    if (equal)
    {
      TUdbgMsg(2, "-> found bucket %d with key.\n", bucket);
      *pbucket = bucket;
      return true;
    }

    bucket = (bucket + 1) % hashtable->numBuckets;
  }
}

TU_ERROR TUlinearhashtableArrayInsertBucketHash(TU* tu, TU_LINEARHASHTABLE_ARRAY* hashtable, const void* keyArray,
  size_t keyLength, TU_LINEARHASHTABLE_BUCKET bucket, TU_LINEARHASHTABLE_HASH hash, const void* value)
{
  assert(tu);
  assert(hashtable);
  assert(keyArray);
  assert(keyLength > 0);
  assert(hash >= 0);

  LinearhashtableArrayBucket* bucketData = &hashtable->buckets[bucket];
  if (bucketData->keyLength)
  {
    /* Key exists already. */
    assert(bucketData->keyIndex < hashtable->freeKeyIndex);
    assert(bucketData->keyLength == keyLength);
    bucketData->value = value;
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

  bucketData->hash = hash;
  bucketData->value = value;
  bucketData->keyIndex = hashtable->freeKeyIndex;
  bucketData->keyLength = keyLength;
  hashtable->freeKeyIndex += keyLength;
  const unsigned char* source = keyArray;
  unsigned char* target = &hashtable->keyStorage[bucketData->keyIndex];
  for (size_t i = keyLength; i; --i)
    *target++ = *source++;
  hashtable->numElements++;

  if (hashtable->numElements > hashtable->numBuckets / 8)
  {
    TUdbgMsg(0, "Enlarging hash table.\n");
    
    /* We now double the size of the hash table. */

    size_t newSize = 2 * hashtable->numBuckets;
    TU_CALL( TUreallocBlockArray(tu, &hashtable->buckets, newSize) );
    for (size_t i = hashtable->numBuckets; i < newSize; ++i)
      hashtable->buckets[i].keyLength = 0;
    size_t oldSize = hashtable->numBuckets;
    hashtable->numBuckets = newSize;

    /* We re-insert each element based on its hash. */
    for (size_t i = 0; i < oldSize; ++i)
    {
      if (!hashtable->buckets[i].keyLength)
        continue;
      
      size_t j = linearhashtableArrayHashToBucket(hashtable, hashtable->buckets[i].hash);
      TUdbgMsg(2, "Hash %ld was at %d before and would like to be at %d.", hashtable->buckets[i].hash, i, j);
      while (j != i && hashtable->buckets[j].keyLength)
        j = (j+1) % hashtable->numBuckets;
      if (j == i)
      {
        TUdbgMsg(1, "-> next available bucket is old one %d.\n", j);
      }
      else
      {
        TUdbgMsg(1, "-> next available bucket is %d. Moving it there.\n", j);
        hashtable->buckets[j].hash = hashtable->buckets[i].hash;
        hashtable->buckets[j].keyIndex = hashtable->buckets[i].keyIndex;
        hashtable->buckets[j].keyLength = hashtable->buckets[i].keyLength;
        hashtable->buckets[j].value = hashtable->buckets[i].value;
        hashtable->buckets[i].keyLength = 0;
      }
    }
  }
  
  return TU_OKAY;
}

TU_ERROR TUlinearhashtableArrayInsert(TU* tu, TU_LINEARHASHTABLE_ARRAY* hashtable, const void* keyArray, size_t keyLength, const void* value)
{
  assert(tu);
  assert(hashtable);
  assert(keyArray);
  assert(keyLength > 0);

  TU_LINEARHASHTABLE_BUCKET bucket;
  TU_LINEARHASHTABLE_HASH hash;
  TU_CALL( TUlinearhashtableArrayFind(hashtable, keyArray, keyLength, &bucket, &hash) );
  TU_CALL( TUlinearhashtableArrayInsertBucketHash(tu, hashtable, keyArray, keyLength, bucket, hash, value) );

  return TU_OKAY;
}


typedef struct
{
  TU_LISTHASHTABLE_HASH hash;
  TU_LISTHASHTABLE_ENTRY next;
  TU_LISTHASHTABLE_VALUE value;
} ListhashtableNode;

typedef ListhashtableNode* _TU_LISTHASHTABLE_NODE;

struct _TU_LISTHASHTABLE
{
  size_t numBuckets;                /**< \brief Number of buckets. */
  TU_LISTHASHTABLE_ENTRY* buckets;  /**< \brief Array with buckets. */

  size_t memNodes;                  /**< \brief Memory allocated for nodes. */
  ListhashtableNode* nodes;         /**< \brief Array with entries. */
  TU_LISTHASHTABLE_ENTRY firstFree; /**< \brief Start of free list. */
};


TU_ERROR TUlisthashtableCreate(TU* tu, TU_LISTHASHTABLE** phashtable, size_t initialNumBuckets, size_t initialMemNodes)
{
  assert(tu);
  assert(phashtable);
  assert(initialNumBuckets > 0);
  assert(initialMemNodes > 0);

  TU_CALL( TUallocBlock(tu, phashtable) );
  TU_LISTHASHTABLE* hashtable = *phashtable;
  assert(hashtable);
  
  TUdbgMsg(6, "Creating listhashtable with %d buckets and memory for %d nodes.\n", initialNumBuckets, initialMemNodes);

  hashtable->numBuckets = initialNumBuckets;
  hashtable->buckets = NULL;
  TU_CALL( TUallocBlockArray(tu, &hashtable->buckets, initialNumBuckets) );
  for (size_t i = 0; i < initialNumBuckets; ++i)
    hashtable->buckets[i] = SIZE_MAX;

  hashtable->memNodes = initialMemNodes;
  hashtable->nodes = NULL;
  TU_CALL( TUallocBlockArray( tu, &hashtable->nodes, initialMemNodes) );
  for (size_t i = 0; i < initialMemNodes-1; ++i)
    hashtable->nodes[i].next = i+1;
  hashtable->nodes[initialMemNodes-1].next = SIZE_MAX;
  hashtable->firstFree = 0;

  return TU_OKAY;
}

TU_ERROR TUlisthashtableFree(TU* tu, TU_LISTHASHTABLE** phashtable)
{
  assert(tu);
  assert(phashtable);

  TU_LISTHASHTABLE* hashtable = *phashtable;
  if (!hashtable)
    return TU_OKAY;

  TU_CALL( TUfreeBlockArray(tu, &hashtable->nodes) );
  TU_CALL( TUfreeBlockArray(tu, &hashtable->buckets) );
  TU_CALL( TUfreeBlock(tu, phashtable) );

  return TU_OKAY;
}

static inline
TU_LISTHASHTABLE_BUCKET listhashtableHashToBucket(TU_LISTHASHTABLE* hashtable, TU_LISTHASHTABLE_HASH hash)
{
  assert(hashtable);

  return hash % hashtable->numBuckets;
}

TU_LISTHASHTABLE_ENTRY TUlisthashtableFindFirst(TU_LISTHASHTABLE* hashtable, TU_LISTHASHTABLE_HASH hash)
{
  assert(hashtable);

  TU_LISTHASHTABLE_BUCKET bucket = listhashtableHashToBucket(hashtable, hash);
  TU_LISTHASHTABLE_ENTRY entry = hashtable->buckets[bucket];
  if (entry == SIZE_MAX)
    return SIZE_MAX;

  ListhashtableNode* node = &hashtable->nodes[entry];
  if (node->hash == hash)
    return entry;

  return TUlisthashtableFindNext(hashtable, hash, entry);
}

TU_LISTHASHTABLE_ENTRY TUlisthashtableFindNext(TU_LISTHASHTABLE* hashtable, TU_LISTHASHTABLE_HASH hash,
  TU_LISTHASHTABLE_ENTRY entry)
{
  assert(hashtable);
  assert(entry < hashtable->memNodes);

  while (true)
  {
    entry = hashtable->nodes[entry].next;
    if (entry == SIZE_MAX)
      return SIZE_MAX;

    if (hashtable->nodes[entry].hash == hash)
      return entry;
  }
}

TU_LISTHASHTABLE_VALUE TUlisthashtableValue(TU_LISTHASHTABLE* hashtable, TU_LISTHASHTABLE_ENTRY entry)
{
  assert(hashtable);

  return hashtable->nodes[entry].value;
}

TU_LISTHASHTABLE_HASH TUlisthashtableHash(TU_LISTHASHTABLE* hashtable, TU_LISTHASHTABLE_ENTRY entry)
{
  assert(hashtable);

  return hashtable->nodes[entry].hash;
}

TU_ERROR TUlisthashtableInsert(TU* tu, TU_LISTHASHTABLE* hashtable, TU_LISTHASHTABLE_HASH hash,
  TU_LISTHASHTABLE_VALUE value, TU_LISTHASHTABLE_ENTRY* pentry)
{
  assert(tu);
  assert(hashtable);

  /* If necessary, reallocate nodes */
  if (hashtable->firstFree == SIZE_MAX)
  {
    size_t newMemNodes = 2 * hashtable->memNodes;
    TU_CALL( TUreallocBlockArray(tu, &hashtable->nodes, newMemNodes) );
    for (size_t i = hashtable->memNodes; i+1 < newMemNodes; ++i)
      hashtable->nodes[i].next = i+1;
    hashtable->nodes[newMemNodes-1].next = SIZE_MAX;
    hashtable->firstFree = hashtable->memNodes;
    hashtable->memNodes = newMemNodes;
  }

  ListhashtableNode* newNode = &hashtable->nodes[hashtable->firstFree];
  TU_LISTHASHTABLE_ENTRY next = newNode->next;
  newNode->hash = hash;
  newNode->value = value;
  TU_LISTHASHTABLE_BUCKET bucket = listhashtableHashToBucket(hashtable, hash);
  newNode->next = hashtable->buckets[bucket];
  if (pentry)
    *pentry = hashtable->firstFree;
  hashtable->buckets[bucket] = hashtable->firstFree;
  hashtable->firstFree = next;

  return TU_OKAY;
}

TU_ERROR TUlisthashtableRemove(TU* tu, TU_LISTHASHTABLE* hashtable, TU_LISTHASHTABLE_ENTRY entry)
{
  assert(tu);
  assert(hashtable);

  TU_LISTHASHTABLE_HASH hash = hashtable->nodes[entry].hash;
  TU_LISTHASHTABLE_ENTRY next = hashtable->nodes[entry].next;
  TU_LISTHASHTABLE_BUCKET bucket = listhashtableHashToBucket(hashtable, hash);
  
  TU_LISTHASHTABLE_ENTRY previous = hashtable->buckets[bucket];
  if (previous == entry)
  {
    hashtable->buckets[bucket] = next;
  }
  else
  {
    TU_LISTHASHTABLE_ENTRY current = hashtable->nodes[previous].next;
    while (current != entry)
    {
      if (current == SIZE_MAX)
        return TU_ERROR_INVALID;

      previous = current;
      current = hashtable->nodes[current].next;
    }
    hashtable->nodes[previous].next = next;
  }

  hashtable->nodes[entry].next = hashtable->firstFree;
  hashtable->firstFree = entry;

  return TU_OKAY;
}
