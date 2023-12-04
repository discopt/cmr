// #define CMR_DEBUG /* Uncomment to debug the hash table. */

#include "hashtable.h"

#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>

#include "env_internal.h"

typedef struct
{
  size_t keyIndex;                /**< \brief Position of first key byte in \ref CMR_HASHTABLE::keyStorage. */
  size_t keyLength;               /**< \brief Length of key array. */
  CMR_LINEARHASHTABLE_HASH hash;  /**< \brief Hash value of key array. */
  const void* value;              /**< \brief Stored value. */
} LinearhashtableArrayBucket;

struct _CMR_LINEARHASHTABLE_ARRAY
{
  size_t numBuckets;                    /**< \brief Size of the hash table. */
  LinearhashtableArrayBucket* buckets;  /**< \brief Actual hash table. */

  unsigned char* keyStorage;            /**< \brief Storage for keys. */
  size_t freeKeyIndex;                  /**< \brief First unused byte in \ref keyStorage. */
  size_t memKeyStorage;                 /**< \brief Length of \ref keyStorage. */

  size_t numElements;                   /**< \brief Number of stored key/value pairs. */
};

static inline
size_t linearhashtableArrayHashToBucket(CMR_LINEARHASHTABLE_ARRAY* hashtable, CMR_LINEARHASHTABLE_HASH hash)
{
  assert(hashtable);

  return hash % hashtable->numBuckets;
}

CMR_ERROR CMRlinearhashtableArrayCreate(CMR* cmr, CMR_LINEARHASHTABLE_ARRAY** phashtable, size_t initialSize, size_t initialKeyMemory)
{
  assert(cmr);
  assert(phashtable);
  assert(!*phashtable);
  assert(initialSize > 0);
  assert(initialKeyMemory > 0);

  CMR_CALL( CMRallocBlock(cmr, phashtable) );
  CMR_LINEARHASHTABLE_ARRAY* hashtable = *phashtable;

  hashtable->numBuckets = initialSize;
  hashtable->buckets = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &hashtable->buckets, initialSize) );
  for (size_t i = 0; i < initialSize; ++i)
    hashtable->buckets[i].keyLength = 0;

  hashtable->freeKeyIndex = 0;
  hashtable->memKeyStorage = initialKeyMemory;
  hashtable->keyStorage = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &hashtable->keyStorage, initialKeyMemory) );

  hashtable->numElements = 0;

  return CMR_OKAY;
}

CMR_ERROR CMRlinearhashtableArrayFree(CMR* cmr, CMR_LINEARHASHTABLE_ARRAY** phashtable)
{
  assert(cmr);
  assert(phashtable);
  assert(*phashtable);

  CMR_LINEARHASHTABLE_ARRAY* hashtable = *phashtable;

  CMR_CALL( CMRfreeBlockArray(cmr, &hashtable->buckets) );
  CMR_CALL( CMRfreeBlockArray(cmr, &hashtable->keyStorage) );
  CMR_CALL( CMRfreeBlock(cmr, phashtable) );
  *phashtable = NULL;

  return CMR_OKAY;
}

const void* CMRlinearhashtableArrayKey(CMR_LINEARHASHTABLE_ARRAY* hashtable, CMR_LINEARHASHTABLE_HASH hash, size_t* pKeyLength)
{
  assert(hashtable);
  assert(pKeyLength);

  LinearhashtableArrayBucket* data = &hashtable->buckets[linearhashtableArrayHashToBucket(hashtable, hash)];
  *pKeyLength = data->keyLength;
  return &hashtable->keyStorage[data->keyIndex];
}

const void* CMRlinearhashtableArrayValue(CMR_LINEARHASHTABLE_ARRAY* hashtable, CMR_LINEARHASHTABLE_HASH hash)
{
  assert(hashtable);

  return hashtable->buckets[linearhashtableArrayHashToBucket(hashtable, hash)].value;
}

bool CMRlinearhashtableArrayFind(CMR_LINEARHASHTABLE_ARRAY* hashtable, const void* keyArray, size_t keyLength,
  CMR_LINEARHASHTABLE_BUCKET* pbucket, CMR_LINEARHASHTABLE_HASH* phash)
{
  assert(hashtable);
  assert(keyArray);
  assert(keyLength > 0);
  assert(pbucket);
  assert(phash);

  CMR_LINEARHASHTABLE_HASH hash = 5381;
  for (size_t i = 0; i < keyLength; ++i)
    hash = ((hash << 5) + hash) + ((unsigned char*)keyArray)[i];

  *phash = hash;
  CMRdbgMsg(0, "CMRhashtableFind computed hash %ld\n", hash);

  CMR_LINEARHASHTABLE_BUCKET bucket = hash % hashtable->numBuckets;
  while (true)
  {
    CMRdbgMsg(2, "Checking bucket %d\n", bucket);
    /* If bucket is empty, then key does not exist. */
    if (!hashtable->buckets[bucket].keyLength)
    {
      CMRdbgMsg(2, "-> found empty bucket.\n");
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
      CMRdbgMsg(2, "-> found bucket %d with key.\n", bucket);
      *pbucket = bucket;
      return true;
    }

    bucket = (bucket + 1) % hashtable->numBuckets;
  }
}

CMR_ERROR CMRlinearhashtableArrayInsertBucketHash(CMR* cmr, CMR_LINEARHASHTABLE_ARRAY* hashtable, const void* keyArray,
  size_t keyLength, CMR_LINEARHASHTABLE_BUCKET bucket, CMR_LINEARHASHTABLE_HASH hash, const void* value)
{
  assert(cmr);
  assert(hashtable);
  assert(keyArray);

  LinearhashtableArrayBucket* bucketData = &hashtable->buckets[bucket];
  if (bucketData->keyLength)
  {
    /* Key exists already. */
    assert(bucketData->keyIndex < hashtable->freeKeyIndex);
    assert(bucketData->keyLength == keyLength);
    bucketData->value = value;
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
    CMRdbgMsg(0, "Enlarging hash table.\n");
    
    /* We now double the size of the hash table. */

    size_t newSize = 2 * hashtable->numBuckets;
    CMR_CALL( CMRreallocBlockArray(cmr, &hashtable->buckets, newSize) );
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
      CMRdbgMsg(2, "Hash %ld was at %d before and would like to be at %d.", hashtable->buckets[i].hash, i, j);
      while (j != i && hashtable->buckets[j].keyLength)
        j = (j+1) % hashtable->numBuckets;
      if (j == i)
      {
        CMRdbgMsg(1, "-> next available bucket is old one %d.\n", j);
      }
      else
      {
        CMRdbgMsg(1, "-> next available bucket is %d. Moving it there.\n", j);
        hashtable->buckets[j].hash = hashtable->buckets[i].hash;
        hashtable->buckets[j].keyIndex = hashtable->buckets[i].keyIndex;
        hashtable->buckets[j].keyLength = hashtable->buckets[i].keyLength;
        hashtable->buckets[j].value = hashtable->buckets[i].value;
        hashtable->buckets[i].keyLength = 0;
      }
    }
  }
  
  return CMR_OKAY;
}

CMR_ERROR CMRlinearhashtableArrayInsert(CMR* cmr, CMR_LINEARHASHTABLE_ARRAY* hashtable, const void* keyArray, size_t keyLength, const void* value)
{
  assert(cmr);
  assert(hashtable);
  assert(keyArray);
  assert(keyLength > 0);

  CMR_LINEARHASHTABLE_BUCKET bucket;
  CMR_LINEARHASHTABLE_HASH hash;
  CMR_CALL( CMRlinearhashtableArrayFind(hashtable, keyArray, keyLength, &bucket, &hash) );
  CMR_CALL( CMRlinearhashtableArrayInsertBucketHash(cmr, hashtable, keyArray, keyLength, bucket, hash, value) );

  return CMR_OKAY;
}


typedef struct
{
  CMR_LISTHASHTABLE_HASH hash;
  CMR_LISTHASHTABLE_ENTRY next;
  CMR_LISTHASHTABLE_VALUE value;
} ListhashtableNode;

typedef ListhashtableNode* _CMR_LISTHASHTABLE_NODE;

struct _CMR_LISTHASHTABLE
{
  size_t numBuckets;                /**< \brief Number of buckets. */
  CMR_LISTHASHTABLE_ENTRY* buckets;  /**< \brief Array with buckets. */

  size_t memNodes;                  /**< \brief Memory allocated for nodes. */
  ListhashtableNode* nodes;         /**< \brief Array with entries. */
  CMR_LISTHASHTABLE_ENTRY firstFree; /**< \brief Start of free list. */
};


CMR_ERROR CMRlisthashtableCreate(CMR* cmr, CMR_LISTHASHTABLE** phashtable, size_t initialNumBuckets, size_t initialMemNodes)
{
  assert(cmr);
  assert(phashtable);
  assert(initialNumBuckets > 0);
  assert(initialMemNodes > 0);

  CMR_CALL( CMRallocBlock(cmr, phashtable) );
  CMR_LISTHASHTABLE* hashtable = *phashtable;
  assert(hashtable);
  
  CMRdbgMsg(6, "Creating listhashtable with %d buckets and memory for %d nodes.\n", initialNumBuckets, initialMemNodes);

  hashtable->numBuckets = initialNumBuckets;
  hashtable->buckets = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &hashtable->buckets, initialNumBuckets) );
  for (size_t i = 0; i < initialNumBuckets; ++i)
    hashtable->buckets[i] = SIZE_MAX;

  hashtable->memNodes = initialMemNodes;
  hashtable->nodes = NULL;
  CMR_CALL( CMRallocBlockArray(cmr, &hashtable->nodes, initialMemNodes) );
  for (size_t i = 0; i < initialMemNodes-1; ++i)
    hashtable->nodes[i].next = i+1;
  hashtable->nodes[initialMemNodes-1].next = SIZE_MAX;
  hashtable->firstFree = 0;

  return CMR_OKAY;
}

CMR_ERROR CMRlisthashtableFree(CMR* cmr, CMR_LISTHASHTABLE** phashtable)
{
  assert(cmr);
  assert(phashtable);

  CMR_LISTHASHTABLE* hashtable = *phashtable;
  if (!hashtable)
    return CMR_OKAY;

  CMR_CALL( CMRfreeBlockArray(cmr, &hashtable->nodes) );
  CMR_CALL( CMRfreeBlockArray(cmr, &hashtable->buckets) );
  CMR_CALL( CMRfreeBlock(cmr, phashtable) );

  return CMR_OKAY;
}

static inline
CMR_LISTHASHTABLE_BUCKET listhashtableHashToBucket(CMR_LISTHASHTABLE* hashtable, CMR_LISTHASHTABLE_HASH hash)
{
  assert(hashtable);

  return hash % hashtable->numBuckets;
}

CMR_LISTHASHTABLE_ENTRY CMRlisthashtableFindFirst(CMR_LISTHASHTABLE* hashtable, CMR_LISTHASHTABLE_HASH hash)
{
  assert(hashtable);

  CMR_LISTHASHTABLE_BUCKET bucket = listhashtableHashToBucket(hashtable, hash);
  CMR_LISTHASHTABLE_ENTRY entry = hashtable->buckets[bucket];
  if (entry == SIZE_MAX)
    return SIZE_MAX;

  ListhashtableNode* node = &hashtable->nodes[entry];
  if (node->hash == hash)
    return entry;

  return CMRlisthashtableFindNext(hashtable, hash, entry);
}

CMR_LISTHASHTABLE_ENTRY CMRlisthashtableFindNext(CMR_LISTHASHTABLE* hashtable, CMR_LISTHASHTABLE_HASH hash,
  CMR_LISTHASHTABLE_ENTRY entry)
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

CMR_LISTHASHTABLE_VALUE CMRlisthashtableValue(CMR_LISTHASHTABLE* hashtable, CMR_LISTHASHTABLE_ENTRY entry)
{
  assert(hashtable);

  return hashtable->nodes[entry].value;
}

CMR_LISTHASHTABLE_HASH CMRlisthashtableHash(CMR_LISTHASHTABLE* hashtable, CMR_LISTHASHTABLE_ENTRY entry)
{
  assert(hashtable);

  return hashtable->nodes[entry].hash;
}

size_t CMRlisthashtableNumBuckets(CMR_LISTHASHTABLE* hashtable)
{
  return hashtable->numBuckets;
}

CMR_ERROR CMRlisthashtableInsert(CMR* cmr, CMR_LISTHASHTABLE* hashtable, CMR_LISTHASHTABLE_HASH hash,
  CMR_LISTHASHTABLE_VALUE value, CMR_LISTHASHTABLE_ENTRY* pentry)
{
  assert(cmr);
  assert(hashtable);

  /* If necessary, reallocate nodes */
  if (hashtable->firstFree == SIZE_MAX)
  {
    size_t newMemNodes = 2 * hashtable->memNodes;
    CMR_CALL( CMRreallocBlockArray(cmr, &hashtable->nodes, newMemNodes) );
    for (size_t i = hashtable->memNodes; i+1 < newMemNodes; ++i)
      hashtable->nodes[i].next = i+1;
    hashtable->nodes[newMemNodes-1].next = SIZE_MAX;
    hashtable->firstFree = hashtable->memNodes;
    hashtable->memNodes = newMemNodes;
  }

  ListhashtableNode* newNode = &hashtable->nodes[hashtable->firstFree];
  CMR_LISTHASHTABLE_ENTRY next = newNode->next;
  newNode->hash = hash;
  newNode->value = value;
  CMR_LISTHASHTABLE_BUCKET bucket = listhashtableHashToBucket(hashtable, hash);
  newNode->next = hashtable->buckets[bucket];
  if (pentry)
    *pentry = hashtable->firstFree;
  hashtable->buckets[bucket] = hashtable->firstFree;
  hashtable->firstFree = next;

  return CMR_OKAY;
}

CMR_ERROR CMRlisthashtableRemove(CMR* cmr, CMR_LISTHASHTABLE* hashtable, CMR_LISTHASHTABLE_ENTRY entry)
{
  CMR_UNUSED(cmr);

  assert(cmr);
  assert(hashtable);

  CMR_LISTHASHTABLE_HASH hash = hashtable->nodes[entry].hash;
  CMR_LISTHASHTABLE_ENTRY next = hashtable->nodes[entry].next;
  CMR_LISTHASHTABLE_BUCKET bucket = listhashtableHashToBucket(hashtable, hash);
  
  CMR_LISTHASHTABLE_ENTRY previous = hashtable->buckets[bucket];
  if (previous == entry)
  {
    hashtable->buckets[bucket] = next;
  }
  else
  {
    CMR_LISTHASHTABLE_ENTRY current = hashtable->nodes[previous].next;
    while (current != entry)
    {
      if (current == SIZE_MAX)
        return CMR_ERROR_INVALID;

      previous = current;
      current = hashtable->nodes[current].next;
    }
    hashtable->nodes[previous].next = next;
  }

  hashtable->nodes[entry].next = hashtable->firstFree;
  hashtable->firstFree = entry;

  return CMR_OKAY;
}
