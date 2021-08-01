#ifndef CMR_HASHTABLE_INTERNAL_H
#define CMR_HASHTABLE_INTERNAL_H

#include "env_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _CMR_HASHTABLE CMR_HASHTABLE;
typedef size_t CMR_HASHTABLE_ENTRY;
typedef size_t CMR_HASHTABLE_HASH;

CMR_ERROR CMRhashtableCreate(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_HASHTABLE** phashtable,  /**< Pointer for storing the hash table. */
  size_t initialSize,         /**< Initial size of hash table. */
  size_t initialKeyMemory     /**< Initial memory for keys. */  
);

CMR_ERROR CMRhashtableFree(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_HASHTABLE** phashtable /**< Pointer to the hash table. */
);

bool CMRhashtableFind(
  CMR_HASHTABLE* hashtable,    /**< Hash table. */
  const void* keyArray,             /**< First byte of key array. */
  size_t keyLength,           /**< Length of key array in bytes. */
  CMR_HASHTABLE_ENTRY* pentry, /**< Pointer for storing the entry in the hash table. */
  CMR_HASHTABLE_HASH* phash    /**< Pointer for storing the entry in the hash table. */
);

const void* CMRhashtableKey(
  CMR_HASHTABLE* hashtable,  /**< Hash table. */
  CMR_HASHTABLE_ENTRY entry, /**< Entry. */
  size_t* pKeyLength        /**< Length of key array. */
);

const void* CMRhashtableValue(
  CMR_HASHTABLE* hashtable,  /**< Hash table. */
  CMR_HASHTABLE_ENTRY entry  /**< Entry. */
);

CMR_ERROR CMRhashtableInsertEntryHash(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_HASHTABLE* hashtable,  /**< Hash table. */
  const void* keyArray,           /**< First byte of key array. */
  size_t keyLength,         /**< Length of key array. */
  CMR_HASHTABLE_ENTRY entry, /**< Known entry of key, determined by \ref CMRhashtableFind. */
  CMR_HASHTABLE_HASH hash,   /**< Known hash of key, determined by \ref CMRhashtableFind. */
  const void* value               /**< Value to be set. */
);

CMR_ERROR CMRhashtableInsert(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_HASHTABLE* hashtable,  /**< Hash table. */
  const void* keyArray,           /**< First byte of key array. */
  size_t keyLength,         /**< Length of key array. */
  const void* value               /**< Value to be set. */
);




#ifdef __cplusplus
}
#endif

#endif /* CMR_HASHTABLE_INTERNAL_H */
