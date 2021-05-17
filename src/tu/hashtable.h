#ifndef TU_HASHTABLE_INTERNAL_H
#define TU_HASHTABLE_INTERNAL_H

#include "env_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _TU_HASHTABLE TU_HASHTABLE;
typedef size_t TU_HASHTABLE_ENTRY;
typedef size_t TU_HASHTABLE_HASH;

TU_ERROR TUhashtableCreate(
  TU* tu,                     /**< \ref TU environment. */
  TU_HASHTABLE** phashtable,  /**< Pointer for storing the hash table. */
  size_t initialSize,         /**< Initial size of hash table. */
  size_t initialKeyMemory     /**< Initial memory for keys. */  
);

TU_ERROR TUhashtableFree(
  TU* tu,                   /**< \ref TU environment. */
  TU_HASHTABLE** phashtable /**< Pointer to the hash table. */
);

bool TUhashtableFind(
  TU_HASHTABLE* hashtable,    /**< Hash table. */
  const void* keyArray,             /**< First byte of key array. */
  size_t keyLength,           /**< Length of key array in bytes. */
  TU_HASHTABLE_ENTRY* pentry, /**< Pointer for storing the entry in the hash table. */
  TU_HASHTABLE_HASH* phash    /**< Pointer for storing the entry in the hash table. */
);

const void* TUhashtableKey(
  TU_HASHTABLE* hashtable,  /**< Hash table. */
  TU_HASHTABLE_ENTRY entry, /**< Entry. */
  size_t* pKeyLength        /**< Length of key array. */
);

const void* TUhashtableValue(
  TU_HASHTABLE* hashtable,  /**< Hash table. */
  TU_HASHTABLE_ENTRY entry  /**< Entry. */
);

TU_ERROR TUhashtableInsertEntryHash(
  TU* tu,                   /**< \ref TU environment. */
  TU_HASHTABLE* hashtable,  /**< Hash table. */
  const void* keyArray,           /**< First byte of key array. */
  size_t keyLength,         /**< Length of key array. */
  TU_HASHTABLE_ENTRY entry, /**< Known entry of key, determined by \ref TUhashtableFind. */
  TU_HASHTABLE_HASH hash,   /**< Known hash of key, determined by \ref TUhashtableFind. */
  const void* value               /**< Value to be set. */
);

TU_ERROR TUhashtableInsert(
  TU* tu,                   /**< \ref TU environment. */
  TU_HASHTABLE* hashtable,  /**< Hash table. */
  const void* keyArray,           /**< First byte of key array. */
  size_t keyLength,         /**< Length of key array. */
  const void* value               /**< Value to be set. */
);




#ifdef __cplusplus
}
#endif

#endif /* TU_HASHTABLE_INTERNAL_H */
