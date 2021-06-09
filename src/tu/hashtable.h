#ifndef TU_HASHTABLE_INTERNAL_H
#define TU_HASHTABLE_INTERNAL_H

#include "env_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _TU_LINEARHASHTABLE_ARRAY TU_LINEARHASHTABLE_ARRAY;
typedef size_t TU_LINEARHASHTABLE_BUCKET;
typedef size_t TU_LINEARHASHTABLE_HASH;

TU_ERROR TUlinearhashtableArrayCreate(
  TU* tu,                     /**< \ref TU environment. */
  TU_LINEARHASHTABLE_ARRAY** phashtable,  /**< Pointer for storing the hash table. */
  size_t initialSize,         /**< Initial size of hash table. */
  size_t initialKeyMemory     /**< Initial memory for keys. */  
);

TU_ERROR TUlinearhashtableArrayFree(
  TU* tu,                   /**< \ref TU environment. */
  TU_LINEARHASHTABLE_ARRAY** phashtable /**< Pointer to the hash table. */
);

bool TUlinearhashtableArrayFind(
  TU_LINEARHASHTABLE_ARRAY* hashtable,  /**< Hash table. */
  const void* keyArray,               /**< First byte of key array. */
  size_t keyLength,                   /**< Length of key array in bytes. */
  TU_LINEARHASHTABLE_BUCKET* pbucket,   /**< Pointer for storing the bucket. */
  TU_LINEARHASHTABLE_HASH* phash        /**< Pointer for storing the hash. */
);

const void* TUlinearhashtableArrayKey(
  TU_LINEARHASHTABLE_ARRAY* hashtable,  /**< Hash table. */
  TU_LINEARHASHTABLE_BUCKET bucket,     /**< Bucket. */
  size_t* pKeyLength                  /**< Length of key array. */
);

const void* TUlinearhashtableArrayValue(
  TU_LINEARHASHTABLE_ARRAY* hashtable,  /**< Hash table. */
  TU_LINEARHASHTABLE_BUCKET bucket      /**< Bucket. */
);

TU_ERROR TUlinearhashtableArrayInsertBucketHash(
  TU* tu,                             /**< \ref TU environment. */
  TU_LINEARHASHTABLE_ARRAY* hashtable,  /**< Hash table. */
  const void* keyArray,               /**< First byte of key array. */
  size_t keyLength,                   /**< Length of key array. */
  TU_LINEARHASHTABLE_BUCKET bucket,     /**< Known bucket of key, determined by \ref TUlinearhashtableArrayFind. */
  TU_LINEARHASHTABLE_HASH hash,         /**< Known hash of key, determined by \ref TUlinearhashtableArrayFind. */
  const void* value                   /**< Value to be set. */
);

TU_ERROR TUlinearhashtableArrayInsert(
  TU* tu,                   /**< \ref TU environment. */
  TU_LINEARHASHTABLE_ARRAY* hashtable,  /**< Hash table. */
  const void* keyArray,     /**< First byte of key array. */
  size_t keyLength,         /**< Length of key array. */
  const void* value         /**< Value to be set. */
);



typedef struct _TU_LISTHASHTABLE TU_LISTHASHTABLE;
typedef size_t TU_LISTHASHTABLE_VALUE;
typedef size_t TU_LISTHASHTABLE_HASH;
typedef size_t TU_LISTHASHTABLE_BUCKET;
typedef size_t TU_LISTHASHTABLE_ENTRY;

TU_ERROR TUlisthashtableCreate(
  TU* tu,                         /**< \ref TU environment. */
  TU_LISTHASHTABLE** phashtable,  /**< Pointer for storing the hash table. */
  size_t initialNumBuckets,       /**< Initial size of hash table. */
  size_t initialMemNodes          /**< Initial memory for actual entries.. */
);

TU_ERROR TUlisthashtableFree(
  TU* tu,                       /**< \ref TU environment. */
  TU_LISTHASHTABLE** phashtable /**< Pointer to the hash table. */
);

TU_LISTHASHTABLE_ENTRY TUlisthashtableFindFirst(
  TU_LISTHASHTABLE* hashtable,  /**< Hash table. */
  TU_LISTHASHTABLE_HASH hash    /**< Hash value. */
);

TU_LISTHASHTABLE_ENTRY TUlisthashtableFindNext(
  TU_LISTHASHTABLE* hashtable,  /**< Hash table. */
  TU_LISTHASHTABLE_HASH hash,   /**< Hash value. */
  TU_LISTHASHTABLE_ENTRY entry  /**< Current entry. */
);

TU_LISTHASHTABLE_VALUE TUlisthashtableValue(
  TU_LISTHASHTABLE* hashtable,  /**< Hash table. */
  TU_LISTHASHTABLE_ENTRY entry  /**< Entry. */
);

TU_LISTHASHTABLE_HASH TUlisthashtableHash(
  TU_LISTHASHTABLE* hashtable,  /**< Hash table. */
  TU_LISTHASHTABLE_ENTRY entry  /**< Entry. */
);

/**
 * \brief Inserts \p value with \p hash into the hash table, even if an entry with the same hash exists.
 */

TU_ERROR TUlisthashtableInsert(
  TU* tu,                         /**< \ref TU environment. */
  TU_LISTHASHTABLE* hashtable,    /**< Hash table. */
  TU_LISTHASHTABLE_HASH hash,     /**< Bucket. */
  TU_LISTHASHTABLE_VALUE value,   /**< Value to be set. */
  TU_LISTHASHTABLE_ENTRY* pentry  /**< Pointer for storing the new entry. */
);

TU_ERROR TUlisthashtableRemove(
  TU* tu,                       /**< \ref TU environment. */
  TU_LISTHASHTABLE* hashtable,  /**< Hash table. */
  TU_LISTHASHTABLE_ENTRY entry  /**< Entry to be removed. */
);


#ifdef __cplusplus
}
#endif

#endif /* TU_HASHTABLE_INTERNAL_H */
