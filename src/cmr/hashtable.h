#ifndef CMR_HASHTABLE_INTERNAL_H
#define CMR_HASHTABLE_INTERNAL_H

#include "env_internal.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _CMR_LINEARHASHTABLE_ARRAY CMR_LINEARHASHTABLE_ARRAY;
typedef size_t CMR_LINEARHASHTABLE_BUCKET;
typedef size_t CMR_LINEARHASHTABLE_HASH;

CMR_ERROR CMRlinearhashtableArrayCreate(
  CMR* cmr,                     /**< \ref CMR environment. */
  CMR_LINEARHASHTABLE_ARRAY** phashtable,  /**< Pointer for storing the hash table. */
  size_t initialSize,         /**< Initial size of hash table. */
  size_t initialKeyMemory     /**< Initial memory for keys. */  
);

CMR_ERROR CMRlinearhashtableArrayFree(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_LINEARHASHTABLE_ARRAY** phashtable /**< Pointer to the hash table. */
);

bool CMRlinearhashtableArrayFind(
  CMR_LINEARHASHTABLE_ARRAY* hashtable,  /**< Hash table. */
  const void* keyArray,               /**< First byte of key array. */
  size_t keyLength,                   /**< Length of key array in bytes. */
  CMR_LINEARHASHTABLE_BUCKET* pbucket,   /**< Pointer for storing the bucket. */
  CMR_LINEARHASHTABLE_HASH* phash        /**< Pointer for storing the hash. */
);

const void* CMRlinearhashtableArrayKey(
  CMR_LINEARHASHTABLE_ARRAY* hashtable,  /**< Hash table. */
  CMR_LINEARHASHTABLE_BUCKET bucket,     /**< Bucket. */
  size_t* pKeyLength                  /**< Length of key array. */
);

const void* CMRlinearhashtableArrayValue(
  CMR_LINEARHASHTABLE_ARRAY* hashtable,  /**< Hash table. */
  CMR_LINEARHASHTABLE_BUCKET bucket      /**< Bucket. */
);

CMR_ERROR CMRlinearhashtableArrayInsertBucketHash(
  CMR* cmr,                             /**< \ref CMR environment. */
  CMR_LINEARHASHTABLE_ARRAY* hashtable,  /**< Hash table. */
  const void* keyArray,               /**< First byte of key array. */
  size_t keyLength,                   /**< Length of key array. */
  CMR_LINEARHASHTABLE_BUCKET bucket,     /**< Known bucket of key, determined by \ref CMRlinearhashtableArrayFind. */
  CMR_LINEARHASHTABLE_HASH hash,         /**< Known hash of key, determined by \ref CMRlinearhashtableArrayFind. */
  const void* value                   /**< Value to be set. */
);

CMR_ERROR CMRlinearhashtableArrayInsert(
  CMR* cmr,                   /**< \ref CMR environment. */
  CMR_LINEARHASHTABLE_ARRAY* hashtable,  /**< Hash table. */
  const void* keyArray,     /**< First byte of key array. */
  size_t keyLength,         /**< Length of key array. */
  const void* value         /**< Value to be set. */
);



typedef struct _CMR_LISTHASHTABLE CMR_LISTHASHTABLE;
typedef size_t CMR_LISTHASHTABLE_VALUE;
typedef size_t CMR_LISTHASHTABLE_HASH;
typedef size_t CMR_LISTHASHTABLE_BUCKET;
typedef size_t CMR_LISTHASHTABLE_ENTRY;

CMR_ERROR CMRlisthashtableCreate(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_LISTHASHTABLE** phashtable,  /**< Pointer for storing the hash table. */
  size_t initialNumBuckets,       /**< Initial size of hash table. */
  size_t initialMemNodes          /**< Initial memory for actual entries.. */
);

CMR_ERROR CMRlisthashtableFree(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_LISTHASHTABLE** phashtable /**< Pointer to the hash table. */
);

CMR_LISTHASHTABLE_ENTRY CMRlisthashtableFindFirst(
  CMR_LISTHASHTABLE* hashtable,  /**< Hash table. */
  CMR_LISTHASHTABLE_HASH hash    /**< Hash value. */
);

CMR_LISTHASHTABLE_ENTRY CMRlisthashtableFindNext(
  CMR_LISTHASHTABLE* hashtable,  /**< Hash table. */
  CMR_LISTHASHTABLE_HASH hash,   /**< Hash value. */
  CMR_LISTHASHTABLE_ENTRY entry  /**< Current entry. */
);

CMR_LISTHASHTABLE_VALUE CMRlisthashtableValue(
  CMR_LISTHASHTABLE* hashtable,  /**< Hash table. */
  CMR_LISTHASHTABLE_ENTRY entry  /**< Entry. */
);

CMR_LISTHASHTABLE_HASH CMRlisthashtableHash(
  CMR_LISTHASHTABLE* hashtable,  /**< Hash table. */
  CMR_LISTHASHTABLE_ENTRY entry  /**< Entry. */
);

size_t CMRlisthashtableNumBuckets(
  CMR_LISTHASHTABLE* hashtable /**< Hash table. */
);

/**
 * \brief Inserts \p value with \p hash into the hash table, even if an entry with the same hash exists.
 */

CMR_ERROR CMRlisthashtableInsert(
  CMR* cmr,                         /**< \ref CMR environment. */
  CMR_LISTHASHTABLE* hashtable,    /**< Hash table. */
  CMR_LISTHASHTABLE_HASH hash,     /**< Bucket. */
  CMR_LISTHASHTABLE_VALUE value,   /**< Value to be set. */
  CMR_LISTHASHTABLE_ENTRY* pentry  /**< Pointer for storing the new entry. */
);

CMR_ERROR CMRlisthashtableRemove(
  CMR* cmr,                       /**< \ref CMR environment. */
  CMR_LISTHASHTABLE* hashtable,  /**< Hash table. */
  CMR_LISTHASHTABLE_ENTRY entry  /**< Entry to be removed. */
);


#ifdef __cplusplus
}
#endif

#endif /* CMR_HASHTABLE_INTERNAL_H */
