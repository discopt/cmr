#include <gtest/gtest.h>

#include "common.h"
#include "../src/cmr/hashtable.h"

TEST(Hashtable, Functionality)
{
  CMR* cmr = NULL;
  ASSERT_CMR_CALL( CMRcreateEnvironment(&cmr) );

  CMR_LINEARHASHTABLE_ARRAY* hashtable = NULL;
  ASSERT_CMR_CALL( CMRlinearhashtableArrayCreate(cmr, &hashtable, 16, 16) );

  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "q", strlen("q"), (void*) 17) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "a", strlen("a"), (void*) 1) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "b", strlen("b"), (void*) 2) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "r", strlen("r"), (void*) 18) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "c", strlen("c"), (void*) 3) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "d", strlen("d"), (void*) 4) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "e", strlen("e"), (void*) 5) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "f", strlen("f"), (void*) 6) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "g", strlen("g"), (void*) 7) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "h", strlen("h"), (void*) 8) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "i", strlen("i"), (void*) 9) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "j", strlen("j"), (void*) 10) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "k", strlen("k"), (void*) 11) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "l", strlen("l"), (void*) 12) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "m", strlen("m"), (void*) 13) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "n", strlen("n"), (void*) 14) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "o", strlen("o"), (void*) 15) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "p", strlen("p"), (void*) 16) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "s", strlen("s"), (void*) 19) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "t", strlen("t"), (void*) 20) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "u", strlen("u"), (void*) 21) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "v", strlen("v"), (void*) 22) );
  ASSERT_CMR_CALL( CMRlinearhashtableArrayInsert(cmr, hashtable, "w", strlen("w"), (void*) 23) );

  ASSERT_CMR_CALL( CMRlinearhashtableArrayFree(cmr, &hashtable) );

  CMRfreeEnvironment(&cmr);
}
