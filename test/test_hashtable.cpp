#include <gtest/gtest.h>

#include "common.h"
#include "../src/cmr/hashtable.h"

TEST(Hashtable, Functionality)
{
  TU* tu = NULL;
  ASSERT_TU_CALL( TUcreateEnvironment(&tu) );

  TU_HASHTABLE* hashtable = NULL;
  ASSERT_TU_CALL( TUhashtableCreate(tu, &hashtable, 16, 16) );

  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "q", strlen("q"), (void*) 17) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "a", strlen("a"), (void*) 1) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "b", strlen("b"), (void*) 2) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "r", strlen("r"), (void*) 18) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "c", strlen("c"), (void*) 3) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "d", strlen("d"), (void*) 4) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "e", strlen("e"), (void*) 5) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "f", strlen("f"), (void*) 6) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "g", strlen("g"), (void*) 7) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "h", strlen("h"), (void*) 8) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "i", strlen("i"), (void*) 9) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "j", strlen("j"), (void*) 10) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "k", strlen("k"), (void*) 11) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "l", strlen("l"), (void*) 12) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "m", strlen("m"), (void*) 13) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "n", strlen("n"), (void*) 14) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "o", strlen("o"), (void*) 15) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "p", strlen("p"), (void*) 16) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "s", strlen("s"), (void*) 19) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "t", strlen("t"), (void*) 20) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "u", strlen("u"), (void*) 21) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "v", strlen("v"), (void*) 22) );
  ASSERT_TU_CALL( TUhashtableInsert(tu, hashtable, "w", strlen("w"), (void*) 23) );

  ASSERT_TU_CALL( TUhashtableFree(tu, &hashtable) );

  TUfreeEnvironment(&tu);
}
