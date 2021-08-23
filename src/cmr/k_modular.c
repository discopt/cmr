#include <cmr/k_modular.h>

#include "env_internal.h"

#include "interface.h"

CMR_EXPORT
CMR_ERROR CMRtestUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisUnimodular)
{
  size_t k;
  /* TODO: Add parameter that only k=1 can succeed in order to skip computations. */
  CMR_CALL( CMRinterfaceKModular(cmr, matrix, &k) );
  *pisUnimodular = k == 1;

  return CMR_OKAY;
}

CMR_EXPORT
CMR_ERROR CMRtestStrongUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisStronglyUnimodular)
{
  CMR_CALL( CMRtestUnimodularity(cmr, matrix, pisStronglyUnimodular) );
  if (!*pisStronglyUnimodular)
    return CMR_OKAY;

  CMR_CHRMAT* transpose = NULL;
  CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );
  CMR_CALL( CMRtestUnimodularity(cmr, transpose, pisStronglyUnimodular) );
  CMR_CALL( CMRchrmatFree(cmr, &transpose) );

  return CMR_OKAY;
}

CMR_ERROR CMRtestKmodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisKmodular, size_t* pk)
{
  assert(cmr);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );
  assert(pisKmodular);

  size_t k;
  CMR_CALL( CMRinterfaceKModular(cmr, matrix, &k) );
  *pisKmodular = k > 0;
  *pk = k;

  return CMR_OKAY;
  
}

CMR_ERROR CMRtestStrongKmodularity(CMR* cmr, CMR_CHRMAT* matrix, bool* pisStronglyKmodular, size_t* pk)
{
  assert(cmr);
  CMRconsistencyAssert( CMRchrmatConsistency(matrix) );
  assert(pisStronglyKmodular);

  size_t k1, k2;
  CMR_CALL( CMRtestKmodularity(cmr, matrix, pisStronglyKmodular, &k1) );
  if (!*pisStronglyKmodular)
    return CMR_OKAY;

  CMR_CHRMAT* transpose = NULL;
  CMR_CALL( CMRchrmatTranspose(cmr, matrix, &transpose) );
  CMR_CALL( CMRtestKmodularity(cmr, transpose, pisStronglyKmodular, &k2) );
  CMR_CALL( CMRchrmatFree(cmr, &transpose) );

  assert(k1 == k2);
  *pk = k1;

  return CMR_OKAY;
}
