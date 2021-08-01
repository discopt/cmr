#include <cmr/matrix.hpp>
#include "matrix_internal.h"

namespace tu
{

  SubmatrixIndices::SubmatrixIndices(const SubmatrixIndices& other)
    : _rows(other._rows), _columns(other._columns)
  {

  }
  
  SubmatrixIndices::SubmatrixIndices(SubmatrixIndices&& other)
    : _rows(std::move(other._rows)), _columns(std::move(other._columns))
  {

  }

  SubmatrixIndices::SubmatrixIndices(const std::vector<std::size_t>& rows,
    const std::vector<std::size_t>& columns)
    : _rows(rows), _columns(columns)
  {

  }

} /* namespace tu */
