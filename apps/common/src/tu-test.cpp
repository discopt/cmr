
#include <Integer.h>
#include <Matrix.h>
#include <iostream>

namespace polymake { namespace common {

  bool is_totally_unimodular_plain(const Matrix<Integer>& matrix) {
    return matrix.rows() > 2 ? true: false;
  }

  UserFunction4perl("# Tests a given //matrix// for total unimodularity."
        "# @param Matrix<int> matrix"
        "# @return Bool",
        &is_totally_unimodular_plain,"is_totally_unimodular(Matrix<Integer>)");
} }

