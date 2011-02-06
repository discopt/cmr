/*
 * smith_normal_form.hpp
 *
 *  Created on: Feb 3, 2011
 *      Author: xammy
 */

#ifndef SMITH_NORMAL_FORM_HPP_
#define SMITH_NORMAL_FORM_HPP_

#include "common.hpp"

namespace unimod
{
  int gcd(int a, int b, int& s, int& t);

  void smith_normal_form(const integer_matrix& matrix, std::vector<int>& diagonal);
}

#endif /* SMITH_NORMAL_FORM_HPP_ */
