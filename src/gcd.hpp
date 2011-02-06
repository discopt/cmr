/*
 * gcd.hpp
 *
 *  Created on: Feb 5, 2011
 *      Author: xammy
 */

#ifndef GCD_HPP_
#define GCD_HPP_

namespace unimod
{
  int gcd_impl(int a, int b, int& s, int& t)
  {
    assert(a >= 0 && b >= 0 && a >= b);

    if (b == 0)
    {
      s = 1;
      t = 0;
      return a;
    }

    int q = a / b;
    int r = a % b;
    int result = gcd_impl(b, r, t, s);
    t -= q * s;
    return result;
  }

  int gcd(int a, int b, int& s, int& t)
  {
    if (a >= 0 && b >= 0)
    {
      if (a >= b)
        return gcd_impl(a, b, s, t);
      else
        return gcd_impl(b, a, t, s);
    }
    else
    {
      bool a_neg = a < 0;
      bool b_neg = b < 0;
      a = a >= 0 ? a : -a;
      b = b >= 0 ? b : -b;
      int result;
      if (a >= b)
        result = gcd_impl(a, b, s, t);
      else
        result = gcd_impl(b, a, t, s);
      if (a_neg)
        s = -s;
      if (b_neg)
        t = -t;
      return result;
    }
  }

  int gcd(int a, int b)
  {
    int s, t;
    return gcd(a, b, s, t);
  }
}

#endif /* GCD_HPP_ */
