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
  template <typename T>
  T gcd_impl(T a, T b, T& s, T& t)
  {
    assert(a >= 0 && b >= 0 && a >= b);

    if (b == 0)
    {
      s = 1;
      t = 0;
      return a;
    }

    T q = a / b;
    T r = a % b;
    T result = gcd_impl(b, r, t, s);
    t -= q * s;
    return result;
  }

  template <typename T>
  T gcd(T a, T b, T& s, T& t)
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
      T result;
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

  template <typename T>
  T gcd(T a, T b)
  {
    T s, t;
    return gcd(a, b, s, t);
  }
}

#endif /* GCD_HPP_ */
