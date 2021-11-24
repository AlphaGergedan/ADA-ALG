#include "extended_euclid.hpp"

/**
 * Extended Euclid Algorithm : O(log n). Naive GCD is in O(n).
 * We modify the algorithm a little to also return two additional
 * integers with the following property.
 *
 * Note (*) that: GCD(a,b) = GDC(b, a mod b)
 *
 * @param a positive integer
 * @param b positive integer
 * @return (GCD(a,b), x, y), where GDC(a,b) is the greatest common divisor of
 *         a and b, x and y are two integers with x*a + y*b = GDC(a,b)
 */
std::tuple<int,int,int> extendedEuclid(int a, int b) {
  /* Base case */
  if (!b) {
    return std::make_tuple(a, 1, 0);
  }
  /* Recursive call with note (*) */
  std::tuple<int,int,int> rec = extendedEuclid(b, a % b);
  int x = std::get<2>(rec);
  int y = std::get<1>(rec) - (a/b)*x;
  return std::make_tuple(std::get<0>(rec), x, y);
}
