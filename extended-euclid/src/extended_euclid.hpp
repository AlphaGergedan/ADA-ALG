#ifndef EXT_EUC_HPP
#define EXT_EUC_HPP

#include <tuple>

/**
 * Extended Euclid Algorithm
 *
 * @param a positive integer
 * @param b positive integer
 * @return (GCD(a,b), x, y), where GDC(a,b) is the greatest common divisor of
 *         a and b, x and y are two integers with x*a + y*b = GDC(a,b)
 */
std::tuple<int,int,int> extendedEuclid(int a, int b);

#endif
