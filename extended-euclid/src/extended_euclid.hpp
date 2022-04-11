#ifndef EXT_EUC_HPP
#define EXT_EUC_HPP

#include <tuple>
#include <gmp.h>

/**
 * Extended Euclid Algorithm : O(log n). Naive GCD is in O(n).
 * We modify the algorithm a little to also return two additional
 * integers with the following property.
 *
 * @param a positive integer
 * @param b positive integer
 * @return (GCD(a,b), x, y), where GCD(a,b) is the greatest common divisor of
 *         a and b, x and y are two integers with x*a + y*b = GCD(a,b)
 */
std::tuple<int,int,int> extendedEuclid(int a, int b);
void extendedEuclid(mpz_t a, mpz_t b, mpz_t gcd, mpz_t x, mpz_t y);

#endif
