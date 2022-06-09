#include "fast_exp.hpp"

uint64_t fastExp(uint64_t a, uint64_t p, uint64_t n) {
  if (p == 0) {
    return 1;
  }
  /* If we assume p is even, we get:
   * a^p mod n = a^(p/2) * a^(p/2) mod n
   *           = ((a^(p/2)) mod n) * ((a^(p/2)) mod n) mod n */
  uint64_t x = fastExp(a, p/2, n);
  uint64_t res = (x*x) % n;
  /* Add the remaining power if p is odd */
  if (!(p % 2 == 0)) {
    res = (a * res) % n;
  }
  return res;
}
