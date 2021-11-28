#include "fast_exp.hpp"

int fastExp(int a, int p, int n) {
  if (p == 0) {
    return 1;
  }
  /* If we assume p is even, we get:
   * a^p mod n = a^(p/2) * a^(p/2) mod n
   *           = ((a^(p/2)) mod n) * ((a^(p/2)) mod n) mod n */
  int x = fastExp(a, p/2, n);
  int res = (x*x) % n;
  /* Add the remaining power if p is odd */
  if (!(p % 2 == 0)) {
    res = (a * res) % n;
  }
  return res;
}
