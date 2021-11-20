#include "extended_euclid.hpp"
#include <iostream>
#include<assert.h>

int main() {
  int p = 31;
  int q = 17;
  int e = 131;

  /* e needs to be relatively prime to (p-1)*(q-1) */

  std::tuple<int,int,int> res = extendedEuclid((p-1)*(q-1), e);

  assert(1 == std::get<0>(res));
  assert(std::get<1>(res)*((p-1)*(q-1)) + std::get<2>(res)*(e) == std::get<0>(res));

  std::cout << "x = " << std::get<1>(res) << std::endl
            << "y = " << std::get<2>(res) << std::endl;
  return 0;
}
