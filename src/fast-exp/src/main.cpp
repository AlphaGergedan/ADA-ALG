#include "fast_exp.hpp"
#include <iostream>

int main() {
  std::cout << "81^131 mod 527 = " << fastExp(81, 131, 527) << std::endl;
  return 0;
}
