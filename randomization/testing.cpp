#include <random>
#include <iostream>
#include <iomanip>
#include <stdlib.h>


/**
 * Generates random number.
 * @param a
 * @param b
 * @return random number in [a,b]
 */
uint64_t generateRandom(uint64_t a, uint64_t b) {
  std::random_device dev;
  std::mt19937 rng(dev()); // A Mersenne Twister pseudo-random generator
  std::uniform_int_distribution<std::mt19937::result_type> range(a,b); // distribution in range [a,b]
  return range(rng);
}



int main(int argc, char **argv) {
  /* creating n random numbers between 0 and 9 to see the probabilities */
  uint64_t n = std::stoul(argv[1]);
  uint64_t size = 100;

  /* occurence numbers */
  uint64_t *t = new uint64_t[size];
  for (int i = 0; i < size; i++) {
    t[i] = 0;
  }

  uint64_t c = n;
  while (c > 0) {
    uint64_t rand = generateRandom(0, size-1);
    t[rand] += 1;
    c--;
  }

  std::cout << std::fixed << std::setprecision(5) << "p = [";
  /* probabilities */
  double *p = new double[size];
  for (int i = 0; i < size; i++) {
    p[i] = static_cast<double>(t[i]) / n;
    std::cout << p[i] << ", ";
  }
  std::cout << std::endl;
  return 0;
}
