#ifndef FAST_EXP_HPP
#define FAST_EXP_HPP

#include <cstdint>

/**
 * Computes a^p mod n
 *
 * @param a
 * @param p
 * @param n
 * @return a^p mod n
 */
uint64_t fastExp(uint64_t a, uint64_t p, uint64_t n);

#endif
