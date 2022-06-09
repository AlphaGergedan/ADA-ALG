#include "../fast-exp/src/fast_exp.hpp"
#include <cmath>
#include <random>
#include <ctime>

/**
 * only used for Miller Rabin primality test
 */
static bool isProbablyPrime_GLOBAL = true;;

/**
 * Naively check if an odd number in the range [1, floor(sqrt(n))] divides n
 * If yes then n is not a prime, if no n is a prime.
 *
 * Running time: O(sqrt(n))
 *
 * @param n >= 2
 * @param bool isPrime(n)
 */
bool isPrime_det_naive(uint64_t n);

/**
 * For any odd (!= 2) prime number p: 2^(p-1) mod p = 1
 *
 * Running time: O(log^3(n)) = #recursive calls * log^2 n for the bit operations
 * This is polynomial in the input length.
 *
 * Check the above condition and return probably prime (true) if it holds.
 *
 * If the algorithm returns false, then n is definitely not prime.
 *
 * @param n >= 2
 * @param bool isProbablyPrime(n)
 */
bool isPrime_nondet_simple(uint64_t n);

/**
 * Fermat's Little Theorem: if p is prime and 0 < a < p then
 *      a^(p-1) mod p = 1
 *
 * Choose a in the range [2, n-1] uniformly at random,
 * check the above condition and return probably prime (true) if it holds.
 *
 * Running time: O(log^3(n)) = #recursive calls * log^2 n for the bit operations
 * This is polynomial in the input length.
 *
 * If the algorithm returns false, then n is definitely not prime.
 *
 * @param n >= 2
 * @param bool isProbablyPrime(n)
 */
bool isPrime_nondet_randomized(uint64_t n);

/**
 * Fermat's Little Theorem: if p is prime and 0 < a < p then
 *      a^(p-1) mod p = 1
 *
 * Choose a in the range [2, n-1] uniformly at random,
 * check the above condition and return probably prime (true) if it holds.
 * While computing the exponent, check if a non-trivial square root mod p exists.
 *
 * A number a is a non-trivial square root mod n if a^2 mod n = 1 and a != 1 or n-1
 * If p is prime and 0 < a < p, then the equation a^2 mod p = 1 has exactly the two
 * solutions a = 1 and a = p-1, so if non-trivial square root mod n exists, then p
 * cannot be prime.
 *
 * Running time: O(log^3(n)) = #recursive calls * log^2 n for the bit operations
 * This is polynomial in the input length.
 *
 * If the algorithm returns false, then n is definitely not prime.
 *
 * This algorithm has success probability >= 3/4 i.e. error probability < 1/4
 *
 * @param n >= 2
 * @param bool isProbablyPrime(n)
 */
bool isPrime_MillerRabin(uint64_t n);

/**
 * @param long a
 * @param long b
 * @return true iff a divides b
 */
bool divides(uint64_t a, uint64_t b) {
  return ( (b % a) == 0);
}

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

bool isPrime_det_naive(uint64_t n) {
  bool isPrime;
  if (n < 2) {
    isPrime = false;
  } else if (n == 2) {
    isPrime = true;
  } else {
    isPrime = true;
    /* n even ? */
    if (n % 2 == 0) {
      isPrime = false;
    } else {
      /* check if odd numbers between 1, .. , sqrt(n) divides n */
      uint64_t root_n_div2 = std::floor(std::sqrt(n) / 2);
      for(uint64_t i = 1; i <= root_n_div2; i++) {
        if (divides(2*i + 1, n)) {
          isPrime = false;
        }
      }
    }
  }
  return isPrime;
}

bool isPrime_nondet_simple(uint64_t n) {
  bool isProbablyPrime;
  if (n < 2) {
    isProbablyPrime = false;
  } else if (n == 2) {
    isProbablyPrime = true;
  } else {
    isProbablyPrime = true;
    /* n even ? */
    if (n % 2 == 0) {
      isProbablyPrime = false;
    } else {
        /* up to this point all the checks can be done in constant time */
        /* 2^n-1 mod n */
        uint64_t z = fastExp(2, n-1, n);
        if (z == 1) {
          isProbablyPrime = true;
        } else {
          isProbablyPrime = false;
        }
    }
  }
  return isProbablyPrime;
}

bool isPrime_nondet_randomized(uint64_t n) {
  bool isProbablyPrime;
  if (n < 2) {
    isProbablyPrime = false;
  } else if (n == 2) {
    isProbablyPrime = true;
  } else {
    isProbablyPrime = true;
    /* n even ? */
    if (n % 2 == 0) {
      isProbablyPrime = false;
    } else {
        /* up to this point all the checks can be done in constant time */
        /* a in the range [2, n-1] uniformly at random */
        uint64_t a = generateRandom(2, n-1);
        /* a^n-1 mod n */
        uint64_t z = fastExp(a, n-1, n);
        if (z == 1) {
          isProbablyPrime = true;
        } else {
          isProbablyPrime = false;
        }
    }
  }
  return isProbablyPrime;
}

/**
 * Same as fast exponent but with checking if
 * non-trivial square root mod n exists.
 * Computes a^p mod n, sets isProbablyPrime to false if
 * non-trivial square root mod n exists.
 *
 * @param a
 * @param p
 * @param n
 * @return a^p mod n, sets isProbablyPrime accordingly
 */
uint64_t power(uint64_t a, uint64_t p, uint64_t n) {
  if (p == 0) {
    return 1;
  }
  /* If we assume p is even, we get:
   * a^p mod n = a^(p/2) * a^(p/2) mod n
   *           = ((a^(p/2)) mod n) * ((a^(p/2)) mod n) mod n */
  uint64_t x = fastExp(a, p/2, n);
  uint64_t res = (x*x) % n;
  /* check if x^2 mod n = 1 and x != 1 and x != n-1 */
  if (res == 1 && x != 1 && x != n-1) {
    isProbablyPrime_GLOBAL = false;
  }
  /* Add the remaining power if p is odd */
  if (!(p % 2 == 0)) {
    res = (a * res) % n;
  }
  return res;
}

bool isPrime_MillerRabin(uint64_t n) {
  bool isProbablyPrime;
  if (n < 2) {
    isProbablyPrime = false;
  } else if (n == 2) {
    isProbablyPrime = true;
  } else {
    isProbablyPrime = true;
    /* n even ? */
    if (n % 2 == 0) {
      isProbablyPrime = false;
    } else {
        /* up to this point all the checks can be done in constant time */
        /* a in the range [2, n-1] uniformly at random */
        uint64_t a = generateRandom(2, n-1);
        /* a^n-1 mod n */
        uint64_t z = power(a, n-1, n);
        if (z != 1 || !isProbablyPrime_GLOBAL) {
          isProbablyPrime = false;
        }
        else {
          isProbablyPrime = true;
        }
    }
  }
  return isProbablyPrime;
}
