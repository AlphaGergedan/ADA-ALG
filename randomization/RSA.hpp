#ifndef _RSA_HPP_
#define _RSA_HPP_

/*  Encryption in a public system
 *  1. A receives B's public key P_b from a public directory or directly from B
 *  2. A computes the ciphertext C = P_b(M) and sends it to B
 *  3. After receiving message C, B decrypts the message using his secret key S_b:
 *     M = S_b(C)
 *
 *  Creating a digital signature
 *  1. A computes the digital signature sigma for M' using her secret key:
 *     sigma = S_a(M')
 *  2. A sends the pair (M', sigma) to B
 *  3. After receiving (M', sigma), B checks the digital signature:
 *     P_a(sigma) = M'
 *  Anybody is able to check sigma using P_a
 *
 *  RSA cryptosystem : Generating the public and secret keys
 *  ----------------
 *  1. Select at random two large primes p and q of l+1 bits (l>2000)
 *  2. Compute n = pq
 *  3. Select a natural number e that is relatively prime to (p-1)(q-1)
 *  4. Compute d = e^-1 (1 mod (p-1)(q-1))
 *  Public key P = (e,n), Secret key S = (d,n)
 *  P(M) = M^e mod n, S(C) = C^d mod n, S(P(M)) = P(S(M)) = M^ed mod n = M,
 *  for any 0 <= M <= 2^2l
 */

/* Extended euclid to compute the inverse d:
 * --> call extendedEuclid( (p-1)(q-1), e) then we have d = y since
 *     algorithm outputs x and y such that (p-1)(q-1)x + ey = 1
 */
#include "../extended-euclid/src/extended_euclid.hpp"
#include <gmp.h>

/**
 *  RSA cryptosystem : Generating the public and secret keys
 *
 *  1. Select at random two large primes p and q of l+1 bits (l>2000)
 *  2. Compute n = pq
 *  3. Select a natural number e that is relatively prime to (p-1)(q-1)
 *  4. Compute d = e^-1 (1 mod (p-1)(q-1))
 *  Public key P = (e,n), Secret key S = (d,n)
 *  P(M) = M^e mod n, S(C) = C^d mod n, S(P(M)) = P(S(M)) = M^ed mod n = M,
 *  for any 0 <= M <= 2^2l
 */
void RSA(mpz_t n, mpz_t e, mpz_t d);

#endif
