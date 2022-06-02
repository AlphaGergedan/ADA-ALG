#include "RSA.hpp"

/* Select at random two large primes p and q of l+1 bits (l>2000) */
void initRand_pq(mpz_t p, mpz_t q) { /* TODO */ };

/* TODO: add bit input to specifiy length of the primes */
void RSA(mpz_t n, mpz_t e, mpz_t d) {
  mpz_t p, q;
  mpz_inits(p, q, n, NULL); // initialize p and q and set their value to 0
  /* Select two large random primes */
  initRand_pq(p, q);

  /* Compute n = pq */


  /* Select a natural number e that is relatively prime to (p-1)(q-1) */

  /* Compute d = e^-1 (1 mod (p-1)(q-1)) */

  /* Important to delete the primes */
  mpz_clears(p, q, NULL);
}
