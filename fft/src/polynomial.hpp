#ifndef POLY_HPP
#define POLY_HPP

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <math.h>

#define EPSILON 0.001f

class Polynomial {
public:
  virtual void toString(int precision);
  int getDegree() { return this->n; };
  virtual ~Polynomial();

protected:
  /** degree of the polynomial */
  int n = -1;
};

class Coef : virtual Polynomial {
public:
  /**
   * Constructor for coefficient representation
   * p(x) = a_n * x^n + ... + a_1 * x^1 + a_0
   *
   * @param deg: degree of the polynomial              (n)
   * @param a: array of coefficients of the polynomial [n+1]
   */
  Coef(int deg, float a[]) {
    this->n = deg; // FIXME should be >= 0
    this->a = new float[deg+1];
    for (int i = 0; i <= deg; i++) {
      this->a[i] = a[i];
    }
  }
  ~Coef() {
    delete this->a;
  }

  /**
   * Multiplication with a constant is in O(n).
   */
  Coef& operator *(float m) {
    float *a = new float[this->getDegree()+1];
    for (int i = 0; i <= this->getDegree(); i++) {
      a[i] *= this->getCoefs()[i] * m;
    }
    Coef *retVal = new Coef(this->getDegree(), a);
    return *retVal;
  }
  Coef& operator -() {
    return (*this) * (-1);
  }

  /**
   * Addition in coefficient form is in O(maxDegree).
   * Simply add the coefficients.
   */
  Coef& operator +(Coef &p) {
    int maxDegree = std::max(this->getDegree(), p.getDegree());
    int minDegree = std::min(this->getDegree(), p.getDegree());
    float *a = new float[maxDegree+1];

    int i = 0;
    /* add the coefficients */
    while (i <= minDegree) {
      a[i] = this->getCoefs()[i] + p.getCoefs()[i];
      i++;
    }
    /* append the remaining coefs */
    if (i <= this->getDegree()) {
      while (i <= this->getDegree()) {
        a[i] = this->getCoefs()[i];
        i++;
      }
    } else {
      while (i <= p.getDegree()) {
        a[i] = p.getCoefs()[i];
        i++;
      }
    }
    Coef *retVal = new Coef(maxDegree, a);
    return *retVal;
  }

  /**
   * Add the negated polynomial. multiplication with a constant and addition
   * of polynomials can be done in linear time, therefore substraction is also
   * in O(maxDegree);
   */
  Coef& operator -(Coef &p) {
    Coef &pNeg = -p;
    Coef &retVal = (*this) + pNeg;
    delete &pNeg;
    return retVal;
  }

  /**
   * Naive poly multiplication is in O(n^2).
   */
  Coef& operator *(Coef &p) {
    /* compute the degree */
    int multDegree = this->getDegree() + p.getDegree();
    float *a = new float[multDegree + 1];
    for (int i = 0; i <= multDegree; i++) {
      a[i] = 0;
    }
    for (int i = 0; i <= p.getDegree(); i++) {
      for (int j = 0; j <= this->getDegree(); j++) {
        a[i+j] += p.getCoefs()[i] * this->getCoefs()[j];
      }
    }
    Coef *retVal = new Coef(multDegree, a);
    return *retVal;
  }

  /**
   * Compares the coefficients and degrees, returns true if they match
   */
  bool operator ==(Coef &p) {
    if (this->getDegree() != p.getDegree()) {
      return false;
    }
    bool isEqual = true;
    for (int i = 0; i <= p.getDegree(); i++) {
      bool isCloseEnough = std::fabs(p.getCoefs()[i] - this->getCoefs()[i]) < EPSILON;
      isEqual = isEqual && isCloseEnough;
    }
    return isEqual;
  }

  /* evaluates the polynomial at n, using Horner's method in O(n)
   * p(x) = a_0 + x(a_1 + x(a_2 + x(a_3 .. (a_n-1 + x a_n ).. )))*/
  float eval(float x) {
    // constant polynomial
    if (this->getDegree() == 0) {
      return this->getCoefs()[0];
    }

    // a_n
    float retVal = this->getCoefs()[this->getDegree()];
    for (int i = this->getDegree() - 1 ; i >= 0; i--) {
      /* a_n-1 + x*(a_n) */
      retVal = this->getCoefs()[i] + x*retVal;
    }
    return retVal;
  }

  void toString(int precision) {
    int i = 0;
    while (i < this->getDegree()) {
      if (this->getCoefs()[i] > 0) {
        std::cout << std::fixed << std::setprecision(precision)
                  << this->getCoefs()[i] << "*x^" << i << " + ";
      }
      i++;
    }
    if (this->getCoefs()[i] > 0) {
      std::cout << std::fixed << std::setprecision(precision)
                << this->getCoefs()[i] << "*x^" << i;
    }
    std::cout << std::flush;
  }

  float* getCoefs() { return this->a; };

private:
  /* coefficient representation,
   * a[i] corresponds to the coefficient a_i * x^i */
  float *a = nullptr;
};

class LinFac : virtual Polynomial {
public:
  /**
   * Constructor for linear factors representation
   * p(x) = b * (x - x_1) * ... * (x - x_n)
   *
   * @param deg: degree of the polynomial (n)
   * @param b: initial coefficient
   * @param x: array of linear factors    [n]
   */
  LinFac(int deg, float b, float x[]);

  /* TODO: addition ?, multiplication O(1), evaluation O(n), substraction ? */

private:
  /* product of linear factors */
  float b = 0, *x = nullptr;
};

class PV : virtual Polynomial {
public:
  /**
   * Constructor for point-value representation
   * Any polynomial of degree n is uniquely defined by
   * n+1 pairs (x_i,p(x_i)), where i = 0,..,n and x_i != x_j for i != j
   *
   * @param deg: degree of the polynomial (n)
   * @param v: array of x_i values        [n+1]
   * @param y: array of p(x_i) values     [n+1]
   */
  PV(int deg, float v[], float y[]);

  /* TODO: addition O(n), multiplication O(n), evaluation ?, substraction O(n) */
  /* operations on polynomials */
  Polynomial& operator +(const Polynomial &p);
  Polynomial& operator -(const Polynomial &p);
  Polynomial& operator *(const Polynomial &p);
  bool operator ==(const Polynomial &p);

  float eval(float x);

  /* interpolate */
  Coef& toCoef();

private:
  /* point-value representation */
  float *v = nullptr, *y = nullptr;

  Coef& interpolate();
};

/**
 * Polynomial multiplication using FFT
 * TODO
 *
 * computes DFT of p = (p(w^0),p(w^1),..,p(w^n-1))
 * Converts coefficient representation to
 * point-value representation
 */
Coef& mult_fft(Coef& p, Coef& q) {
  // TODO
  return *(new Coef(0,nullptr));
}

#endif
