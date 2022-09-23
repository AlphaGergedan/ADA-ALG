#ifndef POLY_H
#define POLY_H

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <complex>

#define EPSILON 0.001f

const float euler = std::exp(1.0);
typedef std::complex<float> comp;

/**
 * Includes polynomial representations and useful functions such as fft
 */
namespace Polynomial {
  /**
   * Returns n-th root of unity, e^(2*pi*i / n)
   */
  comp getRootUnity(int n) {
    comp w_n;
    if (n > 0) {
      w_n = std::exp<float>((2i*M_PI) / (n));
    } else {
      w_n = 1;
    }
    return w_n;
  }

  /**
   * Base class for polynomials
   */
  class Polynomial {
  public:
    virtual void toString(int precision);
    virtual ~Polynomial();
    int getDegree() { return this->n; };

  protected:
    /** degree of the polynomial */
    int n = -1;
  };

  /**
   * Polynomial in coefficient representation.
   * p(x) = a_n * x^n + ... + a_1 * x^1 + a_0
   */
  template<typename T>
  class Coef : public virtual Polynomial {
  public:
    /**
     * Constructor for coefficient representation
     * p(x) = a_n * x^n + ... + a_1 * x^1 + a_0
     *
     * @param deg: degree of the polynomial              (n)
     * @param a: array of coefficients of the polynomial [n+1]
     */
    Coef(int deg, T a[]) {
      this->n = deg; // FIXME should be >= 0
      this->a = new T[deg+1];
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
    Coef<T>& operator *(T m) {
      T *a = new T[this->getDegree()+1];
      for (int i = 0; i <= this->getDegree(); i++) {
        a[i] *= this->getCoefs()[i] * m;
      }
      Coef<T> *retVal = new Coef<T>(this->getDegree(), a);
      return *retVal;
    }
    Coef<T>& operator -() {
      return (*this) * (-1);
    }

    /**
     * Addition in coefficient form is in O(maxDegree).
     * Simply add the coefficients.
     */
    Coef<T>& operator +(Coef<T> &p) {
      int maxDegree = std::max(this->getDegree(), p.getDegree());
      int minDegree = std::min(this->getDegree(), p.getDegree());
      float *a = new T[maxDegree+1];

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
      Coef<T> *retVal = new Coef<T>(maxDegree, a);
      return *retVal;
    }

    /**
     * Add the negated polynomial. multiplication with a constant and addition
     * of polynomials can be done in linear time, therefore substraction is also
     * in O(maxDegree);
     */
    Coef<T>& operator -(Coef<T> &p) {
      Coef<T> &pNeg = -p;
      Coef<T> &retVal = (*this) + pNeg;
      delete &pNeg;
      return retVal;
    }

    /**
     * Naive poly multiplication is in O(n^2).
     */
    Coef<T>& operator *(Coef<T> &p) {
      /* compute the degree */
      int multDegree = this->getDegree() + p.getDegree();
      T *a = new T[multDegree + 1];
      for (int i = 0; i <= multDegree; i++) {
        a[i] = 0;
      }
      for (int i = 0; i <= p.getDegree(); i++) {
        for (int j = 0; j <= this->getDegree(); j++) {
          a[i+j] += p.getCoefs()[i] * this->getCoefs()[j];
        }
      }
      Coef<T> *retVal = new Coef<T>(multDegree, a);
      return *retVal;
    }

    /**
     * Compares the coefficients and degrees, returns true if they match
     */
    bool operator ==(Coef<T> &p) {
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
    T eval(T x) {
      // constant polynomial
      if (this->getDegree() == 0) {
        return this->getCoefs()[0];
      }

      // a_n
      T retVal = this->getCoefs()[this->getDegree()];
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

    T* getCoefs() { return this->a; };

    /**
     * Fast Fourier Transform (FFT)
     *
     * Computes DFT_n(p) = (p(w_n^0),p(w_n^1),..,p(w_n^n-1))
     * where w = exp(2pi*i/n) where i is imaginary number.
     * Running time is in O(nlogn)
     *
     * @param n : power of 2, number of point value pairs
     * @return DFT_n(p)
     */
    comp* fft(int n) {
      /* this is not allowed, TODO also check if power of 2 */
      if (n < this->getDegree()) {
        return nullptr;
      }
      T *a = new T[n];
      /* initialize coefficient array in O(n) */
      for (int i = 0; i < n; i++) {
        if (i < this->getDegree()) {
          a[i] = this->getCoefs()[i];
        } else {
          a[i] = 0;
        }
      }
      return fft_aux(a, n);
    }


  private:
    /* coefficient representation,
     * a[i] corresponds to the coefficient a_i * x^i */
    T*a = nullptr;

    /**
     *  Running time:
     *      T(1) = O(1)
     *      T(n) = 2T(n/2) + O(n)
     *           = O(nlogn) using e.g. Guessing and Induction
     */
    comp* fft_aux(T*a, int n) {
      /* base case */
      if (n == 1) {
        comp *ret_a = new comp[1];
        ret_a[0] = a[0];
        delete[] a;
        return ret_a;
      }
      /* split even and odd coefficients in O(n) */
      T *a_even = new T[n/2];
      T *a_odd = new T[n/2];
      int cnt = 0;
      for (int i = 0; i < n; i++) {
        // odd
        if (i % 2) {
          a_odd[cnt] = this->getCoefs()[i];
          cnt++;
        } else {
          a_even[cnt] = this->getCoefs()[i];
        }
      }
      /* 2 recursive calls */
      comp *d_0 = fft_aux(a_even, n/2);
      comp *d_1 = fft_aux(a_odd, n/2);
      /* merge step */
      comp *d = new comp[n];
      comp w_n = getRootUnity(n);
      comp w = 1;
      /* in O(n) */
      for (int k = 0; k <= (n/2)-1; k++) {
        d[k] = d_0[k] + w*d_1[k];
        d[k + n/2] = d_0[k] - w*d_1[k];
        w *= w_n;
      }
      /* cleanup */
      delete[] a;
      delete[] a_even;
      delete[] a_odd;
      delete[] d_0;
      delete[] d_1;
      return d;
    }
  };

  /**
   * Linear factors representation (fundamental theorem of algebra)
   * p(x) = b * (x - x_1) * ... * (x - x_n)
   *
   * TODO use union find in == overload
   */
  template <typename T>
  class LinFac : public virtual Polynomial {
  public:
    /**
     * Constructor for linear factors representation
     * p(x) = b * (x - x_1) * ... * (x - x_n)
     *
     * @param deg: degree of the polynomial (n)
     * @param b: initial coefficient
     * @param x: array of linear factors    [n]
     */
    LinFac(int deg, T b, T x[]) {
      this->n = deg;
      this->b = b;
      this->x = new T[deg];
      for (int i = 0; i < deg; i++) {
        this->x[i] = x[i];
      }
    }
    ~LinFac() {
      delete[] this->x;
    }

    /**
     * Multiplication with a constant is in O(1).
     */
    LinFac<T>& operator *(T m) {
      T b = m * this->b;
      T *x = new T[this->getDegree()];
      for (int i = 0; i < this->getDegree(); i++) {
        x[i] = this->getFactors()[i];
      }
      LinFac<T> *retVal = new LinFac<T>(this->getDegree(), b, x);
      return *retVal;
    }
    LinFac<T>& operator -() {
      return (*this) * (-1);
    }

    /**
     * Polynomial multiplication is in O(1).
     */
    LinFac<T>& operator *(LinFac<T> &p) {
      int multDeg = p.getDegree() + this->getDegree();
      T *x = new T[multDeg];
      T b = p.getConstantFactor() * this->getConstantFactor();
      LinFac<T> *retVal = new LinFac<T>(multDeg, b, x);
      return *retVal;
    }

    /**
     * Compares the coefficients and degrees, returns true if they match
     */
    bool operator ==(LinFac<T> &p) {
      LinFac<T> &q = *this;
      if (q.getDegree() != p.getDegree() ||
          q.getConstantFactor() != p.getConstantFactor()) {
        return false;
      }
      // TODO union find datastructure implementation for O(n)
      for (int i = 0; i < p.getDegree(); i++) {
        // add to set A (almost constant amortized time)
      }
      for (int i = 0; i < q.getDegree(); i++) {
        // check if in the set A (almost constant amortized time)
      }
      bool isEqual = true; // if all elements of q in the set A

      return isEqual;
    }

    /**
     * Evaluation can be done in O(n) in linear factors form.
     */
    T eval(T x) {
      T retVal = this->getConstantFactor();
      for (int i = 0; i < this->getDegree(); i++) {
        retVal *= (x - this->getFactors()[i]);
      }
      return retVal;
    }

    void toString(int precision) {
      std::cout << std::fixed << std::setprecision(precision)
                << this->getConstantFactor();
      int i = 0;
      while (i < this->getDegree()) {
        if (this->getFactors()[i] > 0) {
          std::cout << std::fixed << std::setprecision(precision)
                    << this->getFactors()[i] << "(x - " << this->getFactors()[i]
                    << ")";
        }
        i++;
      }
      std::cout << std::flush;
    }

    T* getFactors() { return this->x; }
    T getConstantFactor() { return this->b; }

  private:
    /* product of linear factors */
    T b = 0, *x = nullptr;
  };

  /**
   * Point-value representation of a polynomial.
   * Any polynomial of degree n is uniquely defined by
   * n+1 pairs (x_i,p(x_i)), where i = 0,..,n and x_i != x_j for i != j
   */
  template<typename T>
  class PV : public virtual Polynomial {
  public:
    /**
     * Constructor for point-value representation
     * Any polynomial of degree n is uniquely defined by
     * n+1 pairs (x_i,p(x_i)), where i = 0,..,n and x_i != x_j for i != j
     *
     * @param deg: degree of the polynomial (n)
     * @param numOfPoints: number of point-value pairs (should be at least n+1)
     * @param v: array of x_i values        [numOfPoints]
     * @param y: array of p(x_i) values     [numOfPoints]
     */
    PV(int deg, int numOfPoints, T v[], T y[]) {
      if (numOfPoints < deg + 1) {
        std::cout << "ERROR: Not enough point-value pairs! At least " << deg+1
                  << " point-value pairs required to represent a polynomial of degree " << deg
                  << std::endl;
        std::exit(1);
      }
      this->numOfPoints = numOfPoints;
      this->n = deg;
      this->v = new T[numOfPoints];
      this->y = new T[numOfPoints];
      for (int i = 0; i < numOfPoints; i++) {
        this->v[i] = v[i];
        this->y[i] = y[i];
      }
    }
    ~PV() {
      delete[] this->v;
      delete[] this->y;
    }

    /**
     *  Multiplication is in O(n) if degrees are same (same x values)
     *  multiply the y values.
     *
     *  condition: same degree, same points
     */
    PV<T>& operator *(const PV<T> &p) {
      if (this->getDegree() != p.getDegree()) {
        std::cout << "ERROR (PV): different degree!" << std::endl;
        std::exit(1);
      }
      int minNumOfPoints = std::min(this->getNumOfPoints(), p.getNumOfPoints());
      T values_prod = new T[minNumOfPoints];
      T points = new T[minNumOfPoints];
      bool error = false;
      /* O(n) */
      for (int i = 0; i < minNumOfPoints; i++) {
        T p_tmp = p->getX()[i];
        T q_tmp = this->getX()[i];
        error = error || (p_tmp != q_tmp);
        points[i] = p_tmp;
        values_prod[i] = p.getY()[i] + this->getY()[i];
      }
      PV<T> *retVal = new PV<T>(this->getDegree(), minNumOfPoints, points, values_prod);
      return *retVal;
    }

    /**
     * Multiplication with constant is O(n), multiply the y values
     */
    PV<T>& operator *(const T k) {
      T points = new T[this->numOfPoints];
      T values_prod = new T[this->numOfPoints];
      /* O(n) */
      for (int i = 0; i < this->numOfPoints; i++) {
        points[i] = this->getX()[i];
        values_prod[i] = this->getY()[i] * k;
      }
      PV<T> *retVal = new PV<T>(this->getDegree(), this->numOfPoints, points, values_prod);
      return *retVal;
    }

    /**
     * Same as multiplication with -1.
     */
    PV<T>& operator -() {
      return *this * -1;
    }

    /**
     *  Addition with another polynomial is in O(n) if degrees are same (same x values)
     *  add the y values.
     *
     *  Condition: same degree, same points
     */
    PV<T>& operator +(const PV<T> &p) {
      if (this->getDegree() != p.getDegree()) {
        std::cout << "ERROR (PV): different degree!" << std::endl;
        std::exit(1);
      }
      int minNumOfPoints = std::min(this->getNumOfPoints(), p.getNumOfPoints());
      T values_sum = new T[minNumOfPoints];
      T points = new T[minNumOfPoints];
      bool error = false;
      /* O(n) */
      for (int i = 0; i < minNumOfPoints; i++) {
        T p_tmp = p->getX()[i];
        T q_tmp = this->getX()[i];
        error = error || (p_tmp != q_tmp);
        points[i] = p_tmp;
        values_sum[i] = p.getY()[i] + this->getY()[i];
      }
      PV<T> *retVal = new PV<T>(this->getDegree(), minNumOfPoints, points, values_sum);
      return *retVal;
    }

    /**
     * Same as multiplying with constant + addition
     */
    PV& operator -(const PV<T> &p) {
      return (*this + (-p));
    }

    /**
     * Compares degree and values, returns true if polynomials are equal.
     *
     * Condition: same x values
     */
    bool operator ==(const PV<T> &p) {
      bool isEqual = true;
      isEqual = isEqual && (p.getDegree() == this->getDegree());
      int minNumPoints = std::min(this->getNumOfPoints(), p.getNumOfPoints());
      for (int i = 0; i < minNumPoints; i++) {
        isEqual = isEqual && (p.getX[i] == this->getX[i]) && (p.getY()[i] == this->getY()[i]);
      }
      return isEqual;
    }

    /* converts to coefficient representation */
    Coef<T>& toCoef() {
      return this->interpolate();
    }

    T* getX() { return this->v; };
    T* getY() { return this->y; };
    int getNumOfPoints() { return this->numOfPoints; }

  private:
    /* point-value representation */
    T *v = nullptr, *y = nullptr;
    int numOfPoints;

    /**
     * Interpolate using FFT in O(nlogn)
     * TODO
     */
    Coef<T>& interpolate() {
      // TODO
      Coef<T> *retVal;
      return *retVal;
    }
  };


  /**
   * Polynomial multiplication using FFT
   *
   * 1. Computes DFT_2n-2 of p = (p(w^0),p(w^1),..,p(w^2n-1))
   *    which converts coefficient representation to
   *    point-value representation using FFT, O(nlogn)
   * 2. Do multiplication in point-value representation, O(n)
   * 3. Interpolate the resulting polynomial to cofficient form
   *    using FFT, O(nlogn)
   *
   * @param p,q : polynomials in coefficient form
   * @return p*q : multiplication of p and q in coefficient form
   */
  template <typename U>
  Coef<U>& mult_fft(Coef<U> &p, Coef<U> &q) {
    int maxDegree = std::max(p.getDegree(), q.getDegree());
    int eval_number = 2*maxDegree - 2; // required num of evaluations
    /* compute the next power of 2 since fft will be called with a power of 2 */
    int nextPowTwo = std::pow(2, std::ceil(std::log2(eval_number)));
    /* convert p,q to point-value representations using FFT
     * so we get 2n point-value pairs
     *  y_p: values of poly p
     *  y_q: values of poly q
     *  x: corresponding x values
     */
    comp *y_p = new comp[nextPowTwo];
    comp *y_q = new comp[nextPowTwo];
    comp *x = new comp[nextPowTwo];
    y_p = p.fft(nextPowTwo);   // O(2nlog2n)
    y_q = p.fft(nextPowTwo);   // O(2nlog2n)
    comp w_2n = getRootUnity(nextPowTwo);
    comp current = 1;
    /* O(2n) */
    for (int i = 0; i < nextPowTwo; i++) {
      x[0] = current;
      current *= w_2n;
    }
    PV<comp> *pp = new PV<comp>(p.getDegree(), nextPowTwo, x, y_p);
    PV<comp> *qq = new PV<comp>(q.getDegree(), nextPowTwo, x, y_q);
    /* pointwise multiplication in O(2n) */
    PV<comp> res = *pp * *qq;
    /* interpolation via FFT in O(2nlog2n) */
    Coef<U> &retVal = res.toCoef();
    return retVal;
  }
} // namespace Polynomial

#endif // POLY_H
