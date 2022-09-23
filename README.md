# Advanced Algorithms

Implementations of the following advanced algorithms.

## Divide and Conquer Algorithms

### Sorting Problem

Given a list of integers return the sorted list. Implementation:
- Mergesort (TODO)
- Quicksort: Sorting by partitioning `./quicksort`
- Before you use, build it with `bazel build //src/quicksort:quicksort`
- `bazel build //src/quicksort:quicksort-run` for building an example execution.

### Closest Pair Problem

Given a set *S* of *n* points in the plane, find a pair of points with the
smallest distance. Implementation:
- (TODO)

### Line Segment Intersection Problem

### Computation of a Voronoi Diagram

## Fast Fourier Transform

### Polynomials

Operations on polynomials
- Addition
- Multiplication
- Evaluation

Representation of Polynomials
- Coefficient Representation
- Product of Linear Factors
- Point-Value Representation

Fast Fourier Transform (FFT)

- `bazel build //src/fft:fft` (header only lib)

## Randomized Algortihms

### Randomized Quicksort

### Primality Test

Given a natural number greater equal two, tell if it is prime or not. Implementation:
- Deterministic Primality Test (Naive): (TODO)
- Simple Primality Test: (TODO)
- Randomized Primality Test: (TODO)
- Miller Rabin Primality Test: (TODO)

### Verifying Matrix Mulitplication

Application in cryptography

## Algorithms for Treaps

Operations on Treaps:
- Search for element with key k
- Insert a new element x
- Delete an element x

## Minimum-Cut Algorithms

## Algorithms for SuffixTrees

Construction of suffix trees
- naive implementation O(n^2)
- The algorithm MCC O(n)

## Fibonacci Heap

Cut when second child is lost, or when k-th child is lost.

## Dynamic Programming

Fibonacci

weighted interval partitioning

matrix multiplication order

longest common subsequence of two strings

monochromatic sequence

Subset sum / Knapsack


