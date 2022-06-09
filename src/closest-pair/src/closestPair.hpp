#ifndef CLOSEST_HPP
#define CLOSEST_HPP

#include <tuple>
#include <random>
#include <time.h>
#include <iomanip>

#include "point2D.hpp"

/**
 * Finds the closest pair of points by brute force comparison in O(n^2) (compares all pairs).
 * Returns their indices with the distance value.
 *
 * @param s Array of 2d-points
 * @param a Lower index (included)
 * @param b Upper index (included)
 * @return (p1,p2,mindist)
 */
std::tuple<Point2D*,Point2D*,float> findClosestPoints_L2_NAIVE(Point2D *s[], int a, int b);

/**
 * Finds the closest pair of points in O(nlogn) time using D&C.
 * Requires preprocessing: sort in order of increasing x-coordinates O(nlogn).
 *
 * @param s Array of 2d-points
 * @param a Lower index (included)
 * @param b Upper index (included)
 * @return (p1,p2,mindist)
 */
std::tuple<Point2D*,Point2D*,float> findClosestPoints_L2(Point2D *s[], int a, int b);

#endif
