#ifndef CLOSEST_HPP
#define CLOSEST_HPP

#include <tuple>
#include <random>
#include <time.h>

#include "point2D.hpp"

/**
 * Finds the closest pair of points, returns their indices
 */
std::pair<int,int> findClosestPoints_L2_NAIVE(Point2D *s[], int n);

#endif
