#include "point2D.hpp"
#include <math.h>

Point2D::Point2D(float x, float y) {
  this->x = x;
  this->y = y;
}

Point2D Point2D::operator +(const Point2D p) {
  return Point2D(this->x + p.x, this->y + p.y);
}

Point2D Point2D::operator -(const Point2D p) {
  return Point2D(this->x - p.x, this->y - y);
}

void _mergeX(Point2D s[], int a, int b) {
  // TODO Implement merging step for X
}

void _mergeY(Point2D s[], int a, int b) {
  // TODO Implement merging step for Y
}

void mergesort(Point2D s[], int a, int b, char order) {
  // TODO implement merge sort
}

void Point2D::sort_orderX(Point2D s[], int n) {
  mergesort(s, 0, n+1, 'X');
}

void Point2D::sort_orderY(Point2D s[], int n) {
  mergesort(s, 0, n+1, 'Y');
}

float Point2D::dist2D_L1(Point2D p) {
  return abs(this->x - p.x) + abs(this->y - p.y);
}

float Point2D::dist2D_L2(Point2D p) {
  return sqrtf(pow(this->x - p.x, 2) + pow(this->y - p.y, 2));
}

float Point2D::getX() {
  return this->x;
}

float Point2D::getY() {
  return this->y;
}

void Point2D::setX(float x) {
  this->x = x;
}

void Point2D::setY(float y) {
  this->y = y;
}
