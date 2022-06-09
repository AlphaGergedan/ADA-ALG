#ifndef POINT2D_HPP
#define POINT2D_HPP

#include <iostream>

class Point2D {
public:
  Point2D(float x = 0, float y = 0);
  Point2D operator +(const Point2D p);
  Point2D operator -(const Point2D p);

  /* L1 distance to another point */
  float dist2D_L1(Point2D p);
  /* L2 distance to another point */
  float dist2D_L2(Point2D p);

  float getX();
  float getY();
  void setX(float x);
  void setY(float y);

  void toString();
private:
  float x, y;
}; // end of class Point2D

/* Sort a list of 2d points in increasing x or y coordinate */
void sort_orderX(Point2D *s[], int n);
void sort_orderY(Point2D *s[], int n);

#endif
