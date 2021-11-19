#ifndef POINT2D_HPP
#define POINT2D_HPP

#include <math.h>

class Point2D {
public:
  Point2D(float x = 0, float y = 0);

  Point2D operator +(const Point2D p);
  Point2D operator -(const Point2D p);

  float dist2D_L1(Point2D p1, Point2D p2);
  float dist2D_L2(Point2D p1, Point2D p2);

  void sort_orderX(Point2D s[], int n);
  void sort_orderY(Point2D s[], int n);

  float getX();
  float getY();
  void setX(float x);
  void setY(float y);
private:
  float x, y;
}; // end of class Point2D

#endif
