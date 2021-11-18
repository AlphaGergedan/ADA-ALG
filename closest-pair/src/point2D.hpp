#ifndef POINT2D_HPP
#define POINT2D_HPP

#include <math.h>

class Point2D {
public:
  Point2D(float x = 0, float y = 0);

  Point2D operator +(const Point2D p);
  Point2D operator -(const Point2D p);

  void sort_orderX(Point2D s[], int n);
  void sort_orderY(Point2D s[], int n);

  float getX();
  float getY();
  float setX();
  float setY();

private:
  float x, y;

  float dist2D_L2(Point2D p1, Point2D p2);
  float dist2D_L1(Point2D p1, Point2D p2);

  void mergesort(Point2D s[], int a, int b, char order);

  void _mergeX(Point2D s[], int a, int b);
  void _mergeY(Point2D s[], int a, int b);

}; // end of class Point2D

#endif
