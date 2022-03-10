#include "point2D.hpp"
#include <cmath>

Point2D::Point2D(float x, float y) {
  this->x = x;
  this->y = y;
}

Point2D Point2D::operator +(const Point2D p) {
  return Point2D(this->x + p.x, this->y + p.y);
}

Point2D Point2D::operator -(const Point2D p) {
  return Point2D(this->x - p.x, this->y - p.y);
}

void swap(Point2D *s[], int i, int j) {
  Point2D *tmp = s[i];
  s[i] = s[j];
  s[j] = tmp;
}

int partitionX(Point2D *s[], int a, int b, int splitter_index) {
    /* Swap the splitter to the end of the list */
    swap(s, splitter_index, b);
    /* Index of the partition */
    int i = a - 1;
    for (int j = a; j < b; j++) {
        if (s[j]->getX() < s[b]->getX()) {
            swap(s, ++i, j);
        }
    }
    /* Place the splitter element into the correct position */
    swap(s, ++i, b);
    return i;
}

int partitionY(Point2D *s[], int a, int b, int splitter_index) {
    /* Swap the splitter to the end of the list */
    swap(s, splitter_index, b);
    /* Index of the partition */
    int i = a - 1;
    for (int j = a; j < b; j++) {
        if (s[j]->getY() < s[b]->getY()) {
            swap(s, ++i, j);
        }
    }
    /* Place the splitter element into the correct position */
    swap(s, ++i, b);
    return i;
}

void quicksortX(Point2D *s[], int a, int b) {
  if (a < b) {
    int splitter_index = b;
    int splitter_position = partitionX(s, a, b, splitter_index);

    quicksortX(s, a, splitter_position - 1);
    quicksortX(s, splitter_position + 1, b);
  }
}

void quicksortY(Point2D *s[], int a, int b) {
  if (a < b) {
    int splitter_index = b;
    int splitter_position = partitionY(s, a, b, splitter_index);

    quicksortY(s, a, splitter_position - 1);
    quicksortY(s, splitter_position + 1, b);
  }
}

/* Sorts the given array in order of increasing X-coordinate */
void sort_orderX(Point2D *s[], int n) {
  quicksortX(s, 0, n-1);
}

/* Sorts the given array in order of increasing Y-coordinate */
void sort_orderY(Point2D *s[], int n) {
  quicksortY(s, 0, n-1);
}

float Point2D::dist2D_L1(Point2D p) {
    return fabs(this->x - p.x) + fabs(this->y - p.y);
}

float Point2D::dist2D_L2(Point2D p) {
    return sqrtf(powf(this->x - p.x, 2) + powf(this->y - p.y, 2));
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

void Point2D::toString() {
  std::cout << "(" << this->x << "," << this->y << ")" << std::flush;
}
