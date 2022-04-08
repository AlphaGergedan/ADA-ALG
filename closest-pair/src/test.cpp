#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>
#include "point2D.hpp"

/* Print a list of points */
void toString(Point2D *s[], int size) {
  for (int i = 0; i < size; i++) {
    std::cout << "[" << i << "] <- "; s[i]->toString(); std::cout << std::endl;
  }
}


Point2D points[20] = {Point2D(-9.23, 6.3),Point2D(-8.4331, 4.332),Point2D(-5.898, 7.3),Point2D(-5.5, -2),
                      Point2D(-5, -5),Point2D(-5, 8),Point2D(-4.8, 10),Point2D(-4.53, -0.232),
                      Point2D(-9.23, 6.3),Point2D(-9.23, 6.3),Point2D(-9.23, 6.3),Point2D(-9.23, 6.3),
                      Point2D(-9.23, 6.3),Point2D(-9.23, 6.3),Point2D(-9.23, 6.3),Point2D(-9.23, 6.3),
                      Point2D(-9.23, 6.3),Point2D(-9.23, 6.3),Point2D(-9.23, 6.3),Point2D(-9.23, 6.3)};


void testSum() {
    Point2D p1 = Point2D(312,412);
    Point2D p2 = Point2D(32,69);

    Point2D p_sum = p1 + p2;
    Point2D p_sub = p1 - p2;

    std::cout << "p_sum = (" << p_sum.getX() << "," << p_sum.getY() << ")"
              << ", p_sub = (" << p_sub.getX() << "," << p_sub.getY() << ")" << std::endl;

    std::cout << "Should be same as " << "(" << p1.getX() + p2.getX() << "," << p1.getY() + p2.getY() << "), ("
              << p1.getX() - p2.getX() << "," << p1.getY() - p2.getY() << ")" <<std::endl;
}

void test_L1_distance() {
    Point2D p1 = Point2D(312,412);
    Point2D p2 = Point2D(32,69);
    /* L1 distance should be 623 */
    float d1 = p1.dist2D_L1(p2);
    float d2 = p2.dist2D_L1(p1);
    float d3 = p1.dist2D_L1(p1);
    std::cout << d1 << " == " << "623" << " == " << d2 << std::endl;
    std::cout << d3 << " == " << "0";
}

void test_L2_distance() {
    Point2D p1 = Point2D(312,412);
    Point2D p2 = Point2D(32,69);
    /* L2 distance should be 442.7742088 */
    float d1 = p1.dist2D_L2(p2);
    float d2 = p2.dist2D_L2(p1);
    float d3 = p1.dist2D_L2(p1);
    std::cout << std::fixed << std::setprecision(8) << d1 << " == " << "442.7742088" << " == " << d2 << std::endl;
    std::cout << std::fixed << std::setprecision(8) << d3 << " == " << "0" << std::endl ;
}

int main() {
    testSum();
    test_L1_distance();
    test_L2_distance();

    Point2D *s[20];

    /* random point sampling */
    srand (static_cast <unsigned> (time(0)));
    float low = -10, high = 10;
    int size = sizeof(points) / sizeof(points[0]);
    for (int i = 0; i < size; i++) {
        /* generate random points in range [low,high] */
        float x = low + ( (static_cast<float>(rand())) / static_cast<float>(RAND_MAX) ) * (high - low);
        float y = low + ( (static_cast<float>(rand())) / static_cast<float>(RAND_MAX) ) * (high - low);
        s[i] = new Point2D(x,y);
        std::cout << "[" << i << "] <- "; s[i]->toString(); std::cout << std::endl;
    }

    sort_orderX(s, size);
    std::cout << "-- ORDERED BY INCREASING X-COORDINATE -- " << std::endl;
    toString(s, size);
    sort_orderY(s, size);
    std::cout << "-- ORDERED BY INCREASING Y-COORDINATE -- " << std::endl;
    toString(s, size);
}
