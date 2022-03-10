#include "closestPair.hpp"

std::pair<int,int> findClosestPoints_L2_NAIVE(Point2D *s[], int n) {
  /* indices of the smallest distance */
  int l_i, l_j;
  /* the value of the smallest distance */
  float l_d = INFINITY;
  /* For each pair compute the distance */
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j) {
        float d = s[i]->dist2D_L2(*s[j]);
        if (d < l_d) {
          l_i = i;
          l_j = j;
          l_d = d;
        }
      }
    }
  }
  return std::make_pair(l_i, l_d);
}

int main() {
  /* Read user input */
  // Create n 2D-points in range (a,b)
  int n = 100;
  float a = -1000,b = 1000;

  Point2D *s[n];
  srand (static_cast <unsigned> (time(0)));
  for (int i = 0; i < n; i++) {
    /* generate random points in range [low,high] */
    float x = a + ( (static_cast<float>(rand())) / static_cast<float>(RAND_MAX) ) * (b - a);
    float y = a + ( (static_cast<float>(rand())) / static_cast<float>(RAND_MAX) ) * (b - a);
    s[i] = new Point2D(x,y);
    std::cout << "[" << i << "] <- "; s[i]->toString(); std::cout << std::endl;
  }

  std::pair<int,int> p = findClosestPoints_L2_NAIVE(s, n);
  std::cout << "Closest points are "; s[p.first]->toString(); std::cout << " and "; s[p.second]->toString(); std::cout << std::endl;
}
