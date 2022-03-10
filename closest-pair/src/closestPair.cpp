#include "closestPair.hpp"

/*
 * Compare all pairs in O(n^2)
 */
std::tuple<Point2D*,Point2D*,float> findClosestPoints_L2_NAIVE(Point2D *s[], int a, int b) {
  if (a == b) {
    return std::make_tuple(s[a],s[a],INFINITY);
  }
  /* indices of the smallest distance */
  int l_i, l_j;
  /* the value of the smallest distance */
  float l_d = INFINITY;
  /* For each pair compute the distance */
  for (int i = a; i < b; i++) {
    for (int j = i + 1; j <= b; j++) {
      /* never occurs TODO */
      if (i != j) {
        float d = s[i]->dist2D_L2(*s[j]);
        if (d < l_d) {
          l_i = i;
          l_j = j;
          l_d = d;
        }
      } else {
        // FIXME
        exit(1);
      }
    }
  }
  return std::make_tuple(s[l_i], s[l_j], l_d);
}


/*
 * Finds closest pair of points in O(nlogn) time using D&C
 */
std::tuple<Point2D*,Point2D*,float> findClosestPoints_L2(Point2D *s[], int a, int b) {
  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "CALL: (" << a << "," << b << ") BEGIN" << std::endl;
  /* size of the current subproblem */
  int n = b - a + 1;
  /* indices of the smallest distance */
  Point2D *p1, *p2;
  /* the value of the smallest distance */
  float mindist = INFINITY;

  /* determine bisection line O(1) */
  float bi_line = (s[a]->getX() + s[b]->getX()) / 2;
  std::cout << "--> bisection line at " << bi_line << " (" << a << ", " << b << ")" << std::endl;

  /* divide s into two equal sets s_left and s_right, divided by the bisection line */
  int i = a;
  while (i <= b) {
    if (s[i]->getX() > bi_line) {
      break;
    }
    i++;
  }
  i--;
  /* size of the subproblems in left and right */
  int n_l = i - a + 1, n_r = b - (i+1) + 1;

  std::cout << "subproblem sizes = " << n_l << " (" << a << "," << i << ") and " << n_r
            << " (" << i+1 << "," << b << ")" << std::endl;

  /* if there are less than 3 points solve brute force */
  std::tuple<Point2D*,Point2D*,float> l, r;
  if (n_l < 4) {
    l = findClosestPoints_L2_NAIVE(s, a, i);
    std::cout << "d_l = " << std::get<2>(l) << ", subproblem (" << a << "," << i << ")" << std::endl;
  } else {
    l = findClosestPoints_L2(s, a, i);
  }
  if (n_r < 4) {
    r = findClosestPoints_L2_NAIVE(s, i+1, b);
    std::cout << "d_r = " << std::get<2>(r) << ", subproblem (" << i+1 << "," << b << ")" << std::endl;
  } else {
    r = findClosestPoints_L2(s, i+1, b);
  }

  /* consider points within distance d to the bisection line. For each point p consider all points
   * q within y-distance at most d. There are at most 7 such points in L2 distance. */
  float d = std::fmin(std::get<2>(l), std::get<2>(r));
  std::cout << "===============d_min left right is = " << d << std::endl;
  int low_i = a;
  /* determine the lower index of the array */
  while (low_i <= b) {
    if ((s[low_i]->getX() >= (bi_line - d)) && (s[low_i]->getX() <= (bi_line + d))) {
      break;
    }
    low_i++;
  }
  std::cout << "--> low_i = " << low_i << ", " << s[low_i]->getX() << " (" << a << ", " << b << ")" << std::endl;
  int high_i = low_i;
  if (low_i <= b) {
    while (high_i <= b) {
      if (s[high_i]->getX() > (bi_line + d)) {
        break;
      }
      high_i++;
    }
    high_i--;
  } else {
    high_i = b; // there are no such points
  }
  std::cout << "--> high_i = " << high_i << ", " << s[high_i]->getX() << " (" << a << ", " << b << ")" << std::endl;
  std::cout << "--> interval from " << bi_line - d << " to " << bi_line + d << " (" << a << ", " << b << ")" << std::endl;

  /* copy those points into another array that we can sort in order y-coordinates O(n) */
  int n_x = high_i - low_i + 1;
  std::cout << "n_x = " << n_x << " (" << a << ", " << b << ")" << std::endl;
  Point2D *x[n_x];
  for (int i = low_i; i <= high_i; i++) {
    x[i] = s[i];
  }
  sort_orderY(x, n_x);

  /* tuple for smallest distance between points within the bisection line */
  std::tuple<Point2D*,Point2D*,float> m = {NULL,NULL,INFINITY};
  /* computation in O(n) */
  for (int i = 0; i < n_x; i++) {
    int checkedPoints = 1;
    while (checkedPoints <= 7 && (i+checkedPoints) < n_x) {
      float d = s[i]->dist2D_L2(*s[i + checkedPoints]);
      if (d < std::get<2>(m)) {
        m = {x[i],x[i+checkedPoints],d};
      }
      checkedPoints++;
    }
  }

  float mindist_l = std::get<2>(l), mindist_r = std::get<2>(r), mindist_m = std::get<2>(m);
  std::cout << "mindist_l = " << mindist_l << ", mindist_r = " << mindist_r << ", mindist_m = " << mindist_m << ", subprob. (" << a << ", " << b << ")" << std::endl;

  if (mindist_l < mindist_r && mindist_l < mindist_m) {
    p1 = std::get<0>(l);
    p2 = std::get<1>(l);
    mindist = mindist_l;
  } else if (mindist_r < mindist_l && mindist_r < mindist_m) {
    p1 = std::get<0>(r);
    p2 = std::get<1>(r);
    mindist = mindist_r;
  } else if (mindist_m < mindist_l && mindist_m < mindist_r) {
    p1 = std::get<0>(m);
    p2 = std::get<1>(m);
    mindist = mindist_m;
  } else {
    p1 = std::get<0>(l);
    p2 = std::get<1>(l);
    mindist = mindist_l;
  }

  std::cout << "SO we return " << mindist << " for subprob. (" << a << ", " << b << ")" << std::endl;

  return std::make_tuple(p1, p2, mindist);
}

int main() {
  /* Read user input */
  // Create n 2D-points in range (a,b)
  int n = 10;
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

  std::tuple<Point2D*,Point2D*,float> p = findClosestPoints_L2_NAIVE(s, 0, n-1);
  std::cout << "Closest points are "; std::get<0>(p)->toString();
  std::cout << " and "; std::get<1>(p)->toString();
  std::cout << " with distance " << std::setprecision(7) << std::get<2>(p) << std::endl;

  /* D&C approach in O(nlogn) */
  sort_orderX(s, n); // preprocessing
  std::cout << "preprocessing done \u2713" << std::endl;
  p = findClosestPoints_L2(s, 0, n-1);
  std::cout << "Closest points are "; std::get<0>(p)->toString();
  std::cout << " and "; std::get<1>(p)->toString();
  std::cout << " with distance " << std::setprecision(7) << std::get<2>(p) << std::endl;
}
