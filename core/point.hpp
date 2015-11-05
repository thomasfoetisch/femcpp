#ifndef POINT_H
#define POINT_H

namespace core {

typedef ndarray<double> point;

point make_point(double x, double y) {
  point p({2});
  
  p(0) = x;
  p(1) = y;

  return p;
}

point make_point(double x, double y, double z) {
  point p({3});
  
  p(0) = x;
  p(1) = y;
  p(2) = z;

  return p;
}

}  // namespace core

#endif /* POINT_H */
