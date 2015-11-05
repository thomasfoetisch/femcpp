#ifndef MATH_H
#define MATH_H

namespace core {

template<typename T>
inline T pow(T b, unsigned int e) {
  return e == 0 ? 1.0 : b * pow(b, e - 1);
}

inline
std::size_t factorial(std::size_t n) {
  return n == 0 ? 1 : factorial(n - 1);
}

inline
std::size_t nchoosek(std::size_t n, std::size_t k) {
  return factorial(n) / (factorial(k) * factorial(n - k));
}

template<typename T>
T determinant(const ndarray<T>& a);

template<typename T>
ndarray<T> invert_matrix(const ndarray<T>& a);

}

#endif /* MATH_H */
