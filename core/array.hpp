#ifndef ARRAY_H
#define ARRAY_H

#include <cstddef>
#include <algorithm>


namespace core {

namespace operation {

template<typename T>
void compound_add(T& a, const T& b) { a += b; }

template<typename T>
void compound_substract(T& a, const T& b) { a -= b; }

template<typename T>
void compound_multiplie(T& a, const T& b) { a *= b; }

template<typename T>
void compound_divide(T& a, const T& b) { a /= b; }


template<typename T>
T add(const T& a, const T& b) { return a + b; }

template<typename T>
T substract(const T& a, const T& b) { return a - b; }

template<typename T>
T multiplie(const T& a, const T& b) { return a * b; }

template<typename T>
T divide(const T& a, const T& b) { return a / b; }

}

template<typename T>
class array {
 public:
  typedef T value_type;
  
  array(): size(0), elements(NULL) {}
  
  array(const array<T>& op): size(op.size), elements(NULL) {
    if (op.elements) {
      elements = new T[size];
      std::copy(op.elements, op.elements + size, elements);
    }
  }

  array(array<T>&& op): size(0), elements(nullptr) {
    std::swap(elements, op.elements);
    std::swap(size, op.size);
  }

  explicit array(std::size_t _size)
      : size(_size), elements(new T[size]) {}

  virtual ~array() {
    clear();
  }

  void clear() {
    delete [] elements;
    elements = NULL;
    size = 0;
  }

  const array<T>& operator=(const array<T>& op) {
    clear();
    
    if (op.elements) {
      size = op.size;
      elements = new T[size];

      std::copy(op.elements, op.elements + size, elements);
    }

    return *this;
  }

  array<T>& operator=(array<T>&& op) {
    clear();
    
    std::swap(elements, op.elements);
    std::swap(size, op.size);

    return *this;
  }


  std::size_t get_size() const { return size; }
  
  T& operator()(std::size_t i) {
    return elements[i];
  }
  
  const T& operator()(std::size_t i) const {
    return elements[i];
  }

  array<T>& operator+=(const array<T>& op) {
    apply_compound_vector_operation(operation::compound_add<T>, op);
    return *this;
  }

  array<T>& operator+=(const T& s) {
    apply_compound_scalar_operation(operation::compound_add<T>, s);
    return *this;
  }
  
  array<T>& operator-=(const array<T>& op) {
    apply_compound_vector_operation(operation::compound_substract<T>, op);
    return *this;
  }

  array<T>& operator-=(const T& s) {
    apply_compound_scalar_operation(operation::compound_substract<T>, s);
    return *this;
  }
  
  array<T>& operator*=(const array<T>& op) {
    apply_compound_vector_operation(operation::compound_multiplie<T>, op);
    return *this;
  }

  array<T>& operator*=(const T& s) {
    apply_compound_scalar_operation(operation::compound_multiplie<T>, s);
    return *this;
  }
  
  array<T>& operator/=(const array<T>& op) {
    apply_compound_vector_operation(operation::compound_divide<T>, op);
    return *this;
  }

  array<T>& operator/=(const T& s) {
    apply_compound_scalar_operation(operation::compound_divide<T>, s);
    return *this;
  }

  array<T> operator-() const {
    array<T> result(size);
    for (std::size_t i(0); i < size; ++i)
      result.elements[i] = -elements[i];
    return result;
  }

  array<T>& invert() {
    for (std::size_t i(0); i < size; ++i)
      elements[i] = 1.0 / elements[i];
    return *this;
  }
  
  template<typename operation_type>
  void apply_compound_vector_operation(operation_type op, const array<T>& a) {
    for (std::size_t i(0); i < size; ++i)
      op(elements[i], a.elements[i]);
  }

  template<typename operation_type>
  void apply_compound_scalar_operation(operation_type op, const T& a) {
    for (std::size_t i(0); i < size; ++i)
      op(elements[i], a);
  }
  
 private:
  std::size_t size;
  T* elements;
};

}  // namespace core

#endif /* ARRAY_H */
