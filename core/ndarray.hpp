#ifndef NDARRAY_H
#define NDARRAY_H

#include <vector>
#include <algorithm>
#include <functional>

#include <meta.hpp>
#include "array.hpp"
#include "multiindex.hpp"
#include "slice.hpp"


namespace core {

template<typename T>
class data_accessor_base {
 public:
  static void build_offsets(const multiarray& dimensions) {
    if (dimensions.size()) {
      offsets.resize(dimensions.size());
      offsets.back() = 1;

      std::size_t partial_product(dimensions.back());
      for (std::size_t dim(1); dim < dimensions.size(); ++dim) {
        offsets[dimensions.size() - 1 - dim] = partial_product;
        partial_product *= dimensions[dimensions.size() - 1 - dim];
      }
    }
  }
  
  static std::size_t index_to_offset(const multiindex& m,
                                     const multiindex& offsets) const {
    std::size_t result(0);
    const std::size_t ignore(m.size() - get_rank());
    
    for (std::size_t i(0); i < get_rank(); ++i)
      result += m[ignore + i] * offsets[i];
    
    return result;
  }
  
  static std::size_t get_length_from_dimensions(const multiindex& dims) const {
    if (dims.size() == 0)
      return 0;

    return std::accumulate(dims.begin(), dims.end(),
                           1, std::multiplies<std::size_t>());
  }
};

template<typename T>
class default_data_accessor: public data_accessor_base {
 public:
  const T& operator()() const;
  T& operator()();
  
 private:
  array<T> elements;
  multiindex dimensions, offsets;
};

template<typename T>
class slice_data_accessor: public data_accessor_base {
 public:
  const T& operator()() const;
  T& operator()();
  
 private:
  multiindex dimensions, offsets, index_base, strides, extends;
};


/*
 *  Multidimensional array representation
 */
template<typename T, typename data_accessor = default_data_accessor<T> >
class ndarray : public data_accessor {
 public:
  typedef T value_type;
  
  /*
   *  Default construction
   */
  ndarray(): dimensions(0), offsets(0), elements(0) {
    build_offsets();
  }

  /*
   *  Copy construction
   */
  explicit ndarray(const multiindex& _dimensions)
      : dimensions(_dimensions),
        offsets(_dimensions.size(), 0),
        elements(get_length_from_dimensions(_dimensions)) {
    build_offsets();
  }

  /*
   *  Slice construction
   */
  ndarray(const ndarray<T>& a, const std::vector<slice>& ss)
      : dimensions(0),
        offsets(ss.size(), 0),
        elements() {
    multiindex start, stride;
    for (auto s: ss) {
      dimensions.push_back(s.get_size());
      start.push_back(s.get_start());
      stride.push_back(s.get_stride());
    }

    elements = array<T>(get_length_from_dimensions(dimensions));
    build_offsets();
    
    multiindex indices(dimensions.size(), 0);
    do {
      (*this)(indices) = a(start + stride * indices);
    } while (increment(indices, dimensions));
  }


  /*
   *  Default destruction
   */
  ~ndarray() = default;

  
  /*
   *  Multidimensional indexing helper structures
   */
  template<typename element_ref, typename ... args>
  struct slicing_return_type {
    static constexpr bool is_index =
        meta::is_integral_typelist<meta::typelist<args...> >::value;
    static constexpr bool is_slice = not is_index;
  
    typedef
    typename meta::if_<is_slice,
                       core::ndarray<T>,
                       element_ref>::type
    type;
  };


  /*
   *  Factory for the indexing operation return value
   */
  template<bool is_slice, typename ... args>
  struct index_impl;

  /*
   *  Factory specialization for the slicing case
   */
  template<typename ... args>
  struct index_impl<true, args...> {
    
    template<typename head_type, typename ... tail_type>
    struct build_slice_data {
      static void do_it(std::vector<slice>::iterator s,
                        head_type head_value, tail_type ... tail_values) {
        slice_from_index<head_type>(s, head_value);
        build_slice_data<tail_type...>::do_it(++s, tail_values...);
      }
    };

    template<typename tail_type>
    struct build_slice_data<tail_type> {
      static void do_it(std::vector<slice>::iterator s,
                        tail_type tail_value) {
        slice_from_index<tail_type>(s, tail_value);
      }
    };

    /*
     *  Factory method
     */
    static
    core::ndarray<T>
    do_it(const core::ndarray<T>& elements,
          args ... arg_values) {
      std::vector<slice> slices(sizeof...(args));
      build_slice_data<args...>::do_it(slices.begin(),
                                       arg_values...);

      return core::ndarray<T>(elements, slices);
    }
  };

  /*
   *  Factory specialisation for the pure indexing case
   */
  template<typename ... args>
  struct index_impl<false, args...> {
    
    template<typename head_arg, typename ... tail_args>
    struct build_offset {
      static
      void do_it(std::size_t& offset,
                 core::multiindex::const_iterator dim,
                 head_arg i, tail_args ... arg_values) {
        offset += (*dim) * i;
        build_offset<tail_args...>::do_it(offset, ++dim, arg_values...);
      }
    };

    template<typename tail_arg>
    struct build_offset<tail_arg> {
      static
      void do_it(std::size_t& offset,
                 core::multiindex::const_iterator dim,
                 tail_arg i) {
        offset += i;
      }
    };

    /*
     *  Factory method, non const version
     */
    static
    typename slicing_return_type<T&, args...>::type
    do_it(core::ndarray<T>& elements,
          args ... arg_values) {
      std::size_t offset(0);
      build_offset<args...>::do_it(offset,
                                   elements.offsets.begin(),
                                   arg_values...);
    
      return elements.elements(offset);
    }

    /*
     *  Factory method, const version
     */
    static
    typename slicing_return_type<const T&, args...>::type
    do_it(const core::ndarray<T>& elements,
          args ... arg_values) {
      std::size_t offset(0);
      build_offset<args...>::do_it(offset,
                                   elements.offsets.begin(),
                                   arg_values...);
    
      return elements.elements(offset);
    }
  };


  /*
   *  Multidimensional indexing methods
   */
  template<typename ... args>
  auto index(args ... arg_values)
      -> typename slicing_return_type<T&, args...>::type {
    return index_impl<slicing_return_type<T&, args...>::is_slice,
                      args...>::do_it(*this, arg_values...);
  }

  template<typename ... args>
  auto index(args ... arg_values) const
      -> typename slicing_return_type<const T&, args...>::type {
    return index_impl<slicing_return_type<const T&, args...>::is_slice,
                      args...>::do_it(*this, arg_values...);
  }


  /*
   *  Standard indexing methods
   */
  T& operator()(std::size_t offset) {
    return elements(offset);
  }
 
  const T& operator()(std::size_t offset) const {
    return elements(offset);
  }

  T& operator()(const multiindex& m) {
    return elements(index_to_offset(m));
  }
  
  const T& operator()(const multiindex& m) const {
    return elements(index_to_offset(m));
  }


  /*
   *  Operator overloadings
   */
  ndarray<T, data_accessor>& operator=(const ndarray<T, data_accessor>& op) = default;
  ndarray<T>& operator=(const T& op) {
    elements = op;
    return *this;
  }
  
  ndarray<T>& operator+=(const ndarray<T>& op) {
    apply_compound_operation(core::operation::compound_add<T>, op);
    return *this;
  }
  ndarray<T>& operator+=(const T& op) {
    elements += op;
    return *this;
  }
  ndarray<T>& operator-=(const ndarray<T>& op) {
    apply_compound_operation(operation::compound_substract<T>, op);
    return *this;
  }
  ndarray<T>& operator-=(const T& op) {
    elements -= op;
    return *this;
  }

  ndarray<T> operator-() const {
    ndarray<T> result;
    result.elements = -elements;
    result.dimensions = dimensions;
    result.offsets = offsets;
    return result;
  }
  
  ndarray<T>& operator*=(const ndarray<T>& op) {
    apply_compound_operation(operation::compound_multiplie<T>, op);
    return *this;
  }
  ndarray<T>& operator*=(const T& op) {
    elements *= op;
    return *this;
  }

  ndarray<T>& operator/=(const ndarray<T>& op) {
    apply_compound_operation(operation::compound_divide<T>, op);
    return *this;
  }
  ndarray<T>& operator/=(const T& op) {
    elements /= op;
    return *this;
  }

  ndarray<T>& invert() {
    elements.invert();
    return *this;
  }
  
  ndarray<T>& transpose(std::size_t dim_1 = 0, std::size_t dim_2 = 1);
  
  /*
   *  Layout inspection methods
   */
  std::size_t get_rank() const { return dimensions.size(); }
  std::size_t get_dimension(std::size_t i) const { return dimensions[i]; }
  const multiindex& get_dimensions() const { return dimensions; }
  std::size_t get_length() const { return elements.get_size(); }


  /*
   *  Layout alteration methods
   */
  ndarray<T>& squeeze() {
    multiindex new_dimensions(std::count_if(dimensions.begin(),
                                            dimensions.end(),
                                            meta::not_equal<std::size_t>(1ul)),
                              0);
    
    std::copy_if(dimensions.begin(), dimensions.end(),
                 new_dimensions.begin(),
                 meta::not_equal<std::size_t>(1ul));

    dimensions = new_dimensions;
    build_offsets();

    return *this;
  }
  
  void reshape(const multiindex& m) {
    if (get_length_from_dimensions(m) != get_length())
      throw std::exception();
    dimensions = m;
    build_offsets();
  }


  /*
   *  Extern broadcasting operator helpers:
   */
  template<typename operation>
  static ndarray<T> broadcast_operation(operation op,
                                        const ndarray<T>& a,
                                        const ndarray<T>& b) {
    ndarray<T> result(broadcasted_dimensions(a.dimensions, b.dimensions));
    multiindex
        i_1(result.get_rank(), 0),
        i_2(result.get_rank(), 0),
        i_r(result.get_rank(), 0);
    multiindex
        d_1(normalize_index(a.dimensions, result.get_rank())),
        d_2(normalize_index(b.dimensions, result.get_rank()));
    do {
      result(i_r) = op(a(i_1), b(i_2));
    } while (increment(i_r, result.dimensions)
             and
             broadcast_increment(i_1, d_1,
                                 i_2, d_2));

    return result;
  }
  
  template<typename operation>
  static ndarray<T> apply_operation(operation op,
                                    const ndarray<T>& a,
                                    const ndarray<T>& b) {
    if (a.dimensions == b.dimensions) {
      ndarray<T> copy(a);
      copy.elements.apply_compound_vector_operation(op, b.elements);
      return copy;
    }

    if (is_broadcastable(a.dimensions, b.dimensions))
      return core::ndarray<T>::broadcast_operation(op, a, b);
    else
      throw std::exception();
  }

 private:

  /*
   *  Intern broadcasting operator helpers:
   */
  template<typename operation>
  void apply_compound_operation(operation op, const ndarray<T>& a) {
    if (dimensions.size() < a.dimensions.size())
      throw std::exception();

    if (dimensions == a.dimensions) {
      elements.apply_compound_vector_operation(op, a.elements);
      return;
    }

    if (is_compound_broadcastable(dimensions, a.dimensions)) {
      broadcast_compound_operation(op, a);
      return;
    } else {
      throw std::exception();
    }
  }

  template<typename operation>
  void broadcast_compound_operation(operation op, const ndarray<T>& a) {
    multiindex indices_1(get_rank(), 0), indices_2(get_rank(), 0);
    multiindex dim(normalize_index(a.dimensions, get_rank()));
    do {
      op(this->operator()(indices_1), a(indices_2));
    } while (broadcast_increment(indices_1, dimensions,
                                indices_2, dim));
  }
};


/*
 *  Extern operator overloadings
 */

template<typename T>
core::ndarray<T> operator+(const core::ndarray<T>& op1, const T& op2) {
  ndarray<T> copy(op1);
  return copy += op2;
}

template<typename T>
core::ndarray<T> operator+(const T& op1, const core::ndarray<T>& op2) {
  return op2 + op1;
}

template<typename T>
core::ndarray<T> operator+(const core::ndarray<T>& op1,
                           const core::ndarray<T>& op2) {
  return ndarray<T>::apply_operation(operation::add<T>, op1, op2);
}


template<typename T>
core::ndarray<T> operator*(const core::ndarray<T>& op1, const T& op2) {
  ndarray<T> copy(op1);
  return copy *= op2;
}

template<typename T>
core::ndarray<T> operator*(const T& op1, const core::ndarray<T>& op2) {
  return op2 * op1;
}

template<typename T>
core::ndarray<T> operator*(const core::ndarray<T>& op1,
                           const core::ndarray<T>& op2) {
  return ndarray<T>::apply_operation(operation::multiplie<T>, op1, op2);
}

template<typename T>
core::ndarray<T> operator-(const core::ndarray<T>& op1, const T& op2) {
  ndarray<T> copy(op1);
  return copy -= op2;
}

template<typename T>
core::ndarray<T> operator-(const T& op1, const core::ndarray<T>& op2) {
  return - op2 + op1;
}

template<typename T>
core::ndarray<T> operator-(const core::ndarray<T>& op1,
                           const core::ndarray<T>& op2) {
  return core::ndarray<T>::apply_operation(operation::substract<T>, op1, op2);
}

template<typename T>
core::ndarray<T> invert(const core::ndarray<T>& op) {
  core::ndarray<T> copy(op);
  return copy.invert();
}

template<typename T>
core::ndarray<T>& transpose(const core::ndarray<T>& op,
                            std::size_t dim_1 = 0, std::size_t dim_2 = 1) {
  core::ndarray<T> copy(op);
  return copy.transpose(dim_1, dim_2);
}

template<typename T>
core::ndarray<T> operator/(const core::ndarray<T>& op1, const T& op2) {
  ndarray<T> copy(op1);
  return copy /= op2;
}

template<typename T>
core::ndarray<T> operator/(const T& op1, const core::ndarray<T>& op2) {
  return invert(op2) * op1;
}

template<typename T>
core::ndarray<T> operator/(const core::ndarray<T>& op1,
                           const core::ndarray<T>& op2) {
  return ndarray<T>::apply_operation(operation::divide<T>, op1, op2);
}

}  // namespace core

#endif /* NDARRAY_H */
