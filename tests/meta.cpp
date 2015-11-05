// Copyright 2015 Thomas Foetisch <thomas.foetisch@gmail.com>

#include <cstddef>
#include <iostream>

#include <core.hpp>
#include <meta.hpp>

typedef meta::typelist<char, std::string, int> my_list;

template<typename ... args>
struct slicing_return_type {
  static constexpr bool is_index =
      meta::is_homogeneous<meta::typelist<args...> >::value
      and
      meta::equal<typename meta::get<meta::typelist<args...>, 0>::type,
                  std::size_t>::value;
  static constexpr bool is_slice = not is_index;
  
  typedef
  typename meta::if_<is_slice,
                     core::ndarray<double>,
                     double>::type
  type;
};

template<bool is_slice, typename ... args>
struct index_impl;

template<typename ... args>
struct index_impl<true, args...> {
  static
  typename slicing_return_type<args...>::type
  do_it(args ... arg_values) {
    std::cout << "sliced" << std::endl;
    return core::ndarray<double>();
  }
};

template<typename ... args>
struct index_impl<false, args...> {
  template<typename arg1, typename...args_>
  struct accumulate_offsets {
    static
    void do_it(std::size_t& offset,
               core::multiindex::iterator& dim,
               arg1 i1, args_ ... arg_values) {
      offset += (*dim) * i1;
      accumulate_offsets<args_...>::do_it(offset, ++dim, arg_values...);
    }
  };

  template<typename arg>
  struct accumulate_offsets<arg> {
    static
    void do_it(std::size_t& offset,
               core::multiindex::iterator& dim,
               arg i1) {
      offset += i1;
    }
  };
  
  static
  typename slicing_return_type<args...>::type
  do_it(const core::ndarray<double>& elements,
        const core::multiindex& offsets,
        args ... arg_values) {
    std::size_t offset(0);

    accumulate_offsets<args...>::do_it(offset, offsets.begin(), arg_values...);
    
    return elements(offset);
  }
};

template<typename ... args>
auto index(args ... arg_values) -> typename slicing_return_type<args...>::type {
  return index_impl<slicing_return_type<args...>::is_slice, args...>::do_it(arg_values...);
}

int main(int argc, char *argv[]) {

  std::cout << meta::length<my_list>::value << std::endl;
  meta::get<my_list, 1>::type var("hello les copains.");
  std::cout << var << std::endl;

  std::cout << meta::is_homogeneous<meta::typelist<double> >::value << std::endl;
  std::cout << meta::is_homogeneous<meta::typelist<double, double> >::value << std::endl;
  std::cout << meta::is_homogeneous<meta::typelist<double, int> >::value << std::endl;
  std::cout << meta::is_homogeneous<meta::typelist<double, double, double> >::value << std::endl;

  std::cout << index(1lu, 0lu, 0lu) << std::endl;
  
  return 0;
}
