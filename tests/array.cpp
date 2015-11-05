// Copyright 2015 Thomas Foetisch <thomas.foetisch@gmail.com>

#include <iostream>

#include <core.hpp>
#include <core/io.hpp>


int main(int argc, char *argv[]) {

  double id(0.0);
  core::ndarray<double> array({4, 3, 4});
  for (std::size_t k(0); k < array.get_dimension(2); ++k)
    for (std::size_t j(0); j < array.get_dimension(1); ++j)
      for (std::size_t i(0); i < array.get_dimension(0); ++i) {
        array.index(k, j, i) = id;
        id += 1.0;
      }
  std::cout << array << std::endl;

  const core::ndarray<double>& carray(array);
  std::cout << "The slice array(1:2, 1:2, 0:3):" << std::endl << std::endl;
  std::cout << carray.index(core::slice(1, 2),
                            core::slice(1, 2),
                            core::slice(0, 3)) << std::endl;

  core::ndarray<double> my_slice(carray.index(core::slice(1, 2),
                                              core::slice(1, 2),
                                              core::slice(0, 3)));

  std::cout << "Test de broadcasting: " << std::endl;
  core::ndarray<double> singleton({1});
  singleton(0) = 3.;
  
  std::cout << my_slice + singleton << std::endl;

  my_slice += singleton;

  std::cout << my_slice << std::endl;
  
  return 0;
}
