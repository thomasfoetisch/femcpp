// Copyright 2015 Thomas Foetisch <thomas.foetisch@gmail.com>

#ifndef FE_H
#define FE_H

#include <exception>

#include "meta.hpp"

namespace fe {

class error: public std::exception {};

template<typename finite_element_list>
class mixed_fe {
 public:
};

}  // namespace fe

#include "fe/lagrange.hpp"
#include "fe/lagrange_p1bubble.hpp"

#endif /* FE_H */

