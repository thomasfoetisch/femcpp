#ifndef LAGRANGE_P1BULLE_H
#define LAGRANGE_P1BULLE_H

#include <core.hpp>
#include <fe.hpp>

#include "lagrange.hpp"

namespace fe {

template<unsigned int dim>
class lagrange_p1bubble;

template<>
class lagrange_p1bubble<2>: public fe::lagrange_2d<1> {
 public:
};

/*template<>
class lagrange_p1bubble<3>: public fe::lagrange_3d<1> {
 public:
 };*/

}

#endif /* LAGRANGE_P1BULLE_H */
