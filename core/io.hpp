#ifndef IO_H
#define IO_H

#include <vector>
#include <ostream>

#include "multiindex.hpp"
#include "ndarray.hpp"

namespace std {

std::ostream& operator<<(std::ostream& stream, const core::multiindex& id) {
  for (auto i: id)
    stream << i << ", ";
  return stream;
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, const ::core::ndarray<T>& a) {
  if (a.get_rank() > 2) {
    core::multiindex
        indices(a.get_rank() - 2, 0),
        dimensions(a.get_dimensions().begin(),
                   a.get_dimensions().begin() + a.get_rank() - 2);
    std::vector<core::slice> slices(a.get_rank());
    slices[a.get_rank() - 2] =
        core::slice(0, a.get_dimension(a.get_rank() - 2) - 1);
    slices[a.get_rank() - 1] =
        core::slice(0, a.get_dimension(a.get_rank() - 1) - 1);
    do {
      for (std::size_t i(0); i < slices.size() - 2; ++i)
        slices[i] = core::slice(indices[i], indices[i]);
      stream << " showing slice (" << indices << ":, :): " << std::endl;
      stream << core::ndarray<double>(a, slices).squeeze() << std::endl;
    } while (core::increment(indices, dimensions));
  } else {
    for (std::size_t i(0); i < a.get_dimension(0); ++i) {
      stream << "  ";
      for (std::size_t j(0); j < a.get_dimension(1); ++j)
        stream << a.index(i, j) << " ";
      stream << std::endl;
    }
  }
  return stream;
}

}  // namespace std

#endif /* IO_H */
