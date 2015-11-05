#ifndef SLICE_H
#define SLICE_H

#include <cstddef>
#include <vector>


namespace core {

/*
 *  Uniform slice representation
 *    slice(0, 5) == {0, 1, 2, 3, 4, 5}
 *    slice(0, 5, 2) == {0, 2, 4}
 */
class slice {
 public:
  slice(): start(1), end(0), stride(1) {}
  slice(std::size_t _start, std::size_t _end, std::size_t _stride = 1)
      : start(_start), end(_end), stride(_stride) {}

  std::size_t get_size() const {
    return start <= end ? (end - start) / stride + 1 : 0;
  }

  std::size_t operator()(std::size_t i) const {
    return start + i * stride;
  }

  std::size_t get_start() const { return start; }
  std::size_t get_stride() const { return stride; }
  
 private:
  std::size_t start, end, stride;
};


/*
 *  Integer index to slice conversion helper
 */
template<typename index_type>
void slice_from_index(std::vector<slice>::iterator slices,
                      index_type);

template<>
void slice_from_index<std::size_t>(std::vector<slice>::iterator s,
                                   std::size_t i) {
  *s = slice(i, i, 1);
}

template<>
void slice_from_index<core::slice>(std::vector<slice>::iterator s,
                                   core::slice i) {
  *s = i;
}

}  // namespace core

#endif /* SLICE_H */
