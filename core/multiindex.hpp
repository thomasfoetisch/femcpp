#ifndef MULTIINDEX_H
#define MULTIINDEX_H

#include <cstddef>
#include <vector>
#include <algorithm>


namespace core {

typedef std::vector<std::size_t> multiindex;


multiindex operator*(const multiindex& i1, const multiindex& i2) {
  multiindex result(i1);
  auto it1(result.begin());
  auto it2(i2.begin());
  while (it1 != result.end()) {
    *it1 *= *it2;
    ++it1, ++it2;
  }

  return result;
}

multiindex operator+(const multiindex& i1, const multiindex& i2) {
  multiindex result(i1);
  auto it1(result.begin());
  auto it2(i2.begin());
  while (it1 != result.end()) {
    *it1 += *it2;
    ++it1, ++it2;
  }

  return result;
}

bool increment(multiindex& i, const multiindex& d) {
  std::size_t digit(d.size());
  while(digit and (i[digit - 1] + 1 >= d[digit - 1]))
    --digit;

  if(digit) {
    ++i[digit - 1];
    for (std::size_t k(digit); k < d.size(); ++k)
      i[k] = 0;
    return true;
  } else {
    return false;
  }
}

bool is_broadcastable(const multiindex& i1, const multiindex& i2) {
  const std::size_t common_size(std::min(i1.size(), i2.size()));

  for (std::size_t i(0); i < common_size; ++i)
    if (i1[i1.size() - 1 - i] != i2[i2.size() - 1 - i]
        and
        i1[i1.size() - 1 - i] != 1 and i2[i2.size() - 1 - i] != 1)
      return false;

  return true;
}

bool is_compound_broadcastable(const multiindex& i1, const multiindex& i2) {
  const std::size_t common_size(std::min(i1.size(), i2.size()));

  for (std::size_t i(0); i < common_size; ++i)
    if (i1[i1.size() - 1 - i] != i2[i2.size() - 1 - i]
        and
        i2[i2.size() - 1 - i] != 1)
      return false;

  return true;
}

multiindex normalize_index(const multiindex& i, std::size_t length) {
  multiindex result(length, 1);
  std::copy(i.begin(), i.end(), result.begin() + length - i.size());
  return result;
}

multiindex broadcasted_dimensions(const multiindex& d_1,
                                  const multiindex& d_2) {
  const std::size_t max_size(std::max(d_1.size(), d_2.size()));
  multiindex result(max_size, 1);

  for (std::size_t i(0); i < d_1.size(); ++i)
    result[max_size - i - 1] =
        std::max(result[max_size - i - 1],
                 d_1[d_1.size() - i - 1]);

  for (std::size_t i(0); i < d_2.size(); ++i)
    result[max_size - i - 1] =
        std::max(result[max_size - i - 1],
                 d_2[d_2.size() - i - 1]);

  return result;
}

bool broadcast_increment(multiindex& i_1, const multiindex& d_1,
                         multiindex& i_2, const multiindex& d_2) {
  /*
   *  Find which digit can be incremented: (assume d_1.size() == d_2.size())
   */
  std::size_t digit(d_1.size());
  while (digit
         and
         (i_1[digit - 1] + 1 >= d_1[digit - 1]
          and
          i_2[digit - 1] + 1 >= d_2[digit - 1])) {
    --digit;
  }


  /*
   *  Increment the digit as long it's not broadcasted
   */
  if (digit) {
    if (d_1[digit - 1] > 1)
      ++i_1[digit - 1];
    if (d_2[digit - 1] > 1)
      ++i_2[digit - 1];


    /*
     *  Zero out the successive digits:
     */
    for (std::size_t k(digit); k < d_2.size(); ++k) {
      i_1[k] = 0;
      i_2[k] = 0;
    }

    
    /*
     *  Incrementation was successful:
     */
    return true;
  } else {
    /*
     *  All the digits reached the maximum value, hence
     *  the incrementation fails:
     */
    return false;
  }
}

}  // namespace core

#endif /* MULTIINDEX_H */
