#ifndef COMBINATIONS_H
#define COMBINATIONS_H

namespace core {

template<typename T>
std::vector<typename T::value_type>
extract_items(const T& array,
              const std::vector<std::size_t>& selection) {
  std::vector<typename T::value_type> result(selection.size());

  auto index(selection.begin());
  auto item(result.begin());

  while (index != selection.end()) {
    *item = array(*index);
    ++index;
    ++item;
  }

  return result;
}

template<typename T>
std::vector<std::vector<T> >
combination(std::size_t k, const core::ndarray<T>& array) {
  const std::size_t n(array.get_length());

  if (k > n)
    throw std::exception();

  
  std::vector<std::vector<T> > result;
  result.reserve(nchoosek(n, k));

  std::vector<std::size_t> selection(k);
  std::iota(selection.begin(), selection.end(), 0);

  bool done(false);
  while (not done) {
    result.push_back(extract_items(array, selection));

    int digit(k - 1);
    while (digit >= 0) {
      if (selection[digit] < n - (k - digit)) {
        ++selection[digit];
        std::iota(selection.begin() + digit + 1,
                  selection.end(),
                  selection[digit] + 1);
        digit = -1;
      } else {
        --digit;
        if (digit < 0)
          done = true;
      }
    }
  }

  return result;
}

}

#endif /* COMBINATIONS_H */
