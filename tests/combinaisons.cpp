
#include <geometry/mesh.hpp>

/*template<typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<std::vector<T> >& v) {
  for (auto c: v) {
    for (auto item: c)
      stream << item << " ";
    stream << std::endl;
  }
  return stream;
  }*/

int main(int argc, char *argv[]) {
  core::array<std::size_t> set(5);
  for (auto i: {0, 1, 2, 3, 4})
    set(i) = i;

  for (auto i: {1, 2, 3, 4, 5})
    std::cout << combination(i, set) << std::endl;

  return 0;
}
