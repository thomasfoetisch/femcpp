// Copyright 2015 Thomas Foetisch <thomas.foetisch@gmail.com>

#ifndef META_H
#define META_H

namespace meta {

template<typename T>
class not_equal {
 public:
  not_equal(const T& v): value(v) {}

  bool operator()(const T& t) const { return value != t; };
 private:
  const T value;
};


struct null_type {};

template<typename head, typename ... tail>
struct typelist {
  typedef head head_type;
  typedef typelist<tail...> tail_type;
};

template<typename head>
struct typelist<head> {
  typedef head head_type;
  typedef null_type tail_type;
};

template<typename list>
struct length {
  static constexpr std::size_t value = 1 + length<typename list::tail_type>::value;
};

template<>
struct length<null_type>  {
  static constexpr std::size_t value = 0;
};

template<typename list, std::size_t item>
struct get {
  typedef typename get<typename list::tail_type, item - 1>::type type;
};

template<typename list>
struct get<list, 0> {
  typedef typename list::head_type type;
};

template<typename t1, typename t2>
struct equal { static constexpr bool value = false; };

template<typename t>
struct equal<t, t> { static constexpr bool value = true; };

template<bool cond, typename t1, typename t2> struct if_;

template<typename t1, typename t2> struct if_<true, t1, t2> {
  typedef t1 type;
};

template<typename t1, typename t2> struct if_<false, t1, t2> {
  typedef t2 type;
};

template<typename list>
struct is_homogeneous {
  static constexpr bool value
  = equal<typename get<list, 0>::type,
          typename get<list, 1>::type>::value
      and
      is_homogeneous<typename list::tail_type>::value;
};

template<typename t>
struct is_homogeneous< typelist<t> > {
  static constexpr bool value = true;
};

template<typename list>
struct is_integral_typelist {
  static constexpr bool value =
      std::is_integral<typename list::head_type>::value
      and
      is_integral_typelist<typename list::tail_type>::value;
};

template<typename t>
struct is_integral_typelist<typelist<t> > {
  static constexpr bool value = std::is_integral<t>::value;
};


} // namespace meta


#endif /* META_H */
