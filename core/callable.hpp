#ifndef CALLABLE_H
#define CALLABLE_H

namespace core {

template<typename ret, typename ...args>
class callable_impl {
 public:
  virtual ~callable_impl() {}
  virtual ret operator()(args...) = 0;
};

template<typename ret, typename ...args>
class callable_free_function: public callable_impl<ret, args...> {
 public:
  callable_free_function(ret (*_f_ptr)(args...))
  : f_ptr(_f_ptr) {}
  
  ret operator()(args... arg_values) {
    return f_ptr(arg_values...);
  }

 private:
  typedef ret (*func_type)(args...);
  func_type f_ptr;
};

template<typename class_type, typename ret, typename ...args>
class callable_member_function: public callable_impl<ret, args...> {
 public:
  callable_member_function(class_type& inst, ret (class_type::*_mem_f_ptr)(args...))
      : instance(inst), mem_f_ptr(_mem_f_ptr) {}
  
  ret operator()(args ...arg_values) {
    return (instance.*mem_f_ptr)(arg_values...);
  }
  
 private:
  class_type& instance;
  typedef ret (class_type::*mem_func_type)(args...);
  mem_func_type mem_f_ptr;
};

template<typename ret, typename ...args>
class callable {
 public:
  callable(): callable_ptr(nullptr) {}
  
  callable(ret (*free_f)(args...))
  : callable_ptr(new callable_free_function<ret, args...>(free_f)) {}

  template<typename class_type>
  callable(class_type& instance, ret (class_type::*mem_f)(args...))
      : callable_ptr(new callable_member_function<class_type, ret, args...>(instance, mem_f)) {}
  
  ret operator()(args... arg_values) {
    return (*callable_ptr)(arg_values...);
  }
      
 private:
  std::shared_ptr<callable_impl<ret, args...> > callable_ptr;
};

template<typename ret, typename...args>
callable<ret, args...> wrap_callable(ret (*f_ptr)(args...)) {
  return callable<ret, args...>(f_ptr);
}

template<typename class_type, typename ret, typename...args>
callable<ret, args...> wrap_callable(class_type& inst, ret (class_type::*mem_ptr)(args...)) {
  return callable<ret, args...>(inst, mem_ptr);
}


}

#endif /* CALLABLE_H */
