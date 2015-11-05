#ifndef LAGRANGE_H
#define LAGRANGE_H

#include <set>
#include <vector>
#include <iostream>  // tmp

#include <core.hpp>
#include <fe.hpp>

namespace fe {

using core::point;

class finite_element_base {
 public:
  virtual ~finite_element_base() {}

  virtual std::size_t manifold_id(std::size_t local_dof) const = 0;
  virtual std::size_t manifold_local_dof_id(std::size_t local_dof) const = 0;
  virtual unsigned int spatial_dimensions() const = 0;
  virtual unsigned int total_dof_number() const = 0;
  virtual unsigned int dof_number(unsigned int intersection_manifold_dimension) const = 0;
  virtual unsigned int dof_number_per_manifold(unsigned int intersection_manifold_dimension) const = 0;

  virtual double phi(std::size_t dof_id, const point& x) const = 0;
  virtual double phi(std::size_t component_id, std::size_t dof_id, const point& x) const = 0;

  virtual double dphi(std::size_t dof_id, const point& x) const = 0;
  virtual double dphi(std::size_t component_id, std::size_t dof_id, const point& x) const = 0;
  
  virtual unsigned int node_type(std::size_t local_dof) const = 0;
  virtual std::size_t get_component_number() const = 0;
  virtual const finite_element_base* get_component(std::size_t component_id) const = 0;
  virtual const std::size_t get_dof_id(std::size_t component_id, std::size_t component_dof_id) const = 0;

  virtual void add_yourself(std::vector<const finite_element_base*>& comp) const {
    comp.push_back(this);
  }
};

template<unsigned int degree>
class lagrange_2d: public finite_element_base {
 public:
  lagrange_2d(): local_dof_map({total_dof_number(), 3}) {
    // build the local_dof_id to (manifold_id, local_manifold_node_id) map
    typedef std::set<std::size_t> node_collection;

    std::vector<node_collection> vertex_nodes(3);
    std::vector<node_collection> edge_nodes(3);
    std::vector<node_collection> inner_nodes(1);

    for (std::size_t local_dof(0);
         local_dof < total_dof_number();
         ++local_dof) {
      unsigned int i(0), j(0);
      dof_id_to_grid_coordinates(i, j, local_dof);
      const unsigned int type(node_type(i, j));
      switch (type) {
        case 0: {
          if (i == j)
            vertex_nodes[0].insert(local_dof);
          else if (j == 0)
            vertex_nodes[1].insert(local_dof);
          else
            vertex_nodes[2].insert(local_dof);
          break;
        }

        case 1: {
          if (j == 0)
            edge_nodes[0].insert(local_dof);
          else if (i == 0)
            edge_nodes[1].insert(local_dof);
          else
            edge_nodes[2].insert(local_dof);
          break;
        }

        case 2: {
          inner_nodes[0].insert(local_dof);
          break;
        }

        default: throw fe::error();
      }
    }

    {
      unsigned int manifold_id(0), manifold_local_dof(0);
      for (auto vertex_manifold: vertex_nodes) {
        manifold_local_dof = 0;
        for (auto local_dof: vertex_manifold) {
          local_dof_map.index(local_dof, 0) = node_type(local_dof);
          local_dof_map.index(local_dof, 1) = manifold_id;
          local_dof_map.index(local_dof, 2) = manifold_local_dof;

          ++manifold_local_dof;
        }
        ++manifold_id;
      }
    }

    {
      unsigned int manifold_id(0), manifold_local_dof(0);
      for (auto edge_manifold: edge_nodes) {
        manifold_local_dof = 0;
        for (auto local_dof: edge_manifold) {
          local_dof_map.index(local_dof, 0) = node_type(local_dof);
          local_dof_map.index(local_dof, 1) = manifold_id;
          local_dof_map.index(local_dof, 2) = manifold_local_dof;

          ++manifold_local_dof;
        }
        ++manifold_id;
      }
    }

    {
      unsigned int manifold_id(0), manifold_local_dof(0);
      for (auto inner_manifold: inner_nodes) {
        manifold_local_dof = 0;
        for (auto local_dof: inner_manifold) {
          local_dof_map.index(local_dof, 0) = node_type(local_dof);
          local_dof_map.index(local_dof, 1) = manifold_id;
          local_dof_map.index(local_dof, 2) = manifold_local_dof;

          ++manifold_local_dof;
        }
        ++manifold_id;
      }
    }
  }

  std::size_t manifold_id(std::size_t local_dof) const {
    return local_dof_map.index(local_dof, 1);
  }
  
  std::size_t manifold_local_dof_id(std::size_t local_dof) const {
    return local_dof_map.index(local_dof, 2);
  }

  unsigned int spatial_dimensions() const { return 2; }
  unsigned int total_dof_number() const { return (degree + 1) * (degree + 2) / 2; }
  unsigned int dof_number(unsigned int intersection_manifold_dimension) const {
    switch (intersection_manifold_dimension) {
      case 0: return degree < 1 ? 0 : 3;
      case 1: return degree < 2 ? 0 : 3 * (degree - 1);
      case 2: return degree < 3 ? (degree == 0 ? 1 : 0) : (degree - 1) * (degree - 2) / 2;
      default: throw fe::error();
    }
  }

  unsigned int dof_number_per_manifold(unsigned int intersection_manifold_dimension) const {
    switch (intersection_manifold_dimension) {
      case 0: return dof_number(intersection_manifold_dimension) / 3;
      case 1: return dof_number(intersection_manifold_dimension) / 3;
      case 2: return dof_number(intersection_manifold_dimension);
      default: throw fe::error();
    }
  }

  double phi(std::size_t dof_id, const point& x) const {
    const double node_interval(degree == 0 ? 0.0 : 1.0 / degree);
    double result(1.0);
    unsigned int i(0), j(0);

    dof_id_to_grid_coordinates(i, j, dof_id);

    for (unsigned int isoline(0); isoline < i; ++isoline)
      result *= lambda1(x(0), x(1)) - isoline * node_interval;

    for (unsigned int isoline(0); isoline < j; ++isoline)
      result *= lambda2(x(0), x(1)) - isoline * node_interval;

    for (unsigned int isoline(0); isoline < degree - (i + j); ++isoline)
      result *= lambda0(x(0), x(1)) - isoline * node_interval;

    return result;
  }

  double phi(std::size_t component_id, std::size_t dof_id, const point& x) const {
    return phi(dof_id, x);
  }

  double dphi(std::size_t dof_id, const point& x) const {}
  double dphi(std::size_t component_id, std::size_t dof_id, const point& x) const {}

  unsigned int node_type(std::size_t dof_id) const {
    unsigned int i(0), j(0);
    dof_id_to_grid_coordinates(i, j, dof_id);

    return node_type(i, j);
  }

  std::size_t get_component_number() const { return 1; }
  const finite_element_base* get_component(std::size_t component_id) const { return this; }

  const std::size_t get_dof_id(std::size_t component_id, std::size_t component_dof_id) const {
    return component_dof_id;
  }
  
 private:
  core::ndarray<std::size_t> local_dof_map;

  static inline double lambda0(double x, double y) { return 1.0 - x - y; }
  static inline double lambda1(double x, double y) { return x; }
  static inline double lambda2(double x, double y) { return y; }

  void dof_id_to_grid_coordinates(unsigned int& i,
                                  unsigned int& j,
                                  unsigned int dof_id) const {
    if (dof_id >= total_dof_number())
      throw fe::error();

    i = 0;
    j = 0;

    unsigned int node_per_layer(degree + 1);
    while (dof_id) {
      if (dof_id + 1 > node_per_layer) {
        ++j;
        dof_id -= node_per_layer;
      } else {
        i = dof_id;
        dof_id = 0;
      }
      --node_per_layer;
    }
  }

  unsigned int node_type(unsigned int i, unsigned int j) const {
    if (degree == 0)
      return 2;

    if (i == 0 or j == 0 or i + j == degree) {
      if ( i == 0 and j == degree)
        return 0;
      if (j == 0 and i == degree)
        return 0;
      if (i == 0 and j == 0)
        return 0;

      return 1;
    } else {
      return 2;
    }
  }
};


class mixed_finite_element: public finite_element_base {
 public:
  void add_component(const finite_element_base* finite_element) {
    finite_element->add_yourself(components);
  }

  std::size_t manifold_id(std::size_t local_dof) const {
    std::size_t component_id(local_dof_to_component(local_dof));
    return components[component_id]->manifold_id(local_dof);
  }
  
  std::size_t manifold_local_dof_id(std::size_t local_dof) const {
    std::size_t component_id(local_dof_to_component(local_dof));
    std::size_t manifold_dimension(components[component_id]->node_type(local_dof));
    std::size_t cumulative_manifold_dof_number(0);
    for (std::size_t i(0); i < component_id; ++i)
      cumulative_manifold_dof_number += components[i]->dof_number_per_manifold(manifold_dimension);
    
    return cumulative_manifold_dof_number
        + components[component_id]->manifold_local_dof_id(local_dof);
  }
      

  unsigned int spatial_dimensions() const {
    return components.front()->spatial_dimensions();
  }

  unsigned int total_dof_number() const {
    std::size_t cumulative_total_dof_number(0);
    for (auto component: components)
      cumulative_total_dof_number += component->total_dof_number();
    return cumulative_total_dof_number;
  }

  unsigned int dof_number(unsigned int intersection_manifold_dimension) const {
    std::size_t cumulative_dof_number(0);
    for (auto component: components)
      cumulative_dof_number += component->dof_number(intersection_manifold_dimension);
    return cumulative_dof_number;
  }

  unsigned int dof_number_per_manifold(unsigned int intersection_manifold_dimension) const {
    std::size_t cumulative_dof_number_per_manifold(0);
    for (auto component: components)
      cumulative_dof_number_per_manifold
          += component->dof_number_per_manifold(intersection_manifold_dimension);
    return cumulative_dof_number_per_manifold;
  }
  
  double phi(std::size_t dof_id, const point& x) const {
    std::size_t component_id(local_dof_to_component(dof_id));
    return phi(component_id, dof_id, x);
  }

  double phi(std::size_t component_id, std::size_t dof_id, const point& x) const {
    return components[component_id]->phi(dof_id, x);
  }

  double dphi(std::size_t dof_id, const point& x) const {}
  double dphi(std::size_t component_id, std::size_t dof_id, const point& x) const {}

  unsigned int node_type(std::size_t local_dof) const {
    std::size_t component_id(local_dof_to_component(local_dof));
    return components[component_id]->node_type(local_dof);
  }

  std::size_t get_component_number() const {
    std::size_t size(0);
    for (auto component: components)
      size += component->get_component_number();
    return size;
  }

  const finite_element_base* get_component(std::size_t component_id) const {
    return components[component_id];
  }

  const std::size_t get_dof_id(std::size_t component_id,
                               std::size_t component_dof_id) const {
    std::size_t id(0);
    for (std::size_t i(0); i < component_id - 1; ++i)
      id += components[i]->total_dof_number();
    return id + components[component_id]->get_dof_id(0, component_dof_id);
  }

 protected:
  void add_yourself(std::vector<const finite_element_base*>& comp) const {
    for (auto finite_element: components)
      comp.push_back(finite_element);
  }
  
 private:
  std::vector<const finite_element_base*> components;

  std::size_t local_dof_to_component(std::size_t& local_dof) const {
    std::size_t component_id(0);
    while (local_dof >= components[component_id]->total_dof_number()) {
      local_dof -= components[component_id]->total_dof_number();
      ++component_id;
    }
    return component_id;
  }
};


/*
template<unsigned int degree>
class lagrange_3d: public finite_element_base {
 public:
  unsigned int spatial_dimensions() { return 3; }
  unsigned int total_dof_number() { return degree * (degree + 1) * (degree + 2) / 6; }
  unsigned int dof_number(unsigned int intersection_manifold_dimension) {
    switch (intersection_manifold_dimension) {
      case 0: return degree < 1 ? 0 : 4;
      case 1: return degree < 2 ? 0 : 6 * (degree - 1);
      case 2: return degree < 3 ? 0 : 4 * (degree - 1) * (degree - 2) / 2;
      case 3: return degree < 4 ? (degree == 0 ? 1 : 0) : (degree - 1) * (degree - 2) * (degree - 3) / 6;
      default: throw fe::error();
    }
  }


  static double phi(unsigned int dof_id, const point& x) {
    const double node_interval(1.0 / degree);
    double result(1.0);
    unsigned int i(0), j(0), k(0);

    dof_id_to_grid_coordinates(i, j, k, dof_id);

    for (unsigned int isoline(0); isoline < degree - (i + j + k); ++isoline)
      result *= lambda0(x(0), x(1)) - isoline * node_interval;

    for (unsigned int isoline(0); isoline < i; ++isoline)
      result *= lambda1(x(0), x(1)) - isoline * node_interval;

    for (unsigned int isoline(0); isoline < j; ++isoline)
      result *= lambda2(x(0), x(1)) - isoline * node_interval;

    for (unsigned int isoline(0); isoline < k; ++isoline)
      result *= lambda3(x(0), x(1)) - isoline * node_interval;

    return result;
  }

 private:
  lagrange_3d();
  
  static inline double lambda0(double x, double y, double z) { return 1.0 - x - y - z; }
  static inline double lambda1(double x, double y, double z) { return x; }
  static inline double lambda2(double x, double y, double z) { return y; }
  static inline double lambda3(double x, double y, double z) { return z; }

  static void dof_id_to_grid_coordinates(unsigned int& i,
                                         unsigned int& j,
                                         unsigned int& k,
                                         unsigned int dof_id) {

    if (dof_id >= total_dof_number())
      throw fe::error();
    
    i = 0;
    j = 0;
    k = 0;

    unsigned int node_per_slice((degree + 2) * (degree + 1) / 2);
    unsigned int node_per_layer(degree + 1);
    while (dof_id) {
      if (dof_id + 1 > node_per_slice) {
        ++k;
        dof_id -= node_per_slice;
      } else {
        while (dof_id) {
          if (dof_id + 1 > node_per_layer) {
            ++j;
            dof_id -= node_per_layer;
          } else {
            i = dof_id;
            dof_id = 0;
          }
          --node_per_layer;
        }
      }
      node_per_slice -= node_per_layer;
      node_per_layer -= 1;
    }
  }
};
*/

}  // namespace fe

#endif /* LAGRANGE_H */
