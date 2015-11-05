// Copyright 2015 Thomas Foetisch <thomas.foetisch@gmail.com>

#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <numeric>

#include "core.hpp"
#include "fe.hpp"

template<typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& v) {
  for (auto item: v)
    stream << item << " ";
  return stream;
}

template<typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<std::vector<T> >& v) {
  for (auto c: v) {
    stream << c;
    stream << std::endl;
  }
  return stream;
}

namespace geometry {

class mesh {
 public:
  mesh(const core::ndarray<double>& _nodes,
       const core::ndarray<std::size_t>& _elements)
      : nodes(_nodes), elements(_elements) {
    for (std::size_t elem(0); elem < get_elements_number(); ++elem)
      std::sort(&elements.index(elem, 0ul),
                &elements.index(elem, 0ul) + get_dimension() + 1);

    build_metrics();
  }

  unsigned int get_dimension() const { return elements.get_dimension(1) - 1; }
  unsigned int get_embedding_dimension() const { return nodes.get_dimension(1); }
  std::size_t get_nodes_number() const { return nodes.get_dimension(0); }
  std::size_t get_elements_number() const { return elements.get_dimension(0); }

  mesh get_boundary() const;
  mesh get_boundary(core::ndarray<std::size_t>& boundary_to_mesh_injection) const;

  double jac(std::size_t element) const { return jacobian(element); };
  double jmt(std::size_t element, std::size_t i, std::size_t j) const {
    return jmt_matrix({element, i, j});
  }

  const core::ndarray<std::size_t>& get_elements() const { return elements; }
  const core::ndarray<double>& get_nodes() const { return nodes; }
  
 private:
  core::ndarray<double> nodes;
  core::ndarray<std::size_t> elements;

  core::ndarray<double> jacobian, jmt_matrix;
  
  void build_metrics() {
    jacobian = core::ndarray<double>({get_elements_number()});
    jmt_matrix = core::ndarray<double>({get_elements_number(),
                                        get_dimension(), get_dimension()});

    core::slice nodes_slice(0, get_dimension() - 1);
    for (std::size_t element(0); element < get_elements_number(); ++element) {
      core::ndarray<double> change_of_var({get_dimension(), get_dimension()});
      for (std::size_t i(1); i <= get_dimension(); ++i)
        change_of_var.index(nodes_slice, i) =
            nodes.index(elements.index(element, i), nodes_slice)
            - nodes.index(elements.index(element, 0), nodes_slice);

      jmt_matrix.index(element, nodes_slice, nodes_slice) =
          core::invert_matrix(change_of_var).transpose(0, 1);
      
      jacobian(element) = 1.0 / core::determinant(change_of_var);
    }
  }
};


class finite_element_space {
 public:
  explicit finite_element_space(const mesh& m,
                                const fe::finite_element_base* _finite_element)
      : finite_element(_finite_element), dof_map(), total_dof_number(0) {
    build_dof_map(m);
  }

  std::size_t get_dofs_number() const { return total_dof_number; }

  std::size_t get_components_number() const {
    return finite_element->get_component_number();
  }

  std::size_t get_component_local_dof_number(std::size_t component_id) const {
    return finite_element->get_component(component_id)->total_dof_number();
  }

  core::point map_to_reference(std::size_t element_id,
                                       const core::point& x) const;
  core::point map_to_element(std::size_t element_id,
                             const core::point& x) const;
  
  std::size_t get_global_dof_id(std::size_t element_id,
                                std::size_t component_id,
                                std::size_t local_dof_id) const {
    return dof_map.index(element_id, finite_element->get_dof_id(component_id, local_dof_id));
  }
  
  const fe::finite_element_base* get_finite_element(std::size_t component_id) const {
    return finite_element->get_component(component_id);
  }


 private:
  const fe::finite_element_base* finite_element;
  core::ndarray<std::size_t> dof_map;
  std::size_t total_dof_number;

  typedef std::vector<std::size_t> boundary_manifold;
  std::vector<boundary_manifold> build_mesh_manifold(const mesh& m,
                                                     std::size_t manifold_dimension) {
    const core::ndarray<std::size_t>& elements(m.get_elements());
    std::set<boundary_manifold> manifolds;
    for (std::size_t elem(0); elem < m.get_elements_number(); ++elem) {
      core::slice nodes(0, elements.get_dimension(1) - 1);
      std::vector<boundary_manifold>
          element_manifolds(core::combination(manifold_dimension + 1,
                                              elements.index(elem, nodes)));

      for (auto manifold: element_manifolds)
        manifolds.insert(manifold);
    }

    std::vector<boundary_manifold> result(manifolds.size());
    std::copy(manifolds.begin(), manifolds.end(),
              result.begin());

    return result;
  }

  static
  std::size_t compute_global_manifold_id(const std::vector<boundary_manifold>& global_list,
                                         const boundary_manifold& local_manifold) {
    auto position(std::lower_bound(global_list.begin(), global_list.end(), local_manifold));
    return std::distance(global_list.begin(), position);
  }
  
  void build_dof_map(const mesh& m) {
    const core::ndarray<std::size_t>& elements(m.get_elements());

    // first build the intersection manifolds:
    typedef std::vector<boundary_manifold> manifold_collection;
    std::vector<manifold_collection> intersection_manifolds;
    for (std::size_t manifold_dimension(0);
         manifold_dimension <= finite_element->spatial_dimensions();
         ++manifold_dimension)
      intersection_manifolds.push_back(build_mesh_manifold(m, manifold_dimension));


    // build offsets.
    std::vector<std::size_t> dof_offsets(finite_element->spatial_dimensions() + 1);
    dof_offsets[0] = 0;
    for (std::size_t space_dimension(1);
         space_dimension < dof_offsets.size();
         ++space_dimension) {
      dof_offsets[space_dimension] =
          dof_offsets[space_dimension - 1]
          + intersection_manifolds[space_dimension - 1].size()
          * finite_element->dof_number_per_manifold(space_dimension - 1);
    }

    total_dof_number = dof_offsets.back()
        + intersection_manifolds.back().size()
        * finite_element->dof_number_per_manifold(finite_element->spatial_dimensions());


    // allocate dofs:
    dof_map = core::ndarray<std::size_t>({m.get_elements_number(),
                                          finite_element->total_dof_number()});

    // then, for each element we rebuild each manifold, and look for the index:
    for (std::size_t elem(0); elem < m.get_elements_number(); ++elem) {
      std::vector<manifold_collection> element_manifolds;
      for (std::size_t manifold_dimension(0);
           manifold_dimension <= finite_element->spatial_dimensions();
           ++manifold_dimension)
        element_manifolds.push_back(combination(manifold_dimension + 1,
                                                elements.index(elem, core::slice(0, elements.get_dimension(1) - 1))));

      for (std::size_t local_dof(0);
           local_dof < finite_element->total_dof_number();
           ++local_dof) {
        unsigned int manifold_dimension(finite_element->node_type(local_dof));
        
        dof_map.index(elem, local_dof) =
            dof_offsets[manifold_dimension]
            + finite_element->dof_number_per_manifold(manifold_dimension)
            * compute_global_manifold_id(intersection_manifolds[manifold_dimension],
                                         element_manifolds[manifold_dimension][finite_element->manifold_id(local_dof)])
            + finite_element->manifold_local_dof_id(local_dof);
      }
    }
  }
};


class field {
 public:
  explicit field(const finite_element_space* _fes)
      : fes(_fes), coefficients(fes->get_dofs_number()) {}

  //core::array<double> evaluate(const core::point& x) const;
  //double evaluate(std::size_t component_id, const core::point& x) const;
  
  //core::array<double> evaluate_on_element(std::size_t element, const core::point& x) const;
  double evaluate_on_element(std::size_t component_id,
                             std::size_t element_id,
                             const core::point& x) const {
    double result(0.0);
    const core::point& x_hat(fes->map_to_reference(element_id, x));
    for (std::size_t i(0); i < fes->get_component_local_dof_number(component_id); ++i)
      result += coefficients(fes->get_global_dof_id(element_id, component_id, i))
          * fes->get_finite_element(component_id)->phi(i, x_hat);
    return result;
  }

  double operator()(std::size_t element, std::size_t component, std::size_t derivative, const core::point& x) const;

  const double& dof(std::size_t component_id,
                        std::size_t element_id,
                        std::size_t component_dof_id) const {
    return coefficients(fes->get_global_dof_id(component_id,
                                               element_id,
                                               component_dof_id));
  }

  double& dof(std::size_t component_id,
              std::size_t element_id,
              std::size_t component_dof_id) {
    return coefficients(fes->get_global_dof_id(component_id,
                                               element_id,
                                               component_dof_id));
  }
  
 private:
  const finite_element_space* fes;
  core::array<double> coefficients;
};

namespace generate {

mesh square_mesh(double l_x, double l_y, std::size_t n_el_x, std::size_t n_el_y) {
  const double h_x(l_x / n_el_x), h_y(l_y / n_el_y);
  
  core::ndarray<double> nodes({n_el_x + 1, n_el_y + 1, 2});
  core::ndarray<std::size_t> elements({n_el_x, n_el_y, 2, 3});

  core::ndarray<std::size_t> node_ids({n_el_x + 1, n_el_y + 1});
  std::size_t current_node_id(0);
  
    for (std::size_t j(0); j < n_el_y + 1; ++j) 
      for (std::size_t i(0); i < n_el_x + 1; ++i) {

        nodes.index(i, j, 0) = - l_x / 2.0 + i * h_x;
        nodes.index(i, j, 1) = - l_y / 2.0 + j * h_y;

        node_ids.index(i, j) = current_node_id;
        ++current_node_id;
    }

  for (std::size_t i(0); i < n_el_x; ++i)
    for (std::size_t j(0); j < n_el_y; ++j) {
      elements.index(i, j, 0, 0) = node_ids.index(i, j);
      elements.index(i, j, 0, 1) = node_ids.index(i + 1, j);
      elements.index(i, j, 0, 2) = node_ids.index(i + 1, j + 1);

      elements.index(i, j, 1, 0) = node_ids.index(i, j);
      elements.index(i, j, 1, 1) = node_ids.index(i + 1, j + 1);
      elements.index(i, j, 1, 2) = node_ids.index(i, j + 1);
    }

  nodes.reshape({(n_el_x + 1) * (n_el_y + 1), 2});
  elements.reshape({n_el_x * n_el_y * 2, 3});

  return mesh(nodes, elements);
}

}


}  // namespace geometry

#endif /* MESH_H */
