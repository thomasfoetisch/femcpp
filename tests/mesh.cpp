// Copyright 2015 Thomas Foetisch <thomas.foetisch@gmail.com>

#include <utility>
#include <memory>
#include <map>

#include <geometry/mesh.hpp>


class current_finite_element {
 public:
  current_finite_element() {}
  void set(const geometry::mesh& m,
           const fe::finite_element_base* fe,
           std::size_t element, const core::point& x) {
    core::ndarray<double>
        phi({fe->total_dof_number()}),
        dphi({fe->spatial_dimensions(),
                fe->total_dof_number()});
    
    for (std::size_t j(0); j < phi.get_dimension(0); ++j) {
      phi(j) = fe->phi(j, x);
      for (std::size_t d(0); d < fe->spatial_dimensions(); ++d) {
        dphi.index(d, j) = 0.0;
        for (std::size_t p(0); p < fe->spatial_dimensions(); ++p)
          dphi.index(d, j) += m.jmt(element, d, p)
              * fe->dphi(p, j, x);
      }
    }
  }

  const core::ndarray<double>& phi() const;
  const core::ndarray<double>& dphi() const;

  std::size_t get_dof_number() const { return ev_phi.get_dimension(0); }

 private:
  core::ndarray<double> ev_phi, ev_dphi;
  
  
};

template<std::size_t space_dim>
class navierstokes_smagorinsky {
 public:
  navierstokes_smagorinsky(geometry::finite_element_space* sol_fes)
      : rho(1.0),
        laminar_viscosity(0.1),
        smagorinsky_constant(0.1), smagorinsky_caracteristic_length(0.1),
        epsilon_regularization(1e-13), kronecker({space_dim, space_dim}),
        sol(sol_fes) {
    for (std::size_t i(0); i < space_dim; ++i)
      for (std::size_t j(0); j < space_dim; ++j)
        kronecker.index(i, j) = i == j ? 1.0 : 0.0;
  }
  
  core::ndarray<double>
  evaluate_velocity_velocity(const core::point& x, std::size_t element,
                             std::size_t test_comp,
                             std::size_t trial_comp,
                             const current_finite_element& test_fe,
                             const current_finite_element& trial_fe) {
      
    core::ndarray<double> result({test_fe.get_dof_number(),
                                  trial_fe.get_dof_number()});
                                   

    core::ndarray<double> epsilon({space_dim, space_dim});
    double epsilon_norm(0.0);
    for (std::size_t i(0); i < space_dim; ++i)
      for (std::size_t j(0); j < space_dim; ++j) {
        epsilon.index(i, j) = 0.5 * (sol(element, j, i, x)
                                     + sol(element, i, j, x));
        epsilon_norm += std::pow(epsilon.index(i, j), 2);
      }
    epsilon_norm = std::sqrt(epsilon_norm);

    const double mu(laminar_viscosity
                    + (2 * smagorinsky_constant
                       * smagorinsky_caracteristic_length
                       * epsilon_norm));
    
    for (std::size_t i(0); i < test_fe.get_dof_number(); ++i)
      for (std::size_t j(0); j < test_fe.get_dof_number(); ++j) {
        result.index(i, j) += rho
            * sol(element, trial_comp, test_comp, x)
            * test_fe.phi()(j)
            * trial_fe.phi()(i);
    
        for (std::size_t p(0); p < space_dim; ++p)
          result.index(i, j) += rho
              * kronecker.index(trial_comp, test_comp)
              * trial_fe.dphi().index(p, j)
              * test_fe.dphi().index(p, i);

        if (epsilon_norm > epsilon_regularization)
          for (std::size_t p(0); p < space_dim; ++p)
            for (std::size_t q(0); q < space_dim; ++q)
              result.index(i, j) += 2.0 * rho
                  * smagorinsky_constant * smagorinsky_caracteristic_length
                  * epsilon.index(trial_comp, p) * trial_fe.dphi().index(p, j)
                  * epsilon.index(test_comp, q) * test_fe.dphi().index(q, i)
                  / epsilon_norm;

        result.index(i, j) += mu
            * trial_fe.dphi().index(trial_comp, i) * test_fe.dphi().index(test_comp, j);
        
        for (std::size_t p(0); p < space_dim; ++p)
          result.index(i, j) += mu * kronecker.index(trial_comp, test_comp)
              * trial_fe.dphi().index(p, j) * test_fe.dphi().index(p, i);
      }

    return result;
  }

  
  core::ndarray<double>
  evaluate_pressure_velocity(const core::point& x, std::size_t element_id,
                             std::size_t test_comp, std::size_t trial_comp,
                             const current_finite_element& test_fe,
                             const current_finite_element& trial_fe) {
    core::ndarray<double> result({test_fe.get_dof_number(),
                                  trial_fe.get_dof_number()});
    for (std::size_t i(0); i < test_fe.get_dof_number(); ++i)
      for (std::size_t j(0); j < trial_fe.get_dof_number(); ++j)
        result.index(i, j) += - trial_fe.phi().index(j) * test_fe.dphi().index(test_comp, i);

    return result;
  }

  
  core::ndarray<double>
  evaluate_velocity_pressure(const core::point& x, std::size_t element_id,
                             std::size_t trial_comp, std::size_t test_comp,
                             const current_finite_element& test_fe,
                             const current_finite_element& trial_fe) {
    core::ndarray<double> result({test_fe.get_dof_number(),
                                  trial_fe.get_dof_number()});

    for (std::size_t i(0); i < test_fe.get_dof_number(); ++i)
      for (std::size_t j(0); j < trial_fe.get_dof_number(); ++j)
        result.index(i, j) += - trial_fe.phi().index(i) * test_fe.dphi().index(trial_comp, j);

    return result;
  }

  
  core::ndarray<double>
  evaluate_velocity(const core::point& x, std::size_t element_id,
                    std::size_t test_comp,
                    const current_finite_element& test_fe) {
    // todo
  }


  core::ndarray<double>
  evaluate_pressure(const core::point& x, std::size_t element_id,
                    std::size_t test_comp,
                    const current_finite_element& test_fe) {
    // todo
  }

 private:
  const double rho, laminar_viscosity,
    smagorinsky_constant, smagorinsky_caracteristic_length;
  const double epsilon_regularization;
  core::ndarray<double> kronecker;
  
  geometry::field sol;
};


class quadrature_formula {
 public:
  const core::point& node(std::size_t id) const {return nodes[id]; }
  double weight(std::size_t id) const { return weights[id]; }

  std::size_t get_nodes_number() const { return nodes.size(); }
  std::size_t get_polynomial_order() const;

 private:
  std::vector<core::point> nodes;
  std::vector<double> weights;
};


class sparse_matrix {
 public:
  sparse_matrix(std::size_t _n_lines, std::size_t _n_columns)
      : n_lines(_n_lines), n_columns(_n_columns) {

  }
  const double& operator()(std::size_t i, std::size_t j) const {
    auto component(components.find(std::make_pair(i, j)));
    if (component == components.end())
      return zero;
    else
      return component->second;
  }
  double& operator()(std::size_t i, std::size_t j) {
    auto component(components.find(std::make_pair(i, j)));
    if (component == components.end())
      component = components.insert(std::pair<coordinate, double>(std::make_pair(i, j), 0.0)).first;
    return component->second;
  }

 private:
  std::size_t n_lines, n_columns;
  typedef std::pair<std::size_t, std::size_t> coordinate;
  std::map<coordinate, double> components;

  static constexpr double zero = 0.0;
};


class linear_system {
 public:
  linear_system(const geometry::mesh& d,
                geometry::finite_element_space *_test_space,
                geometry::finite_element_space *_trial_space)
      : domain(d),
        test_space(_test_space),
        trial_space(_trial_space),
        m(test_space->get_dofs_number(),
          trial_space->get_dofs_number()) {}

  typedef core::callable<core::ndarray<double>,
                         const core::point&, std::size_t,
                         std::size_t, std::size_t,
                         const current_finite_element&,
                         const current_finite_element&> system_assembler_t;
  typedef core::callable<core::ndarray<double>,
                         const core::point&, std::size_t,
                         std::size_t,
                         const current_finite_element&> rhs_assembler_t;
  
  system_assembler_t get_system_block_assembler(std::size_t i, std::size_t j) {
    return system_block_assemblers.find(std::make_pair(i, j))->second;
  }
  rhs_assembler_t get_rhs_block_assembler(std::size_t j) {
    return rhs_block_assemblers.find(j)->second;
  }
  const quadrature_formula& get_block_quadrature(std::size_t i, std::size_t j) {
    return quadratures.find(std::make_pair(i, j))->second;
  }

  void set_quadrature_formula(const quadrature_formula& q) {
    for (std::size_t j(0); j < trial_space->get_components_number(); ++j)
      for (std::size_t i(0); i < test_space->get_components_number(); ++i)
        quadratures[std::make_pair(i, j)] = q;
  }
  void set_quadrature_formula(std::size_t i, std::size_t j, const quadrature_formula& q) {
    quadratures[std::make_pair(i, j)] = q;
  }

  void set_system_assembler(std::size_t i, std::size_t j, system_assembler_t a) {
    system_block_assemblers[std::make_pair(i, j)] = a;
  }
  void set_system_assembler(system_assembler_t a) {
    for (std::size_t j(0); j < trial_space->get_components_number(); ++j)
      for (std::size_t i(0); i < test_space->get_components_number(); ++i)
        system_block_assemblers[std::make_pair(i, j)] = a;
  }

  void set_rhs_assembler(std::size_t i, rhs_assembler_t a) {
    rhs_block_assemblers[i] = a;
  }
  void set_rhs_assembler(rhs_assembler_t a) {
    for (std::size_t j(0); j < trial_space->get_components_number(); ++j)
      rhs_block_assemblers[j] = a;
  }
  
  void assemble_linear_system() {
    for (std::size_t comp_i(0); comp_i < test_space->get_components_number(); ++comp_i) {
      for (std::size_t comp_j(0); comp_j < trial_space->get_components_number(); ++comp_j) {
        current_finite_element
            current_test_fe,
            current_trial_fe;

        
        system_assembler_t assembler(get_system_block_assembler(comp_i, comp_j));
        const quadrature_formula& q(get_block_quadrature(comp_i, comp_j));

        for (std::size_t element(0); element < domain.get_elements_number(); ++element) {

          core::ndarray<double> ke({test_space->get_finite_element(comp_i)->total_dof_number(),
                                    trial_space->get_finite_element(comp_j)->total_dof_number()});

          for (std::size_t q_i(0); q_i < q.get_nodes_number(); ++q_i) {

            current_test_fe.set(domain, test_space->get_finite_element(comp_i), element, q.node(q_i));
            current_trial_fe.set(domain, trial_space->get_finite_element(comp_j), element, q.node(q_i));

            ke += domain.jac(element)
                * q.weight(q_i) * assembler(test_space->map_to_element(element, q.node(q_i)),
                                            element,
                                            comp_i, comp_j,
                                            current_test_fe, current_trial_fe);
          }

          for (std::size_t i(0); i < test_space->get_finite_element(comp_i)->total_dof_number(); ++i)
            for (std::size_t j(0); j < trial_space->get_finite_element(comp_j)->total_dof_number(); ++j)
              m(test_space->get_global_dof_id(comp_i, element, i),
                trial_space->get_global_dof_id(comp_j, element, j)) += ke.index(i, j);
        }
      }
    }
  }

  void assemble_right_hand_side() {
    for (std::size_t comp_j(0); comp_j < test_space->get_components_number(); ++comp_j) {
      current_finite_element current_test_fe;

        
      rhs_assembler_t assembler(get_rhs_block_assembler(comp_j));
      const quadrature_formula& q(get_block_quadrature(comp_j, 0));

      for (std::size_t element(0); element < domain.get_elements_number(); ++element) {

        core::ndarray<double> ke({test_space->get_finite_element(comp_j)->total_dof_number()});

        for (std::size_t q_i(0); q_i < q.get_nodes_number(); ++q_i) {

          current_test_fe.set(domain, test_space->get_finite_element(comp_j), element, q.node(q_i));
        
          ke += domain.jac(element)
              * q.weight(q_i) * assembler(test_space->map_to_element(element, q.node(q_i)),
                                          element,
                                          comp_j,
                                          current_test_fe);
        }

        for (std::size_t i(0); i < test_space->get_finite_element(comp_j)->total_dof_number(); ++i)
          rhs(test_space->get_global_dof_id(comp_j, element, i)) += ke(i);
      }
    }
  }
  

 private:
  geometry::mesh domain;
  geometry::finite_element_space *test_space, *trial_space;

  typedef std::pair<std::size_t, std::size_t> block_coordinate;
  std::map<block_coordinate, system_assembler_t> system_block_assemblers;
  std::map<std::size_t, rhs_assembler_t> rhs_block_assemblers;
  std::map<block_coordinate, quadrature_formula> quadratures;
  
  sparse_matrix m;
  core::array<double> rhs;
};


  

int main(int argc, char *argv[]) {
  geometry::mesh square(geometry::generate::square_mesh(1.0, 1.0, 1, 1));

  fe::lagrange_2d<3> p3;
  fe::lagrange_2d<1> p1;

  fe::mixed_finite_element p3_p1;
  p3_p1.add_component(&p3);
  p3_p1.add_component(&p1);
  
  geometry::finite_element_space v_h(square, &p3_p1);
  geometry::field u(&v_h);

  navierstokes_smagorinsky<2> ns(&v_h);

  linear_system finite_element_system(square, &v_h, &v_h);
  
  auto uu_block(core::wrap_callable(ns, &navierstokes_smagorinsky<2>::evaluate_velocity_velocity));
  auto up_block(core::wrap_callable(ns, &navierstokes_smagorinsky<2>::evaluate_velocity_pressure));
  auto pu_block(core::wrap_callable(ns, &navierstokes_smagorinsky<2>::evaluate_pressure_velocity));

  auto u_block(core::wrap_callable(ns, &navierstokes_smagorinsky<2>::evaluate_velocity));
  auto p_block(core::wrap_callable(ns, &navierstokes_smagorinsky<2>::evaluate_pressure));

  for (std::size_t i(0); i < 2; ++i) {
    for (std::size_t j(0); j < 2; ++j)
      finite_element_system.set_system_assembler(i, j, uu_block);
    finite_element_system.set_system_assembler(2, i, pu_block);
    finite_element_system.set_system_assembler(i, 2, up_block);
  }

  finite_element_system.set_quadrature_formula(quadrature_formula());

  finite_element_system.assemble_linear_system();
  finite_element_system.assemble_right_hand_side();
    
  return 0;
}
