// Copyright 2015 Thomas Foetisch <thomas.foetisch@gmail.com>

#include <fe.hpp>

int main(int argc, char *argv[]) {
  {
    fe::lagrange_2d<0> l_p02d;

    assert(l_p02d.spatial_dimensions() == 2);
    assert(l_p02d.total_dof_number() == 1);

    assert(l_p02d.dof_number(0) == 0);
    assert(l_p02d.dof_number(1) == 0);
    assert(l_p02d.dof_number(2) == 1);

    assert(l_p02d.phi(0, core::make_2d_point(0, 0)) == 1.0);

    assert(l_p02d.node_type(0) == 2);
  }

  {
    fe::lagrange_2d<1> l_p12d;
    assert(l_p12d.total_dof_number() == 3);

    assert(l_p12d.dof_number(0) == 3);
    assert(l_p12d.dof_number(1) == 0);
    assert(l_p12d.dof_number(2) == 0);

    assert(abs(l_p12d.phi(0, core::make_2d_point(0.5, 0.5)) - 0.5) < 1.e-16);
    assert(abs(l_p12d.phi(2, core::make_2d_point(0.0, 1.0)) - 1.0) < 1.e-16);

    assert(l_p12d.node_type(0) == 0);
    assert(l_p12d.node_type(1) == 0);
    assert(l_p12d.node_type(2) == 0);
  }

  {
    fe::lagrange_2d<4> l_p42d;
    assert(l_p42d.total_dof_number() == 15);

    assert(l_p42d.dof_number(0) == 3);
    assert(l_p42d.dof_number(1) == 9);
    assert(l_p42d.dof_number(2) == 3);

    assert(l_p42d.node_type(0) == 0);
    assert(l_p42d.node_type(1) == 1);
    assert(l_p42d.node_type(4) == 0);
    assert(l_p42d.node_type(5) == 1);
    assert(l_p42d.node_type(6) == 2);
    assert(l_p42d.node_type(14) == 0);
  }

  {
    fe::lagrange_2d<1> p1;
    fe::lagrange_2d<3> p3;

    fe::mixed_finite_element p1_p3;
    p1_p3.add_component(&p1);
    p1_p3.add_component(&p3);

    assert(p1_p3.manifold_id(0) == 0);
    assert(p1_p3.manifold_id(4) == 0);
    assert(p1_p3.manifold_id(6) == 1);
    assert(p1_p3.manifold_id(8) == 0);
    assert(p1_p3.manifold_id(9) == 2);
    assert(p1_p3.manifold_id(12) == 2);

    assert(p1_p3.spatial_dimensions() == 2);

    assert(p1_p3.total_dof_number() == 13);

    assert(p1_p3.dof_number(0) == 6);
    assert(p1_p3.dof_number(1) == 6);
    assert(p1_p3.dof_number(2) == 1);

    assert(p1_p3.dof_number_per_manifold(0) == 2);
    assert(p1_p3.dof_number_per_manifold(1) == 2);
    assert(p1_p3.dof_number_per_manifold(2) == 1);
    
    assert(p1_p3.node_type(0) == 0);
    assert(p1_p3.node_type(4) == 1);
    assert(p1_p3.node_type(6) == 0);
    assert(p1_p3.node_type(8) == 2);
    assert(p1_p3.node_type(9) == 1);
    assert(p1_p3.node_type(12) == 0);
  }
  
  return 0;
}
