// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_UTILS_HPP
#define FOUR_C_SCATRA_UTILS_HPP

#include "4C_config.hpp"

#include "4C_fem_condition.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_fixedsizematrix.hpp"


FOUR_C_NAMESPACE_OPEN

namespace ScaTra::ScaTraUtils
{
  //! check, if s2i condition definition is consistent
  void check_consistency_of_s2_i_conditions(
      std::shared_ptr<Core::FE::Discretization> discretization);

  //! check, if nodes of input conditions equal s2i kinetics condition
  void check_consistency_with_s2_i_kinetics_condition(const std::string& condition_to_be_tested,
      std::shared_ptr<Core::FE::Discretization> discretization);

  //! determine the s2i kinetics butler volmer conditions modeling
  //! simplified growth
  std::vector<const Core::Conditions::Condition*>
  get_s2i_kinetics_butler_volmer_simplified_growth_conditions(
      std::shared_ptr<Core::FE::Discretization> discretization);

  //! Calculate the reconstructed nodal gradient at a node by means of mean value averaging
  template <const int dim>
  std::shared_ptr<Core::LinAlg::MultiVector<double>> compute_gradient_at_nodes_mean_average(
      Core::FE::Discretization& discret, const Core::LinAlg::Vector<double>& state,
      const int scatra_dofid);

  //! Calculate the reconstructed nodal gradient at a node by means of mean value averaging
  template <const int dim, Core::FE::CellType distype>
  Core::LinAlg::Matrix<dim, 1> do_mean_value_averaging_of_element_gradient_node(
      Core::FE::Discretization& discret, std::vector<const Core::Elements::Element*> elements,
      Core::LinAlg::Vector<double>& phinp_node, const int nodegid, const int scatra_dofid);

}  // namespace ScaTra::ScaTraUtils
FOUR_C_NAMESPACE_CLOSE

#endif
