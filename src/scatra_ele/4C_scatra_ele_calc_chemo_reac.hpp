// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_ELE_CALC_CHEMO_REAC_HPP
#define FOUR_C_SCATRA_ELE_CALC_CHEMO_REAC_HPP

#include "4C_config.hpp"

#include "4C_mat_scatra_chemotaxis.hpp"
#include "4C_mat_scatra_reaction.hpp"
#include "4C_scatra_ele_calc.hpp"
#include "4C_scatra_ele_calc_advanced_reaction.hpp"
#include "4C_scatra_ele_calc_chemo.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret
{
  namespace Elements
  {
    template <Core::FE::CellType distype, int probdim = Core::FE::dim<distype>>
    class ScaTraEleCalcChemoReac : public ScaTraEleCalcChemo<distype, probdim>,
                                   public ScaTraEleCalcAdvReac<distype, probdim>
    {
     private:
      //! private constructor, since we are a Singleton.
      ScaTraEleCalcChemoReac(
          const int numdofpernode, const int numscal, const std::string& disname);

      using my = ScaTraEleCalc<distype, probdim>;
      using chemo = ScaTraEleCalcChemo<distype, probdim>;
      using advreac = ScaTraEleCalcAdvReac<distype, probdim>;

     public:
      //! Singleton access method
      static ScaTraEleCalcChemoReac<distype, probdim>* instance(
          const int numdofpernode, const int numscal, const std::string& disname);

     protected:
      //! get the material parameters
      void get_material_params(
          const Core::Elements::Element* ele,  //!< the element we are dealing with
          std::vector<double>& densn,          //!< density at t_(n)
          std::vector<double>& densnp,         //!< density at t_(n+1) or t_(n+alpha_F)
          std::vector<double>& densam,         //!< density at t_(n+alpha_M)
          double& visc,                        //!< fluid viscosity
          const int iquad = -1                 //!< id of current gauss point (default = -1)
          ) override;

    };  // end class ScaTraEleCalcChemoReac

  }  // namespace Elements

}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
