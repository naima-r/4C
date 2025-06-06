// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_POROFLUID_PRESSURE_BASED_ELAST_DYN_HPP
#define FOUR_C_POROFLUID_PRESSURE_BASED_ELAST_DYN_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*! entry point for the solution of porofluid-elasticity problems */
void porofluid_elast_dyn(int restart /* do we have to perform a restart?  */
);

FOUR_C_NAMESPACE_CLOSE

#endif
