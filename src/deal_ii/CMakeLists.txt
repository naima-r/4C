# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

four_c_auto_define_module(NO_CYCLES)

# deal.II support is a pure extension of the core module and may not depend on any other module
four_c_add_dependency(${AUTO_DEFINED_MODULE_NAME} core)
