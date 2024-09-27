// Copyright (c) 2020-2024 Simons Foundation
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0.txt
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Authors: Thomas Hahn, Olivier Parcollet, Nils Wentzell

/**
 * @file
 * @brief Provides a cross product for 3-dimensional vectors or other arrays/views of rank 1.
 */

#pragma once

#include "../declarations.hpp"
#include "../macros.hpp"
#include "../traits.hpp"

namespace nda::linalg {

  /**
   * @ingroup linalg_tools
   * @brief Compute the cross product of two 3-dimensional vectors.
   *
   * @tparam V Vector type.
   * @param x Left hand side vector.
   * @param y Right hand side vector.
   * @return nda::array of rank 1 containing the cross product of the two vectors.
   */
  template <typename V>
  auto cross_product(V const &x, V const &y) {
    EXPECTS_WITH_MESSAGE(x.shape()[0] == 3, "nda::linalg::cross_product: Only defined for 3-dimensional vectors");
    EXPECTS_WITH_MESSAGE(y.shape()[0] == 3, "nda::linalg::cross_product: Only defined for 3-dimensional vectors");
    array<get_value_t<V>, 1> r(3);
    r(0) = x(1) * y(2) - y(1) * x(2);
    r(1) = -x(0) * y(2) + y(0) * x(2);
    r(2) = x(0) * y(1) - y(0) * x(1);
    return r;
  }

} // namespace nda::linalg
