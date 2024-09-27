// Copyright (c) 2019-2024 Simons Foundation
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
 * @brief Provides automatic assignment for std::vector.
 */

#pragma once

#include "../clef.hpp"

#include <utility>
#include <vector>

namespace nda::clef {

  /**
   * @addtogroup clef_autoassign
   * @{
   */

  namespace detail {

    // Helper function to auto assign to a std::vector object.
    template <typename T, typename RHS>
    void clef_auto_assign_std_vector_impl(T &x, RHS &&rhs) {
      x = std::forward<RHS>(rhs);
    }

    // Helper function to auto assign to a std::vector object.
    template <typename Expr, int... Is, typename T>
    void clef_auto_assign_std_vector_impl(T &x, make_fun_impl<Expr, Is...> &&rhs) { // NOLINT (why rvalue reference?)
      clef_auto_assign_subscript(x, std::forward<make_fun_impl<Expr, Is...>>(rhs));
    }

  } // namespace detail

  /**
   * @brief Overload of `clef_auto_assign_subscript` function for std::vector.
   *
   * @tparam T Value type of the std::vector.
   * @tparam F Callable type.
   * @param v std::vector object.
   * @param f Callable object.
   */
  template <typename T, typename F>
  void clef_auto_assign_subscript(std::vector<T> &v, F f) {
    for (size_t i = 0; i < v.size(); ++i) detail::clef_auto_assign_std_vector_impl(v[i], f(i));
  }

  /** @} */

} // namespace nda::clef
