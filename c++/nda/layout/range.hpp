// Copyright (c) 2018 Commissariat à l'énergie atomique et aux énergies alternatives (CEA)
// Copyright (c) 2018 Centre national de la recherche scientifique (CNRS)
// Copyright (c) 2018-2024 Simons Foundation
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
 * @brief Includes the itertools header and provides some additional utilities.
 */

#pragma once

#include "../traits.hpp"

#include <itertools/itertools.hpp>

#include <ostream>
#include <type_traits>

namespace nda {

  /**
   * @addtogroup layout_utils
   * @{
   */

  /// Using declaration for itertools::range.
  using itertools::range;

  /**
   * @brief Mimics Python's `...` syntax.
   *
   * @details While itertools's `range::all_t` mimics Python's `:`, `ellipsis` mimics Python's `...`. It is repeated as
   * much as necessary to match the number of dimensions of an array/view when used to access elements/slices.
   */
  struct ellipsis : range::all_t {};

  /**
   * @brief Write `nda::range::all_t` to a std::ostream as `_`.
   *
   * @param os Output stream.
   * @return Reference to the output stream.
   */
  inline std::ostream &operator<<(std::ostream &os, range::all_t) noexcept { return os << "_"; }

  /**
   * @brief Write nda::ellipsis to a std::ostream as `___`.
   *
   * @param os Output stream.
   * @return Reference to the output stream.
   */
  inline std::ostream &operator<<(std::ostream &os, ellipsis) noexcept { return os << "___"; }

  /// Constexpr variable that is true if the parameter pack `Args` contains an nda::ellipsis.
  template <typename... Args>
  constexpr bool ellipsis_is_present = is_any_of<ellipsis, std::remove_cvref_t<Args>...>;

  /**
   * @brief Constexpr variable that is true if the type `T` is either an `nda::range`, an `nda::range::all_t` or an
   * nda::ellipsis.
   */
  template <typename T>
  constexpr bool is_range_or_ellipsis = is_any_of<std::remove_cvref_t<T>, range, range::all_t, ellipsis>;

  /** @} */

} // namespace nda
