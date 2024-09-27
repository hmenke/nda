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
// Authors: Thomas Hahn, Sergei Iskakov, Olivier Parcollet, Nils Wentzell

/**
 * @file
 * @brief Provides additional operators for std::complex and other arithmetic types.
 */

#ifndef STDUTILS_COMPLEX_H
#define STDUTILS_COMPLEX_H

#include <complex>
#include <concepts>
#include <type_traits>

using namespace std::literals::complex_literals;

namespace std { // has to be in the right namespace for ADL

  /**
   * @addtogroup utils_std
   * @{
   */

// define operators (+,-,*,/) for std::complex and various other arithmetic types
#define IMPL_OP(OP)                                                                                                                                  \
  /** @brief Implementation of operator `OP` for std::complex and some other arithmetic type. */                                                     \
  template <typename T, typename U>                                                                                                                  \
    requires(std::is_arithmetic_v<T> and std::is_arithmetic_v<U> and std::common_with<T, U>)                                                         \
  auto operator OP(std::complex<T> const &x, U y) {                                                                                                  \
    using C = std::complex<std::common_type_t<T, U>>;                                                                                                \
    return C(x.real(), x.imag()) OP C(y);                                                                                                            \
  }                                                                                                                                                  \
                                                                                                                                                     \
  /** @brief Implementation of operator `OP` for some other arithmetic type and std::complex. */                                                     \
  template <typename T, typename U>                                                                                                                  \
    requires(std::is_arithmetic_v<T> and std::is_arithmetic_v<U> and std::common_with<T, U>)                                                         \
  auto operator OP(T x, std::complex<U> const &y) {                                                                                                  \
    using C = std::complex<std::common_type_t<T, U>>;                                                                                                \
    return C(x) OP C(y.real(), y.imag());                                                                                                            \
  }                                                                                                                                                  \
                                                                                                                                                     \
  /** @brief Implementation of operator `OP` for two std::complex types with different value types. */                                               \
  template <typename T, typename U>                                                                                                                  \
    requires(std::is_arithmetic_v<T> and std::is_arithmetic_v<U> and std::common_with<T, U> and !std::is_same_v<T, U>)                               \
  auto operator OP(std::complex<T> const &x, std::complex<U> const &y) {                                                                             \
    using C = std::complex<std::common_type_t<T, U>>;                                                                                                \
    return C(x.real(), x.imag()) OP C(y.real(), y.imag());                                                                                           \
  }

  IMPL_OP(+)
  IMPL_OP(-)
  IMPL_OP(*)
  IMPL_OP(/)

#undef IMPL_OP

  /** @}  */

} // namespace std

#endif // STDUTILS_COMPLEX_H
