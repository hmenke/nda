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
 * @brief Provides some custom implementations of standard mathematical functions used for lazy, coefficient-wise array
 * operations.
 */

#pragma once

#include "./concepts.hpp"
#include "./map.hpp"
#include "./traits.hpp"

#include <cmath>
#include <complex>
#include <type_traits>
#include <utility>

namespace nda {

  /**
   * @addtogroup av_math
   * @{
   */

  /**
   * @brief Get the real part of a scalar.
   *
   * @tparam T Scalar type.
   * @param t Scalar value.
   * @return Real part of the scalar.
   */
  template <typename T>
  auto real(T t)
    requires(nda::is_scalar_v<T>)
  {
    if constexpr (is_complex_v<T>) {
      return std::real(t);
    } else {
      return t;
    }
  }

  /**
   * @brief Get the complex conjugate of a scalar.
   *
   * @tparam T Scalar type.
   * @param t Scalar value.
   * @return The given scalar if it is not complex, otherwise its complex conjugate.
   */
  template <typename T>
  auto conj(T t)
    requires(nda::is_scalar_v<T>)
  {
    if constexpr (is_complex_v<T>) {
      return std::conj(t);
    } else {
      return t;
    }
  }

  /**
   * @brief Get the squared absolute value of a double.
   *
   * @param x Double value.
   * @return Squared absolute value of the given double.
   */
  inline double abs2(double x) { return x * x; }

  /**
   * @brief Get the squared absolute value of a std::complex<double>.
   *
   * @param z std::complex<double> value.
   * @return Squared absolute value of the given complex number.
   */
  inline double abs2(std::complex<double> z) { return (conj(z) * z).real(); }

  /**
   * @brief Check if a std::complex<double> is NaN.
   *
   * @param z std::complex<double> value.
   * @return True if either the real or imaginary part of the given complex number is `NaN`, false otherwise.
   */
  inline bool isnan(std::complex<double> const &z) { return std::isnan(z.real()) or std::isnan(z.imag()); }

  /**
   * @brief Calculate the integer power of an integer.
   *
   * @tparam T Integer type.
   * @param x Base value.
   * @param n Exponent value.
   * @return The result of the base raised to the power of the exponent.
   */
  template <typename T>
  T pow(T x, int n)
    requires(std::is_integral_v<T>)
  {
    T r = 1;
    for (int i = 0; i < n; ++i) r *= x;
    return r;
  }

  /**
   * @brief Lazy, coefficient-wise power function for nda::Array types.
   *
   * @tparam A nda::Array type.
   * @param a nda::Array object.
   * @param p Exponent value.
   * @return A lazy nda::expr_call object.
   */
  template <Array A>
  auto pow(A &&a, double p) {
    return nda::map([p](auto const &x) {
      using std::pow;
      return pow(x, p);
    })(std::forward<A>(a));
  }

  /// Wrapper for nda::conj.
  struct conj_f {
    /// Function call operator that forwards the call to nda::conj.
    auto operator()(auto const &x) const { return conj(x); };
  };

  /**
   * @brief Lazy, coefficient-wise complex conjugate function for nda::Array types.
   *
   * @tparam A nda::Array type.
   * @param a nda::Array object.
   * @return A lazy nda::expr_call object if the array is complex valued, otherwise the array itself.
   */
  template <Array A>
  decltype(auto) conj(A &&a) {
    if constexpr (is_complex_v<get_value_t<A>>)
      return nda::map(conj_f{})(std::forward<A>(a));
    else
      return std::forward<A>(a);
  }

  /** @} */

} // namespace nda
