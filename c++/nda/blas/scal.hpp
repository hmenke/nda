// Copyright (c) 2023-2024 Simons Foundation
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
// Authors: Thomas Hahn, Miguel Morales, Nils Wentzell

/**
 * @file
 * @brief Provides a generic interface to the BLAS `scal` routine.
 */

#pragma once

#include "./interface/cxx_interface.hpp"
#include "./tools.hpp"
#include "../concepts.hpp"
#include "../device.hpp"
#include "../mem/address_space.hpp"
#include "../traits.hpp"

namespace nda::blas {

  /**
   * @ingroup linalg_blas
   * @brief Interface to the BLAS `scal` routine.
   *
   * @details Scales a vector by a constant. This function calculates
   * \f[
   *   \mathbf{x} \leftarrow \alpha \mathbf{x} ;,
   * \f]
   * where \f$ \alpha \f$ is a scalar constant and \f$ \mathbf{x} \f$ is a vector.
   *
   * @tparam X nda::MemoryVector or a conjugate array expression.
   * @param alpha Input scalar.
   * @param x Input/Output vector to be scaled.
   */
  template <typename X>
    requires(MemoryVector<X> or is_conj_array_expr<X>)
  void scal(get_value_t<X> alpha, X &&x) { // NOLINT (temporary views are allowed here)
    static_assert(is_blas_lapack_v<get_value_t<X>>, "Error in nda::blas::scal: Value type of vector is incompatible with blas");

    if constexpr (mem::on_host<X>) {
      f77::scal(x.size(), alpha, x.data(), x.indexmap().strides()[0]);
    } else {
#if defined(NDA_HAVE_DEVICE)
      device::scal(x.size(), alpha, x.data(), x.indexmap().strides()[0]);
#else
      compile_error_no_gpu();
#endif
    }
  }

} // namespace nda::blas
