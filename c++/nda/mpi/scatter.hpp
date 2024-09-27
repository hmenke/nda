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
 * @brief Provides an MPI scatter function for nda::Array types.
 */

#pragma once

#include "../basic_functions.hpp"
#include "../concepts.hpp"
#include "../exceptions.hpp"
#include "../traits.hpp"

#include <mpi/mpi.hpp>

#include <type_traits>
#include <utility>
#include <vector>

/**
 * @ingroup av_mpi
 * @brief Specialization of the `mpi::lazy` class for nda::Array types and the `mpi::tag::scatter` tag.
 *
 * @details An object of this class is returned when scattering nda::Array objects across multiple MPI processes.
 *
 * It models an nda::ArrayInitializer, that means it can be used to initialize and assign to nda::basic_array and
 * nda::basic_array_view objects. The input array will be a chunked along its first dimension using `mpi::chunk_length`.
 *
 * See nda::mpi_scatter for an example.
 *
 * @tparam A nda::Array type to be scattered.
 */
template <nda::Array A>
struct mpi::lazy<mpi::tag::scatter, A> {
  /// Value type of the array/view.
  using value_type = typename std::decay_t<A>::value_type;

  /// Const view type of the array/view stored in the lazy object.
  using const_view_type = decltype(std::declval<const A>()());

  /// View of the array/view to be scattered.
  const_view_type rhs;

  /// MPI communicator.
  mpi::communicator comm;

  /// MPI root process.
  const int root{0}; // NOLINT (const is fine here)

  /// Should all processes receive the result. (doesn't make sense for scatter)
  const bool all{false}; // NOLINT (const is fine here)

  /**
   * @brief Compute the shape of the target array.
   *
   * @details The target shape will be the same as the input shape, except that the first dimension of the input array
   * is chunked into equal (as much as possible) parts using `mpi::chunk_length` and assigned to each MPI process.
   *
   * @warning This makes an MPI call.
   *
   * @return Shape of the target array.
   */
  [[nodiscard]] auto shape() const {
    auto dims = rhs.shape();
    long dim0 = dims[0];
    mpi::broadcast(dim0, comm, root);
    dims[0] = mpi::chunk_length(dim0, comm.size(), comm.rank());
    return dims;
  }

  /**
   * @brief Execute the lazy MPI operation and write the result to a target array/view.
   *
   * @tparam T nda::Array type of the target array/view.
   * @param target Target array/view.
   */
  template <nda::Array T>
  void invoke(T &&target) const { // NOLINT (temporary views are allowed here)
    if (not target.is_contiguous() or not target.has_positive_strides())
      NDA_RUNTIME_ERROR << "Error in MPI scatter for nda::Array: Target array needs to be contiguous with positive strides";

    static_assert(std::decay_t<A>::layout_t::stride_order_encoded == std::decay_t<T>::layout_t::stride_order_encoded,
                  "Error in MPI scatter for nda::Array: Incompatible stride orders");

    // special case for non-mpi runs
    if (not mpi::has_env) {
      target = rhs;
      return;
    }

    // get target shape and resize or check the target array
    auto dims = shape();
    resize_or_check_if_view(target, dims);

    // compute send counts, receive counts and memory displacements
    auto dim0       = rhs.extent(0);
    auto stride0    = rhs.indexmap().strides()[0];
    auto sendcounts = std::vector<int>(comm.size());
    auto displs     = std::vector<int>(comm.size() + 1, 0);
    int recvcount   = mpi::chunk_length(dim0, comm.size(), comm.rank()) * stride0;
    for (int r = 0; r < comm.size(); ++r) {
      sendcounts[r] = mpi::chunk_length(dim0, comm.size(), r) * stride0;
      displs[r + 1] = sendcounts[r] + displs[r];
    }

    // scatter the data
    auto mpi_value_type = mpi::mpi_type<value_type>::get();
    MPI_Scatterv((void *)rhs.data(), &sendcounts[0], &displs[0], mpi_value_type, (void *)target.data(), recvcount, mpi_value_type, root, comm.get());
  }
};

namespace nda {

  /**
   * @ingroup av_mpi
   * @brief Implementation of an MPI scatter for nda::basic_array or nda::basic_array_view types.
   *
   * @details Since the returned `mpi::lazy` object models an nda::ArrayInitializer, it can be used to initialize/assign
   * to nda::basic_array and nda::basic_array_view objects:
   *
   * @code{.cpp}
   * // create an array on all processes
   * nda::array<int, 2> arr(10, 4);
   *
   * // ...
   * // fill array on root process
   * // ...
   *
   * // scatter the array to all processes
   * nda::array<int, 2> res = mpi::scatter(arr);
   * @endcode
   *
   * Here, the array `res` will have a shape of `(10 / comm.size(), 4)`.
   *
   * @tparam A nda::basic_array or nda::basic_array_view type.
   * @param a Array or view to be scattered.
   * @param comm `mpi::communicator` object.
   * @param root Rank of the root process.
   * @param all Should all processes receive the result of the scatter (not used).
   * @return An `mpi::lazy` object modelling an nda::ArrayInitializer.
   */
  template <typename A>
  ArrayInitializer<std::remove_reference_t<A>> auto mpi_scatter(A &&a, mpi::communicator comm = {}, int root = 0, bool all = false)
    requires(is_regular_or_view_v<A>)
  {
    if (not a.is_contiguous() or not a.has_positive_strides())
      NDA_RUNTIME_ERROR << "Error in MPI scatter for nda::Array: Array needs to be contiguous with positive strides";
    return mpi::lazy<mpi::tag::scatter, A>{std::forward<A>(a), comm, root, all};
  }

} // namespace nda
