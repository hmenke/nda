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
 * @brief Provides an MPI reduce function for nda::Array types.
 */

#pragma once

#include "../basic_functions.hpp"
#include "../concepts.hpp"
#include "../exceptions.hpp"
#include "../map.hpp"
#include "../traits.hpp"

#include <mpi/mpi.hpp>

#include <cstdlib>
#include <type_traits>
#include <utility>

/**
 * @ingroup av_mpi
 * @brief Specialization of the `mpi::lazy` class for nda::Array types and the `mpi::tag::reduce` tag.
 *
 * @details An object of this class is returned when reducing nda::Array objects across multiple MPI processes.
 *
 * It models an nda::ArrayInitializer, that means it can be used to initialize and assign to nda::basic_array and
 * nda::basic_array_view objects. The target array will have the same shape as the input arrays.
 *
 * See nda::mpi_reduce for an example.
 *
 * @tparam A nda::Array type to be reduced.
 */
template <nda::Array A>
struct mpi::lazy<mpi::tag::reduce, A> {
  /// Value type of the array/view.
  using value_type = typename std::decay_t<A>::value_type;

  /// Const view type of the array/view stored in the lazy object.
  using const_view_type = decltype(std::declval<const A>()());

  /// View of the array/view to be reduced.
  const_view_type rhs;

  /// MPI communicator.
  mpi::communicator comm;

  /// MPI root process.
  const int root{0}; // NOLINT (const is fine here)

  /// Should all processes receive the result.
  const bool all{false}; // NOLINT (const is fine here)

  /// MPI reduction operation.
  const MPI_Op op{MPI_SUM}; // NOLINT (const is fine here)

  /**
   * @brief Compute the shape of the target array.
   * @note It is assumed that the shape of the input array is the same for all MPI processes.
   * @return Shape of the input array.
   */
  [[nodiscard]] auto shape() const { return rhs.shape(); }

  /**
   * @brief Execute the lazy MPI operation and write the result to a target array/view.
   *
   * @details If the target array/view is the same as the input array/view, i.e. if their data pointers are the same,
   * the reduction is performed in-place.
   *
   * @tparam T nda::Array type of the target array/view.
   * @param target Target array/view.
   */
  template <nda::Array T>
  void invoke(T &&target) const { // NOLINT (temporary views are allowed here)
    // check if the arrays can be used in the MPI call
    if (not target.is_contiguous() or not target.has_positive_strides())
      NDA_RUNTIME_ERROR << "Error in MPI reduce for nda::Array: Target array needs to be contiguous with positive strides";

    static_assert(std::decay_t<A>::layout_t::stride_order_encoded == std::decay_t<T>::layout_t::stride_order_encoded,
                  "Error in MPI reduce for nda::Array: Incompatible stride orders");

    // special case for non-mpi runs
    if (not mpi::has_env) {
      target = rhs;
      return;
    }

    // perform the reduction
    if constexpr (not mpi::has_mpi_type<value_type>) {
      // if the value type cannot be reduced directly, we call mpi::reduce for each element
      target = nda::map([this](auto const &x) { return mpi::reduce(x, this->comm, this->root, this->all, this->op); })(rhs);
    } else {
      // value type has a corresponding MPI type
      bool in_place = (target.data() == rhs.data());
      if (in_place) {
        if (rhs.size() != target.size())
          NDA_RUNTIME_ERROR << "Error in MPI reduce for nda::Array: In-place reduction requires arrays of the same size";
      } else {
        if ((comm.rank() == root) || all) nda::resize_or_check_if_view(target, shape());
        if (std::abs(target.data() - rhs.data()) < rhs.size()) NDA_RUNTIME_ERROR << "Error in MPI reduce for nda::Array: Overlapping arrays";
      }

      void *target_ptr    = (void *)target.data();
      void *rhs_ptr       = (void *)rhs.data();
      auto count          = rhs.size();
      auto mpi_value_type = mpi::mpi_type<value_type>::get();
      if (!all) {
        if (in_place)
          MPI_Reduce((comm.rank() == root ? MPI_IN_PLACE : rhs_ptr), rhs_ptr, count, mpi_value_type, op, root, comm.get());
        else
          MPI_Reduce(rhs_ptr, target_ptr, count, mpi_value_type, op, root, comm.get());
      } else {
        if (in_place)
          MPI_Allreduce(MPI_IN_PLACE, rhs_ptr, count, mpi_value_type, op, comm.get());
        else
          MPI_Allreduce(rhs_ptr, target_ptr, count, mpi_value_type, op, comm.get());
      }
    }
  }
};

namespace nda {

  /**
   * @ingroup av_mpi
   * @brief Implementation of an MPI reduce for nda::basic_array or nda::basic_array_view types.
   *
   * @details Since the returned `mpi::lazy` object models an nda::ArrayInitializer, it can be used to initialize/assign
   * to nda::basic_array and nda::basic_array_view objects:
   *
   * @code{.cpp}
   * // create an array on all processes
   * nda::array<int, 2> arr(3, 4);
   *
   * // ...
   * // fill array on each process
   * // ...
   *
   * // reduce the array to the root process
   * nda::array<int, 2> res = mpi::reduce(arr);
   * @endcode
   *
   * @tparam A nda::basic_array or nda::basic_array_view type.
   * @param a Array or view to be gathered.
   * @param comm `mpi::communicator` object.
   * @param root Rank of the root process.
   * @param all Should all processes receive the result of the gather.
   * @param op MPI reduction operation.
   * @return An `mpi::lazy` object modelling an nda::ArrayInitializer.
   */
  template <typename A>
  ArrayInitializer<std::remove_reference_t<A>> auto mpi_reduce(A &&a, mpi::communicator comm = {}, int root = 0, bool all = false,
                                                               MPI_Op op = MPI_SUM)
    requires(is_regular_or_view_v<A>)
  {
    if (not a.is_contiguous() or not a.has_positive_strides())
      NDA_RUNTIME_ERROR << "Error in MPI reduce for nda::Array: Array needs to be contiguous with positive strides";
    return mpi::lazy<mpi::tag::reduce, A>{std::forward<A>(a), comm, root, all, op};
  }

} // namespace nda
