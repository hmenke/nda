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

#include "./test_common.hpp"

#include <nda/gtest_tools.hpp>
#include <nda/nda.hpp>
#include <nda/mpi.hpp>

#include <itertools/itertools.hpp>
#include <mpi/mpi.hpp>

#include <numeric>

// Test fixture for testing various algorithms.
struct NDAMpi : public ::testing::Test {
  protected:
  NDAMpi() {
    A.resize(shape_3d);
    std::iota(A.begin(), A.end(), 0);
    M = nda::matrix<std::complex<double>>::zeros(shape_2d);
    for (long i = 0; i < shape_2d[0]; ++i) {
      for (long j = 0; j < shape_2d[1]; ++j) {
        auto x  = static_cast<double>(i * shape_2d[1] + j);
        M(i, j) = std::complex<double>(x, x + 1.0);
      }
    }
    A2 = A;
    M2 = M;
  }

  std::array<long, 3> shape_3d{6, 4, 2};
  nda::array<long, 3> A;
  nda::array<long, 3, nda::basic_layout<0, nda::encode(std::array{1, 2, 0}), nda::layout_prop_e::contiguous>> A2;
  std::array<long, 2> shape_2d{4, 4};
  nda::matrix<std::complex<double>> M;
  nda::matrix<std::complex<double>, nda::F_layout> M2;
  const int root = 0;
  mpi::communicator comm;
};

TEST_F(NDAMpi, ExpectEqualRanks) {
  // test the expect_equal_ranks function
  if (comm.size() > 1) {
    EXPECT_TRUE(nda::detail::have_mpi_equal_ranks(A, comm));
    if (comm.rank() == 1) {
      auto B = nda::array<long, 2>::zeros(6, 4);
      EXPECT_FALSE(nda::detail::have_mpi_equal_ranks(B, comm));
    } else {
      EXPECT_FALSE(nda::detail::have_mpi_equal_ranks(A, comm));
    }
  }
}

TEST_F(NDAMpi, ExpectEqualShapes) {
  // test the expect_equal_shapes function
  if (comm.size() > 1) {
    EXPECT_TRUE(nda::detail::have_mpi_equal_shapes(A, comm));
    if (comm.rank() == 1) A = nda::array<long, 3>::zeros(6, 5, 2);
    EXPECT_FALSE(nda::detail::have_mpi_equal_shapes(A, comm));
  }
}

TEST_F(NDAMpi, ExpectEqualShapeSaveFirst) {
  // test the expect_equal_shape_save_first function
  if (comm.size() > 1) {
    if (comm.rank() == 1) A = nda::array<long, 3>::zeros(5, 4, 2);
    EXPECT_TRUE(nda::detail::have_mpi_equal_shapes(A(nda::range{1}, nda::ellipsis{}), comm));
    if (comm.rank() == 1) A = nda::array<long, 3>::zeros(5, 5, 2);
    EXPECT_FALSE(nda::detail::have_mpi_equal_shapes(A(nda::range{1}, nda::ellipsis{}), comm));
  }
}

TEST_F(NDAMpi, BroadcastCLayout) {
  // broadcast to arrays with same dimensions
  auto A_bcast = A;
  if (comm.rank() == root) {
    EXPECT_ARRAY_EQ(A, A_bcast);
  } else {
    A_bcast = 0;
    EXPECT_ARRAY_ZERO(A_bcast);
  }
  mpi::broadcast(A_bcast, comm, root);
  EXPECT_ARRAY_EQ(A, A_bcast);

  // broadcast to arrays with different dimensions
  decltype(A) B_bcast;
  if (comm.rank() == root) {
    B_bcast = A;
    EXPECT_ARRAY_EQ(A, B_bcast);
  } else {
    EXPECT_NE(A.shape(), B_bcast.shape());
  }
  mpi::broadcast(B_bcast, comm, root);
  EXPECT_ARRAY_EQ(A, B_bcast);

  // broadcast a matrix into an array view
  if (comm.rank() == root) {
    mpi::broadcast(M, comm, root);
  } else {
    auto C_bcast = nda::array<std::complex<double>, 3>::zeros(2, shape_2d[0], shape_2d[1]);
    EXPECT_ARRAY_ZERO(C_bcast);
    auto C_view = C_bcast(1, nda::ellipsis{});
    mpi::broadcast(C_view, comm, root);
    EXPECT_ARRAY_EQ(M, C_view);
    EXPECT_ARRAY_ZERO(C_bcast(0, nda::ellipsis{}));
  }

  // broadcast a view to views
  if (comm.rank() == root) {
    mpi::broadcast(M(0, nda::range::all), comm, root);
  } else {
    auto M_bcast = M;
    M_bcast      = 0;
    mpi::broadcast(M_bcast(0, nda::range::all), comm, root);
    EXPECT_ARRAY_ZERO(M_bcast(nda::range(1, shape_2d[0]), nda::range::all));
    EXPECT_ARRAY_EQ(M_bcast(0, nda::range::all), M(0, nda::range::all));
  }
}

TEST_F(NDAMpi, BroadcastOtherLayouts) {
  // broadcast to arrays with same dimensions
  auto A2_bcast = A2;
  if (comm.rank() == root) {
    EXPECT_ARRAY_EQ(A2, A2_bcast);
  } else {
    A2_bcast = 0;
    EXPECT_ARRAY_ZERO(A2_bcast);
  }
  mpi::broadcast(A2_bcast, comm, root);
  EXPECT_ARRAY_EQ(A2, A2_bcast);

  // broadcast to arrays with different dimensions
  decltype(A2) B2_bcast;
  if (comm.rank() == root) {
    B2_bcast = A2;
    EXPECT_ARRAY_EQ(A2, B2_bcast);
  } else {
    EXPECT_NE(A2.shape(), B2_bcast.shape());
  }
  mpi::broadcast(B2_bcast, comm, root);
  EXPECT_ARRAY_EQ(A2, B2_bcast);

  // broadcast a matrix into an array view
  if (comm.rank() == root) {
    mpi::broadcast(M2, comm, root);
  } else {
    auto C2_bcast = nda::array<std::complex<double>, 3, nda::F_layout>::zeros(shape_2d[0], shape_2d[1], 2);
    EXPECT_ARRAY_ZERO(C2_bcast);
    auto C2_view = C2_bcast(nda::ellipsis{}, 1);
    mpi::broadcast(C2_view, comm, root);
    EXPECT_ARRAY_EQ(M2, C2_view);
    EXPECT_ARRAY_ZERO(C2_bcast(nda::ellipsis{}, 0));
  }

  // broadcast a view to views
  if (comm.rank() == root) {
    mpi::broadcast(M2(nda::range::all, 0), comm, root);
  } else {
    auto M2_bcast = M2;
    M2_bcast      = 0;
    mpi::broadcast(M2_bcast(nda::range::all, 0), comm, root);
    EXPECT_ARRAY_ZERO(M2_bcast(nda::range::all, nda::range(1, shape_2d[1])));
    EXPECT_ARRAY_EQ(M2_bcast(nda::range::all, 0), M2(nda::range::all, 0));
  }

  // broadcast a C layout matrix into an F layout matrix
  if (comm.rank() == root) {
    mpi::broadcast(M, comm, root);
  } else {
    mpi::broadcast(M2, comm, root);
    EXPECT_ARRAY_EQ(nda::transpose(M), M2);
  }
}

TEST_F(NDAMpi, Gather1DArray) {
  // allgather 1-dimensional arrays of different sizes
  auto C        = nda::array<int, 1>(comm.rank() + 1);
  C             = comm.rank();
  auto C_gather = mpi::all_gather(C, comm);
  EXPECT_EQ(C_gather.size(), comm.size() * (comm.size() + 1) / 2);
  for (int i = 0; i < comm.size(); ++i) {
    auto view = C_gather(nda::range(i * (i + 1) / 2, (i + 1) * (i + 2) / 2));
    auto exp  = nda::array<int, 1>(i + 1, i);
    EXPECT_ARRAY_EQ(exp, view);
  }
}

TEST_F(NDAMpi, GatherCLayout) {
  // allgather C-layout arrays
  auto B        = nda::make_regular(A * (comm.rank() + 1));
  auto B_gather = mpi::all_gather(B, comm);
  EXPECT_EQ(B_gather.shape(), (std::array{comm.size() * shape_3d[0], shape_3d[1], shape_3d[2]}));
  for (int i = 0; i < comm.size(); ++i) {
    auto view = B_gather(nda::range(i * shape_3d[0], (i + 1) * shape_3d[0]), nda::range::all, nda::range::all);
    EXPECT_ARRAY_EQ(nda::make_regular(A * (i + 1)), view);
  }

  // gather C-layout matrix views
  auto rg       = nda::range(0, 2);
  decltype(M) N = nda::make_regular(M * (comm.rank() + 1));
  auto N_gather = mpi::gather(N(rg, nda::range::all), comm);
  if (comm.rank() == root) {
    EXPECT_EQ(N_gather.shape(), (std::array<long, 2>{comm.size() * 2l, shape_2d[1]}));
    for (long i = 0; i < comm.size(); ++i) {
      auto view = N_gather(nda::range(i * 2, (i + 1) * 2), nda::range::all);
      auto exp  = nda::make_regular(M * (i + 1));
      EXPECT_ARRAY_EQ(exp(rg, nda::range::all), view);
    }
  } else {
    EXPECT_TRUE(N_gather.empty());
    EXPECT_TRUE(N_gather.size() == 0);
  }
}

TEST_F(NDAMpi, GatherOtherLayouts) {
  // allgather non C-layout arrays by first reshaping it
  constexpr auto perm     = decltype(A2)::layout_t::stride_order;
  constexpr auto inv_perm = nda::permutations::inverse(perm);

  decltype(A2) B2  = nda::make_regular(A2 * (comm.rank() + 1));
  auto B2_gather   = mpi::all_gather(nda::permuted_indices_view<nda::encode(inv_perm)>(B2), comm);
  auto B2_gather_v = nda::permuted_indices_view<nda::encode(perm)>(B2_gather);
  for (int i = 0; i < comm.size(); ++i) {
    auto view = B2_gather_v(nda::range::all, nda::range(i * shape_3d[1], (i + 1) * shape_3d[1]), nda::range::all);
    EXPECT_ARRAY_EQ(nda::make_regular(A2 * (i + 1)), view);
  }
}

TEST_F(NDAMpi, GatherCustomType) {
  // allgather an array of matrices
  using matrix_t = nda::matrix<int>;
  nda::array<matrix_t, 1> B(2);
  nda::array<matrix_t, 1> exp(2 * comm.size());
  for (int i = 0; i < comm.size(); ++i) {
    exp(i * 2)     = matrix_t::zeros(shape_2d);
    exp(i * 2 + 1) = matrix_t::zeros(shape_2d);
    exp(i * 2)     = i;
    exp(i * 2 + 1) = i + 1;
    if (i == comm.rank()) {
      B(0) = exp(i * 2);
      B(1) = exp(i * 2 + 1);
    }
  }

  auto B_gathered = mpi::all_gather(B, comm);
  EXPECT_EQ(B_gathered.shape(), (std::array<long, 1>{2l * comm.size()}));
  for (int i = 0; i < B_gathered.size(); ++i) { EXPECT_ARRAY_EQ(exp(i), B_gathered(i)); }
}

TEST_F(NDAMpi, LazyGather) {
  // lazy-allgather 1-dimensional arrays of different sizes
  auto C               = nda::array<int, 1>(comm.rank() + 1);
  C                    = comm.rank();
  decltype(C) C_gather = nda::lazy_mpi_gather(C, comm, root, true);
  EXPECT_EQ(C_gather.size(), comm.size() * (comm.size() + 1) / 2);
  for (int i = 0; i < comm.size(); ++i) {
    auto view = C_gather(nda::range(i * (i + 1) / 2, (i + 1) * (i + 2) / 2));
    auto exp  = nda::array<int, 1>(i + 1, i);
    EXPECT_ARRAY_EQ(exp, view);
  }
}

TEST_F(NDAMpi, Reduce) {
  // reduce an array
  decltype(A) A_sum = mpi::reduce(A, comm);
  if (comm.rank() == 0) { EXPECT_ARRAY_EQ(nda::make_regular(A * comm.size()), A_sum); }

  // all reduce an array view
  auto B                    = nda::make_regular(A * (comm.rank() + 1));
  nda::array<long, 1> B_max = mpi::all_reduce(B(0, 0, nda::range::all), comm, MPI_MAX);
  nda::array<long, 1> B_min = mpi::all_reduce(B(0, 0, nda::range::all), comm, MPI_MIN);
  EXPECT_ARRAY_EQ(B_max, A(0, 0, nda::range::all) * comm.size());
  EXPECT_ARRAY_EQ(B_min, A(0, 0, nda::range::all));
}

TEST_F(NDAMpi, ReduceCustomType) {
  using namespace nda::clef::literals;

  // reduce an array of matrices
  using matrix_t = nda::matrix<double>;
  nda::array<matrix_t, 1> B(7);
  nda::array<matrix_t, 1> exp_sum(7);

  for (int i = 0; i < B.extent(0); ++i) {
    B(i) = matrix_t{4, 4};
    B(i)(k_, l_) << i * (comm.rank() + 1) * (k_ + l_);

    exp_sum(i) = matrix_t{4, 4};
    exp_sum(i)(k_, l_) << i * (comm.size() + 1) * comm.size() / 2 * (k_ + l_);
  }

  nda::array<matrix_t, 1> B_sum = mpi::all_reduce(B, comm);

  EXPECT_ARRAY_EQ(B_sum, exp_sum);
}

TEST_F(NDAMpi, Scatter) {
  // scatter an array
  decltype(A) A_scatter = mpi::scatter(A, comm);
  auto chunked_rg       = itertools::chunk_range(0, A.shape()[0], comm.size(), comm.rank());
  auto exp_shape        = A.shape();
  exp_shape[0]          = chunked_rg.second - chunked_rg.first;
  EXPECT_EQ(exp_shape, A_scatter.shape());
  EXPECT_ARRAY_EQ(A(nda::range(chunked_rg.first, chunked_rg.second), nda::ellipsis{}), A_scatter);
}

TEST_F(NDAMpi, BroadcastTransposedMatrix) {
  nda::matrix<std::complex<double>> M_t = transpose(M);
  nda::matrix<std::complex<double>> N;
  if (comm.rank() == 0) N = M_t;
  mpi::broadcast(N, comm, 0);
  EXPECT_ARRAY_EQ(M_t, N);
}

TEST_F(NDAMpi, BroadcastTransposedArray) {
  nda::array<long, 3> A_t = transpose(A);
  nda::array<long, 3> B(2, 4, 6);
  if (comm.rank() == 0) B = A_t;
  mpi::broadcast(B, comm, 0);
  EXPECT_ARRAY_EQ(A_t, B);
}

TEST_F(NDAMpi, VariousCollectiveCommunications) {
  using arr_t = nda::array<std::complex<double>, 2>;

  arr_t A(7, 3);
  for (int i = 0; i < A.extent(0); ++i)
    for (int j = 0; j < A.extent(1); ++j) A(i, j) = i + 10 * j;

  // scatter an array
  arr_t B;
  B       = mpi::scatter(A, comm);
  arr_t C = mpi::scatter(A, comm);
  auto rg = itertools::chunk_range(0, 7, comm.size(), comm.rank());
  EXPECT_ARRAY_EQ(B, A(nda::range(rg.first, rg.second), nda::range::all));
  EXPECT_ARRAY_NEAR(C, B);

  // gather an array
  B *= -1;
  arr_t D = mpi::gather(B, comm);
  if (comm.rank() == 0) { EXPECT_ARRAY_NEAR(D, -A); }

  // broadcast an array
  mpi::broadcast(D, comm);
  EXPECT_ARRAY_NEAR(D, -A);

  // all gather an array
  D() = 0;
  D   = mpi::all_gather(B, comm);
  EXPECT_ARRAY_NEAR(D, -A);

  // reduce an array
  arr_t R1 = mpi::reduce(A, comm);
  if (comm.rank() == 0) { EXPECT_ARRAY_NEAR(R1, comm.size() * A); }

  // all reduce an array
  arr_t R2 = mpi::all_reduce(A, comm);
  EXPECT_ARRAY_NEAR(R2, comm.size() * A);
}

MPI_TEST_MAIN
