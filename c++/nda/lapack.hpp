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
 * @brief Includes all LAPACK relevant headers.
 */

#pragma once

#include "./lapack/interface/cxx_interface.hpp"
#include "./lapack/gelss.hpp"
#include "./lapack/geqp3.hpp"
#include "./lapack/gesvd.hpp"
#include "./lapack/getrf.hpp"
#include "./lapack/getri.hpp"
#include "./lapack/getrs.hpp"
#include "./lapack/gtsv.hpp"
#include "./lapack/orgqr.hpp"
#include "./lapack/ungqr.hpp"
