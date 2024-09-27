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
 * @brief Defines accessors for nda::array objects (cf. std::default_accessor).
 */

#pragma once

#include "./macros.hpp"

#include <cstddef>

namespace nda {

  /**
   * @addtogroup av_utils
   * @{
   */

  /// Default accessor for various array and view types.
  struct default_accessor {
    /**
     * @brief Accessor type of the nda::default_accessor.
     * @tparam T Value type of the data.
     */
    template <typename T>
    struct accessor {
      /// Value type of the data.
      using element_type = T;

      /// Pointer type to the data.
      using pointer = T *;

      /// Reference type to the data.
      using reference = T &;

      /**
       * @brief Access a specific element of the data.
       *
       * @param p Pointer to the data.
       * @param i Index of the element to access.
       * @return Reference to the element.
       */
      FORCEINLINE static reference access(pointer p, std::ptrdiff_t i) noexcept {
        EXPECTS(p != nullptr);
        return p[i];
      }

      /**
       * @brief Offset the pointer by a certain number of elements.
       *
       * @param p Pointer to the data.
       * @param i Number of elements to offset the pointer by.
       * @return Pointer after applying the offset.
       */
      FORCEINLINE static T *offset(pointer p, std::ptrdiff_t i) noexcept { return p + i; }
    };
  };

  /// Accessor for array and view types with no aliasing.
  struct no_alias_accessor {
    /**
     * @brief Accessor type of the nda::no_alias_accessor.
     * @tparam T Value type of the data.
     */
    template <typename T>
    struct accessor {
      /// Value type of the data.
      using element_type = T;

      /// Restricted pointer type to the data.
      using pointer = T *__restrict;

      /// Reference type to the data.
      using reference = T &;

      /**
       * @brief Access a specific element of the data.
       *
       * @param p Pointer to the data.
       * @param i Index of the element to access.
       * @return Reference to the element.
       */
      FORCEINLINE static reference access(pointer p, std::ptrdiff_t i) noexcept { return p[i]; }

      /**
       * @brief Offset the pointer by a certain number of elements.
       *
       * @param p Pointer to the data.
       * @param i Number of elements to offset the pointer by.
       * @return Pointer after applying the offset.
       */
      FORCEINLINE static T *offset(pointer p, std::ptrdiff_t i) noexcept { return p + i; }
    };
  };

  /** @} */

} // namespace nda
