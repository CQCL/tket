// Copyright Quantinuum
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Adapted from the corresponding files in:
// - https://github.com/pybind/pybind11_json
// - https://github.com/ianhbell/nanobind_json
// both of which are licensed under the BSD 3-Clause license reproduced below.

// BSD 3-Clause License

// Copyright (c) 2019,
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.

// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.

// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>

#include <nlohmann/json.hpp>
#include <string>

namespace nb = nanobind;
namespace nl = nlohmann;

namespace pyjson {
inline nb::object from_json(const nl::json& j) {
  if (j.is_null()) {
    return nb::none();
  } else if (j.is_boolean()) {
    return nb::bool_(j.get<bool>());
  } else if (j.is_number_unsigned()) {
    return nb::int_(j.get<nl::json::number_unsigned_t>());
  } else if (j.is_number_integer()) {
    return nb::int_(j.get<nl::json::number_integer_t>());
  } else if (j.is_number_float()) {
    return nb::float_(j.get<double>());
  } else if (j.is_string()) {
    const std::string s = j.get<std::string>();
    return nb::str(s.c_str(), s.size());
  } else if (j.is_array()) {
    nb::list obj;
    for (const auto& el : j) {
      obj.append(from_json(el));
    }
    return obj;
  } else  // Object
  {
    nb::dict obj;
    for (nl::json::const_iterator it = j.cbegin(); it != j.cend(); ++it) {
      obj[nb::str(it.key().c_str())] = from_json(it.value());
    }
    return obj;
  }
}

inline nl::json to_json(const nb::handle& obj) {
  if (obj.ptr() == nullptr || obj.is_none()) {
    return nullptr;
  }
  if (nb::isinstance<nb::bool_>(obj)) {
    return nb::cast<bool>(obj);
  }
  if (nb::isinstance<nb::int_>(obj)) {
    try {
      nl::json::number_integer_t s = nb::cast<nl::json::number_integer_t>(obj);
      if (nb::int_(s).equal(obj)) {
        return s;
      }
    } catch (...) {
    }
    try {
      nl::json::number_unsigned_t u =
          nb::cast<nl::json::number_unsigned_t>(obj);
      if (nb::int_(u).equal(obj)) {
        return u;
      }
    } catch (...) {
    }
    throw std::runtime_error(
        "to_json received an integer out of range for both "
        "nl::json::number_integer_t and nl::json::number_unsigned_t type: " +
        nb::cast<std::string>(nb::repr(obj)));
  }
  if (nb::isinstance<nb::float_>(obj)) return nb::cast<double>(obj);

  if (nb::isinstance<nb::bytes>(obj)) {
    nb::module_ base64 = nb::module_::import_("base64");
    return nb::cast<std::string>(
        base64.attr("b64encode")(obj).attr("decode")("utf-8"));
  }
  if (nb::isinstance<nb::str>(obj)) {
    return nb::cast<std::string>(obj);
  }
  if (nb::isinstance<nb::tuple>(obj) || nb::isinstance<nb::list>(obj)) {
    auto out = nl::json::array();
    for (const nb::handle value : obj) {
      out.push_back(to_json(value));
    }
    return out;
  }
  if (nb::isinstance<nb::dict>(obj)) {
    auto out = nl::json::object();
    for (const nb::handle key : obj) {
      out[nb::cast<std::string>(nb::str(key))] = to_json(obj[key]);
    }
    return out;
  }
  throw std::runtime_error(
      "to_json not implemented for this type of object: " +
      nb::cast<std::string>(nb::repr(obj)));
}
}  // namespace pyjson

// nlohmann_json serializers
namespace nlohmann {
#define MAKE_NLJSON_SERIALIZER_DESERIALIZER(T)          \
  template <>                                           \
  struct adl_serializer<T> {                            \
    inline static void to_json(json& j, const T& obj) { \
      j = pyjson::to_json(obj);                         \
    }                                                   \
                                                        \
    inline static T from_json(const json& j) {          \
      return nb::cast<T>(pyjson::from_json(j));         \
    }                                                   \
  };

#define MAKE_NLJSON_SERIALIZER_ONLY(T)                  \
  template <>                                           \
  struct adl_serializer<T> {                            \
    inline static void to_json(json& j, const T& obj) { \
      j = pyjson::to_json(obj);                         \
    }                                                   \
  };

MAKE_NLJSON_SERIALIZER_DESERIALIZER(nb::object);

MAKE_NLJSON_SERIALIZER_DESERIALIZER(nb::bool_);
MAKE_NLJSON_SERIALIZER_DESERIALIZER(nb::int_);
MAKE_NLJSON_SERIALIZER_DESERIALIZER(nb::float_);
MAKE_NLJSON_SERIALIZER_DESERIALIZER(nb::str);

MAKE_NLJSON_SERIALIZER_DESERIALIZER(nb::list);
MAKE_NLJSON_SERIALIZER_DESERIALIZER(nb::tuple);
MAKE_NLJSON_SERIALIZER_DESERIALIZER(nb::dict);

MAKE_NLJSON_SERIALIZER_ONLY(nb::handle);
// MAKE_NLJSON_SERIALIZER_ONLY(nb::detail::item_accessor);
// MAKE_NLJSON_SERIALIZER_ONLY(nb::detail::list_accessor);
// MAKE_NLJSON_SERIALIZER_ONLY(nb::detail::tuple_accessor);
// MAKE_NLJSON_SERIALIZER_ONLY(nb::detail::sequence_accessor);
// MAKE_NLJSON_SERIALIZER_ONLY(nb::detail::str_attr_accessor);
// MAKE_NLJSON_SERIALIZER_ONLY(nb::detail::obj_attr_accessor);

#undef MAKE_NLJSON_SERIALIZER_DESERIALIZER
#undef MAKE_NLJSON_SERIALIZER_ONLY
}  // namespace nlohmann

// nanobind caster
namespace nanobind {
namespace detail {
template <>
struct type_caster<nl::json> {
 public:
  NB_TYPE_CASTER(nl::json, const_name("JSON"))

  bool from_python(handle src, uint8_t, cleanup_list*) noexcept {
    try {
      auto value = pyjson::to_json(src);
      return true;
    } catch (...) {
      return false;
    }
  }

  static handle from_cpp(nl::json src, rv_policy, cleanup_list*) noexcept {
    object obj = pyjson::from_json(src);
    return obj.release();
  }
};
}  // namespace detail
}  // namespace nanobind
