// Copyright 2019-2022 Cambridge Quantum Computing
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

#pragma once

#if !defined(_MSC_VER)
#pragma GCC diagnostic push

#if defined(__clang__)
#pragma GCC diagnostic ignored "-Wunknown-warning-option"
#if __has_warning("-Wdeprecated-copy")
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif
#else
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif

#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wuninitialized"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#if !defined(_MSC_VER)
#pragma GCC diagnostic pop
#endif
