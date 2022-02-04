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

#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index_container.hpp>

namespace tket {

struct TagKey {};
struct TagValue {};
struct TagSeq {};

template <typename A, typename B>
using sequenced_bimap_t = boost::multi_index::multi_index_container<
    std::pair<A, B>,
    boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
            boost::multi_index::tag<TagKey>,
            boost::multi_index::member<
                std::pair<A, B>, A, &std::pair<A, B>::first>>,
        boost::multi_index::ordered_unique<
            boost::multi_index::tag<TagValue>,
            boost::multi_index::member<
                std::pair<A, B>, B, &std::pair<A, B>::second>>,
        boost::multi_index::sequenced<boost::multi_index::tag<TagSeq>>>>;

template <typename A, typename B>
using sequenced_map_t = boost::multi_index::multi_index_container<
    std::pair<A, B>,
    boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
            boost::multi_index::tag<TagKey>,
            boost::multi_index::member<
                std::pair<A, B>, A, &std::pair<A, B>::first>>,
        boost::multi_index::sequenced<boost::multi_index::tag<TagSeq>>>>;

template <typename T>
using sequence_set_t = boost::multi_index::multi_index_container<
    T,
    boost::multi_index::indexed_by<
        boost::multi_index::ordered_unique<
            boost::multi_index::tag<TagKey>, boost::multi_index::identity<T>>,
        boost::multi_index::sequenced<boost::multi_index::tag<TagSeq>>>>;

}  // namespace tket
