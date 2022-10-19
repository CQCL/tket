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

#include <boost/range/join.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>

#include "Architecture/Architecture.hpp"
#include "Circuit/CircPool.hpp"
#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Circuit/Command.hpp"
#include "CircuitsForTesting.hpp"
#include "Converters/PhasePoly.hpp"
#include "Gate/SymTable.hpp"
#include "Mapping/LexiLabelling.hpp"
#include "Mapping/LexiRoute.hpp"
#include "Mapping/RoutingMethod.hpp"
#include "MeasurementSetup/MeasurementSetup.hpp"
#include "OpType/OpType.hpp"
#include "Ops/OpPtr.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Predicates/PassLibrary.hpp"
#include "Transformations/OptimisationPass.hpp"
#include "Transformations/PauliOptimisation.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/Json.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_json {

template <class T>
bool serialize_deserialize(const T& obj) {
  nlohmann::json j = obj;
  auto new_obj = j.get<T>();
  return obj == new_obj;
}

template <class T>
void check_cases(const std::vector<T>& cases) {
  for (const auto& test : cases) {
    CHECK(serialize_deserialize(test));
  }
}

bool check_circuit(const Circuit& c) {
  nlohmann::json j = c;
  auto new_c = j.get<Circuit>();
  return c.circuit_equality(new_c);
}


SCENARIO("TEST PROBLEM"){

  Architecture arc = SquareGrid(2, 4);
  Placement::Ptr la_place = std::make_shared<LinePlacement>(arc);

    Circuit circ = CircuitsForTesting::get().uccsd;    
    CompilationUnit cu{circ};                          
    CompilationUnit copy = cu;                         
    PassPtr pp = gen_placement_pass(la_place);                                 
    nlohmann::json j_pp = pp;                          
    PassPtr loaded = j_pp.get<PassPtr>();              
    pp->apply(cu);                                     
    loaded->apply(copy);                               
    REQUIRE(cu.get_circ_ref() == copy.get_circ_ref()); 
    nlohmann::json j_loaded = loaded;         
    std::cout << "JSON: " << j_loaded << std::endl;         
    REQUIRE(j_pp == j_loaded);                         
}



}
}
I'm no necessarily shocked by this result. Frame randomisation is designed to spin out coherent errors to make them stochastic, rather than reduce the errors. As such I would expect the errors in the 11111 and 00000 strings to become comparable but would not necessarily expect the total error in both strings to reduce. Have. I understood that correctly 
@Silas Dilkes
? I've not used frame randomisation on the IBM devices so I can't say for sure if that's what's happening.