// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "Transform.hpp"

namespace tket {

namespace Transforms {

/**
 * Commute all measurement gates to the end of the circuit.
 * Throws a CircuitInvalidity exception if it is not possible to delay.
 */
Transform delay_measures();

namespace DelayMeasures {

/** Commute all measurement gates to the end of the circuit.
 * @param circ The circuit to delay measurements in.
 * @param dry_run If true, do not modify the circuit, just check if it is
 * possible to delay.
 * @throws CircuitInvalidity if it is not possible to delay and dry_run is
 * false.
 * @return A pair of booleans. The first indicates whether the circuit was
 * changed, and the second indicates whether it was possible to delay (i.e.
 * there where no errors).
 **/
std::pair<bool, bool> run_delay_measures(Circuit& circ, bool dry_run = false);

/**
 * Gathers all end-measurements, and adds the measured units to the list.
 * Rejects gates acting on measured_units, and terminates early.
 * Applies recursively for CircBoxes and Conditionals.
 * @param com The command to check.
 * @param measured_units The list of measured units to add to, initialy
 * populated with previously-measured units.
 * @return Whether there are no mid-circuit measurements.
 **/
bool check_only_end_measures(const Command& com, unit_set_t& measured_units);

}  // namespace DelayMeasures

}  // namespace Transforms

}  // namespace tket
