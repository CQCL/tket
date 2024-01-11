// Copyright 2019-2024 Cambridge Quantum Computing
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
 * @param allow_partial Whether to allow measurements that cannot be commuted to
 * the end, and delay them as much as possible instead.
 * @throws CircuitInvalidity if it is not possible to delay a measurement, and
 * \p allow_partial is false.
 */
Transform delay_measures(bool allow_partial = false);

namespace DelayMeasures {

/** Commute all measurement gates to the end of the circuit.
 * @param circ The circuit to delay measurements in.
 * @param allow_partial Whether to allow measurements that cannot be commuted to
 * the end, and delay them as much as possible instead.
 * @param dry_run If true, do not modify the circuit, just check if it is
 * possible to delay.
 * @throws CircuitInvalidity if it is not possible to delay and both \p
 * allow_partial and \p dry_run are false.
 * @return A pair of booleans. The first indicates when the circuit was changed,
 * and the second indicates if the run found no errors (i.e. it was possible to
 * delay all measures to the end, or \p allow_partial was true).
 **/
std::pair<bool, bool> run_delay_measures(
    Circuit& circ, bool allow_partial = false, bool dry_run = false);

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
