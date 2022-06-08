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

#include "CompilerPass.hpp"

namespace tket {

const PassPtr &SynthesiseTK();
const PassPtr &SynthesiseTket();
const PassPtr &SynthesiseHQS();
const PassPtr &SynthesiseOQC();
const PassPtr &SynthesiseUMD();

const PassPtr &PeepholeOptimise2Q();
const PassPtr &RemoveRedundancies();
const PassPtr &CommuteThroughMultis();
const PassPtr &DecomposeArbitrarilyControlledGates();
// Expects: CX and any single-qubit gates,
// but does not break if it encounters others
const PassPtr &DecomposeMultiQubitsCX();
const PassPtr &DecomposeSingleQubitsTK1();
const PassPtr &DecomposeBoxes();

/**
 * converts a circuit containing all possible gates to a circuit containing only
 * phase poly boxes + H gates (and measure + reset + collapse + barrier)
 * @param min_size value for the minimal number of CX in each box, groups with
 * less than min_size CX gates are not converted to a PhasePolyBox, default
 * value is 0
 * @return PassPtr to perform the conversion
 */
PassPtr ComposePhasePolyBoxes(unsigned min_size = 0);

/** Squash sequences of single-qubit gates to TK1 gates. */
const PassPtr &SquashTK1();

const PassPtr &RebaseTket();
const PassPtr &RebaseUFR();

const PassPtr &DecomposeBridges();
const PassPtr &FlattenRegisters();

/** Remove all& \ref OpType::Barrier from the circuit. */
const PassPtr &RemoveBarriers();

/** Commutes measurements to the end of the circuit. */
const PassPtr &DelayMeasures();

/**
 * Remove all operations that have no @ref OpType::Output or
 * @ref OpType::ClOutput in their causal future.
 */
const PassPtr &RemoveDiscarded();

/**
 * Replace all measured classical maps that are followed by Measure operations
 * whose quantum output is discarded with classical operations following the
 * Measure.
 */
const PassPtr &SimplifyMeasured();

}  // namespace tket
