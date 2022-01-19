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

/**
 * @file
 * @brief Generally useful typedefs and constants
 */

#include <complex>

namespace tket {

// Typedefs

/** Complex number */
typedef std::complex<double> Complex;

// Numeric constants
/** A fixed square root of -1 */
constexpr Complex i_(0, 1);
/** Complex zero */
constexpr Complex czero(0, 0);

/** Default tolerance for floating-point comparisons */
constexpr double EPS = 1e-11;

/** \f$ \pi \f$ */
constexpr double PI =
    3.141'592'653'589'793'238'462'643'383'279'502'884'197'169'399'375'105'820'974;

/**
 * Enum for readout basis and ordering
 *
 * Readouts are viewed in increasing lexicographic order (ILO) of the bit's
 * UnitID. This is our default convention for column indexing for ALL readout
 * forms (shots, counts, statevector, and unitaries). e.g. |abc> corresponds to
 * the readout: ('c', 0) --> a, ('c', 1) --> b, ('d', 0) --> c
 *
 * For statevector and unitaries, the string abc is interpreted as an index
 * in a big-endian (BE) fashion.
 * e.g. the statevector (a_{00}, a_{01}, a_{10}, a_{11})
 *
 * Some backends (Qiskit, ProjectQ, Quest, etc.) use a DLO-BE
 * (decreasing lexicographic order, big-endian) convention. This is the same
 * as ILO-LE (little-endian) for statevectors and unitaries, but gives shot
 * tables/readouts in a counter-intuitive manner.
 *
 * Every backend and matrix-based box should have a BasisOrder option which
 * can toggle between ILO-BE (ilo) and DLO-BE (dlo).
 */
enum class BasisOrder {
  ilo, /**< increasing lexicographic order (big-endian) */
  dlo  /**< decreasing lexicographic order (big-endian) */
};

}  // namespace tket
