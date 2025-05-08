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

#include "tket/Circuit/PauliExpBoxes.hpp"

#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/ConjugationBox.hpp"
#include "tket/Converters/PhasePoly.hpp"
#include "tket/Diagonalisation/Diagonalisation.hpp"
#include "tket/Ops/OpJsonFactory.hpp"
#include "tket/Transformations/CliffordOptimisation.hpp"
#include "tket/Transformations/GreedyPauliOptimisation.hpp"

namespace tket {

PauliExpBox::PauliExpBox(
    const SymPauliTensor &paulis, CXConfigType cx_config_type)
    : Box(OpType::PauliExpBox,
          op_signature_t(paulis.size(), EdgeType::Quantum)),
      paulis_(paulis),
      cx_config_(cx_config_type) {}

PauliExpBox::PauliExpBox(const PauliExpBox &other)
    : Box(other), paulis_(other.paulis_), cx_config_(other.cx_config_) {}

PauliExpBox::PauliExpBox() : PauliExpBox({{}, 0}) {}

bool PauliExpBox::is_clifford() const {
  return equiv_0(4 * paulis_.coeff) || paulis_.string.empty();
}

SymSet PauliExpBox::free_symbols() const { return paulis_.free_symbols(); }

Op_ptr PauliExpBox::dagger() const {
  return std::make_shared<PauliExpBox>(
      SymPauliTensor(paulis_.string, -paulis_.coeff), cx_config_);
}

Op_ptr PauliExpBox::transpose() const {
  SymPauliTensor tr = paulis_;
  tr.transpose();
  return std::make_shared<PauliExpBox>(tr, cx_config_);
}

Op_ptr PauliExpBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  return std::make_shared<PauliExpBox>(
      this->paulis_.symbol_substitution(sub_map), this->cx_config_);
}

void PauliExpBox::generate_circuit() const {
  // paulis_ gets cast to a sparse form, so circuit from pauli_gadget will only
  // contain qubits with {X, Y, Z}; appending it to a blank circuit containing
  // all qubits makes the size of the circuit fixed
  Circuit circ(paulis_.size());
  circ.append(pauli_gadget(paulis_, cx_config_));
  circ_ = std::make_shared<Circuit>(circ);
}

bool PauliExpBox::is_equal(const Op &op_other) const {
  const PauliExpBox &other = dynamic_cast<const PauliExpBox &>(op_other);
  if (id_ == other.get_id()) return true;
  return cx_config_ == other.cx_config_ && paulis_.equiv_mod(other.paulis_, 4);
}

nlohmann::json PauliExpBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const PauliExpBox &>(*op);
  nlohmann::json j = core_box_json(box);
  // Serialise paulis and phase separately for backwards compatibility with
  // before templated PauliTensor
  j["paulis"] = box.get_paulis();
  j["phase"] = box.get_phase();
  j["cx_config"] = box.get_cx_config();
  return j;
}

Op_ptr PauliExpBox::from_json(const nlohmann::json &j) {
  PauliExpBox box = PauliExpBox(
      SymPauliTensor(
          j.at("paulis").get<std::vector<Pauli>>(), j.at("phase").get<Expr>()),
      j.at("cx_config").get<CXConfigType>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(PauliExpBox, PauliExpBox)

PauliExpPairBox::PauliExpPairBox(
    const SymPauliTensor &paulis0, const SymPauliTensor &paulis1,
    CXConfigType cx_config_type)
    : Box(OpType::PauliExpPairBox,
          op_signature_t(paulis0.size(), EdgeType::Quantum)),
      paulis0_(paulis0),
      paulis1_(paulis1),
      cx_config_(cx_config_type) {
  if (paulis0.size() != paulis1.size()) {
    throw PauliExpBoxInvalidity(
        "Pauli strings within PauliExpPairBox must be of same length (pad with "
        "identities if necessary)");
  }
}

PauliExpPairBox::PauliExpPairBox(const PauliExpPairBox &other)
    : Box(other),
      paulis0_(other.paulis0_),
      paulis1_(other.paulis1_),
      cx_config_(other.cx_config_) {}

PauliExpPairBox::PauliExpPairBox() : PauliExpPairBox({{}, 0}, {{}, 0}) {}

bool PauliExpPairBox::is_clifford() const {
  auto is_clifford0 = equiv_0(4 * paulis0_.coeff) || paulis0_.string.empty();
  auto is_clifford1 = equiv_0(4 * paulis1_.coeff) || paulis1_.string.empty();
  return is_clifford0 && is_clifford1;
}

SymSet PauliExpPairBox::free_symbols() const {
  return expr_free_symbols({paulis0_.coeff, paulis1_.coeff});
}

Op_ptr PauliExpPairBox::dagger() const {
  return std::make_shared<PauliExpPairBox>(
      SymPauliTensor(paulis1_.string, -paulis1_.coeff),
      SymPauliTensor(paulis0_.string, -paulis0_.coeff), cx_config_);
}

Op_ptr PauliExpPairBox::transpose() const {
  SymPauliTensor tr0 = paulis0_;
  tr0.transpose();
  SymPauliTensor tr1 = paulis1_;
  tr1.transpose();
  return std::make_shared<PauliExpPairBox>(tr1, tr0, cx_config_);
}

Op_ptr PauliExpPairBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  return std::make_shared<PauliExpPairBox>(
      this->paulis0_.symbol_substitution(sub_map),
      this->paulis1_.symbol_substitution(sub_map), this->cx_config_);
}

void PauliExpPairBox::generate_circuit() const {
  // paulis0_ and paulis1_ gets cast to a sparse form, so circuit from
  // pauli_gadget_pair will only contain qubits with {X, Y, Z} on at least one;
  // appending it to a blank circuit containing all qubits makes the size of the
  // circuit fixed
  Circuit circ(paulis0_.size());
  circ.append(pauli_gadget_pair(paulis0_, paulis1_, cx_config_));
  circ_ = std::make_shared<Circuit>(circ);
}

bool PauliExpPairBox::is_equal(const Op &op_other) const {
  const PauliExpPairBox &other =
      dynamic_cast<const PauliExpPairBox &>(op_other);
  if (id_ == other.get_id()) return true;
  return cx_config_ == other.cx_config_ &&
         paulis0_.equiv_mod(other.paulis0_, 4) &&
         paulis1_.equiv_mod(other.paulis1_, 4);
}

nlohmann::json PauliExpPairBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const PauliExpPairBox &>(*op);
  nlohmann::json j = core_box_json(box);
  // Encode pauli strings and phases separately for backwards compatibility from
  // before templated PauliTensor
  auto paulis_pair = box.get_paulis_pair();
  // use vector to avoid serialising into a dictionary if the Pauli strings are
  // of length 2
  std::vector<std::vector<Pauli>> paulis_vec{
      paulis_pair.first, paulis_pair.second};
  j["paulis_pair"] = paulis_vec;
  j["phase_pair"] = box.get_phase_pair();
  j["cx_config"] = box.get_cx_config();
  return j;
}

Op_ptr PauliExpPairBox::from_json(const nlohmann::json &j) {
  const auto [paulis0, paulis1] =
      j.at("paulis_pair")
          .get<std::tuple<std::vector<Pauli>, std::vector<Pauli>>>();
  const auto [phase0, phase1] =
      j.at("phase_pair").get<std::tuple<Expr, Expr>>();
  PauliExpPairBox box = PauliExpPairBox(
      SymPauliTensor(paulis0, phase0), SymPauliTensor(paulis1, phase1),
      j.at("cx_config").get<CXConfigType>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(PauliExpPairBox, PauliExpPairBox)

PauliExpCommutingSetBox::PauliExpCommutingSetBox(
    const std::vector<SymPauliTensor> &pauli_gadgets,
    CXConfigType cx_config_type)
    : Box(OpType::PauliExpCommutingSetBox),
      pauli_gadgets_(pauli_gadgets),
      cx_config_(cx_config_type) {
  // check at least one gadget
  if (pauli_gadgets.empty()) {
    throw PauliExpBoxInvalidity(
        "PauliExpCommutingSetBox requires at least one Pauli string");
  }
  // check all gadgets have same Pauli string length
  auto n_qubits = pauli_gadgets[0].size();
  for (const auto &gadget : pauli_gadgets) {
    if (gadget.size() != n_qubits) {
      throw PauliExpBoxInvalidity(
          "the Pauli strings within PauliExpCommutingSetBox must all be the "
          "same length");
    }
  }
  if (!this->paulis_commute()) {
    throw PauliExpBoxInvalidity(
        "Pauli strings used to define PauliExpCommutingSetBox must all "
        "commute");
  }
  signature_ = op_signature_t(n_qubits, EdgeType::Quantum);
}

PauliExpCommutingSetBox::PauliExpCommutingSetBox(
    const PauliExpCommutingSetBox &other)
    : Box(other),
      pauli_gadgets_(other.pauli_gadgets_),
      cx_config_(other.cx_config_) {}

PauliExpCommutingSetBox::PauliExpCommutingSetBox()
    : PauliExpCommutingSetBox({{{}, 0}}) {}

bool PauliExpCommutingSetBox::is_clifford() const {
  return std::all_of(
      pauli_gadgets_.begin(), pauli_gadgets_.end(),
      [](const SymPauliTensor &pauli_exp) {
        return equiv_0(4 * pauli_exp.coeff) || pauli_exp.string.empty();
      });
}

bool PauliExpCommutingSetBox::paulis_commute() const {
  for (auto string0 = pauli_gadgets_.begin(); string0 != pauli_gadgets_.end();
       ++string0) {
    for (auto string1 = std::next(string0); string1 != pauli_gadgets_.end();
         ++string1) {
      if (!string0->commutes_with(*string1)) {
        return false;
      }
    }
  }
  return true;
}

SymSet PauliExpCommutingSetBox::free_symbols() const {
  std::vector<Expr> angles;
  for (const auto &pauli_exp : pauli_gadgets_) {
    angles.push_back(pauli_exp.coeff);
  }
  return expr_free_symbols(angles);
}

Op_ptr PauliExpCommutingSetBox::dagger() const {
  std::vector<SymPauliTensor> dagger_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    dagger_gadgets.emplace_back(pauli_exp.string, -pauli_exp.coeff);
  }
  return std::make_shared<PauliExpCommutingSetBox>(dagger_gadgets, cx_config_);
}

Op_ptr PauliExpCommutingSetBox::transpose() const {
  std::vector<SymPauliTensor> transpose_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    SymPauliTensor tr = pauli_exp;
    tr.transpose();
    transpose_gadgets.push_back(tr);
  }
  return std::make_shared<PauliExpCommutingSetBox>(
      transpose_gadgets, cx_config_);
}

Op_ptr PauliExpCommutingSetBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  std::vector<SymPauliTensor> symbol_sub_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    symbol_sub_gadgets.push_back(pauli_exp.symbol_substitution(sub_map));
  }
  return std::make_shared<PauliExpCommutingSetBox>(
      symbol_sub_gadgets, this->cx_config_);
}

void PauliExpCommutingSetBox::generate_circuit() const {
  unsigned n_qubits = pauli_gadgets_[0].size();
  Circuit circ(n_qubits);

  std::list<SpSymPauliTensor> gadgets;
  for (const auto &pauli_gadget : pauli_gadgets_) {
    gadgets.push_back((SpSymPauliTensor)pauli_gadget);
  }
  std::set<Qubit> qubits;
  for (unsigned i = 0; i < n_qubits; i++) qubits.insert(Qubit(i));

  Circuit cliff_circ = mutual_diagonalise(gadgets, qubits, cx_config_);

  Circuit phase_poly_circ(n_qubits);

  for (const SpSymPauliTensor &pgp : gadgets) {
    phase_poly_circ.append(pauli_gadget(pgp, CXConfigType::Snake));
  }
  phase_poly_circ.decompose_boxes_recursively();
  PhasePolyBox ppbox(phase_poly_circ);
  Circuit after_synth_circ = *ppbox.to_circuit();

  ConjugationBox box(
      std::make_shared<CircBox>(cliff_circ),
      std::make_shared<CircBox>(after_synth_circ));

  circ.add_box(box, circ.all_qubits());

  circ_ = std::make_shared<Circuit>(circ);
}

bool PauliExpCommutingSetBox::is_equal(const Op &op_other) const {
  const PauliExpCommutingSetBox &other =
      dynamic_cast<const PauliExpCommutingSetBox &>(op_other);
  if (id_ == other.get_id()) return true;
  if (cx_config_ != other.cx_config_) return false;
  return std::equal(
      pauli_gadgets_.begin(), pauli_gadgets_.end(),
      other.pauli_gadgets_.begin(), other.pauli_gadgets_.end(),
      [](const SymPauliTensor &a, const SymPauliTensor &b) {
        return a.equiv_mod(b, 4);
      });
}

nlohmann::json PauliExpCommutingSetBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const PauliExpCommutingSetBox &>(*op);
  nlohmann::json j = core_box_json(box);
  // Encode SymPauliTensor as unlabelled pair of Pauli vector and Expr for
  // backwards compatibility from before templated PauliTensor
  std::vector<std::pair<std::vector<Pauli>, Expr>> gadget_encoding;
  for (const SymPauliTensor &g : box.get_pauli_gadgets())
    gadget_encoding.push_back({g.string, g.coeff});
  j["pauli_gadgets"] = gadget_encoding;
  j["cx_config"] = box.get_cx_config();
  return j;
}

Op_ptr PauliExpCommutingSetBox::from_json(const nlohmann::json &j) {
  std::vector<std::pair<std::vector<Pauli>, Expr>> gadget_encoding =
      j.at("pauli_gadgets")
          .get<std::vector<std::pair<std::vector<Pauli>, Expr>>>();
  std::vector<SymPauliTensor> gadgets;
  for (const std::pair<std::vector<Pauli>, Expr> &g : gadget_encoding)
    gadgets.push_back(SymPauliTensor(g.first, g.second));
  PauliExpCommutingSetBox box =
      PauliExpCommutingSetBox(gadgets, j.at("cx_config").get<CXConfigType>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

REGISTER_OPFACTORY(PauliExpCommutingSetBox, PauliExpCommutingSetBox)

TermSequenceBox::TermSequenceBox(
    const std::vector<SymPauliTensor> &pauli_gadgets,
    Transforms::PauliSynthStrat synth_strategy,
    PauliPartitionStrat partition_strategy, GraphColourMethod graph_colouring,
    CXConfigType cx_configuration, double depth_weight)
    : Box(OpType::TermSequenceBox),
      pauli_gadgets_(pauli_gadgets),
      synth_strategy_(synth_strategy),
      partition_strategy_(partition_strategy),
      graph_colouring_(graph_colouring),
      cx_configuration_(cx_configuration),
      depth_weight_(depth_weight) {
  // check at least one gadget
  if (pauli_gadgets.empty()) {
    signature_ = op_signature_t(0, EdgeType::Quantum);
  } else {
    // check all gadgets have same Pauli string length
    auto n_qubits = pauli_gadgets[0].size();
    for (const auto &gadget : pauli_gadgets) {
      if (gadget.size() != n_qubits) {
        throw PauliExpBoxInvalidity(
            "the Pauli strings within TermSequenceBox must all be the "
            "same length - add Pauli.I to pad strings to required length.");
      }
      if (synth_strategy == Transforms::PauliSynthStrat::Greedy &&
          gadget.size() ==
              std::count(
                  gadget.string.begin(), gadget.string.end(), Pauli::I)) {
        throw PauliExpBoxInvalidity(
            "A TermSequenceBox using the Greedy strategy does not support "
            "rotations with an identity Pauli string. "
            "Note: global phase is not preserved with the Greedy synthesis.");
      }
    }
    signature_ = op_signature_t(n_qubits, EdgeType::Quantum);
  }
}

TermSequenceBox::TermSequenceBox(const TermSequenceBox &other)
    : Box(other),
      pauli_gadgets_(other.pauli_gadgets_),
      synth_strategy_(other.synth_strategy_),
      partition_strategy_(other.partition_strategy_),
      graph_colouring_(other.graph_colouring_),
      cx_configuration_(other.cx_configuration_),
      depth_weight_(other.depth_weight_) {}

TermSequenceBox::TermSequenceBox() : TermSequenceBox({{{}, 0}}) {}

bool TermSequenceBox::is_clifford() const {
  return std::all_of(
      pauli_gadgets_.begin(), pauli_gadgets_.end(),
      [](const SymPauliTensor &pauli_exp) {
        return equiv_0(4 * pauli_exp.coeff) || pauli_exp.string.empty();
      });
}

SymSet TermSequenceBox::free_symbols() const {
  std::vector<Expr> angles;
  for (const auto &pauli_exp : pauli_gadgets_) {
    angles.push_back(pauli_exp.coeff);
  }
  return expr_free_symbols(angles);
}

Op_ptr TermSequenceBox::dagger() const {
  std::vector<SymPauliTensor> dagger_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    dagger_gadgets.emplace_back(pauli_exp.string, -pauli_exp.coeff);
  }
  return std::make_shared<TermSequenceBox>(
      dagger_gadgets, synth_strategy_, partition_strategy_, graph_colouring_,
      cx_configuration_, depth_weight_);
}

Op_ptr TermSequenceBox::transpose() const {
  std::vector<SymPauliTensor> transpose_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    SymPauliTensor tr = pauli_exp;
    tr.transpose();
    transpose_gadgets.push_back(tr);
  }
  return std::make_shared<TermSequenceBox>(
      transpose_gadgets, synth_strategy_, partition_strategy_, graph_colouring_,
      cx_configuration_, depth_weight_);
}

Op_ptr TermSequenceBox::symbol_substitution(
    const SymEngine::map_basic_basic &sub_map) const {
  std::vector<SymPauliTensor> symbol_sub_gadgets;
  for (const auto &pauli_exp : pauli_gadgets_) {
    symbol_sub_gadgets.push_back(pauli_exp.symbol_substitution(sub_map));
  }
  return std::make_shared<TermSequenceBox>(
      symbol_sub_gadgets, synth_strategy_, partition_strategy_,
      graph_colouring_, cx_configuration_, depth_weight_);
}

bool TermSequenceBox::is_equal(const Op &op_other) const {
  const TermSequenceBox &other =
      dynamic_cast<const TermSequenceBox &>(op_other);
  if (id_ == other.get_id()) return true;
  if (synth_strategy_ != other.synth_strategy_) return false;
  if (partition_strategy_ != other.partition_strategy_) return false;
  if (graph_colouring_ != other.graph_colouring_) return false;
  if (cx_configuration_ != other.cx_configuration_) return false;
  if (depth_weight_ != other.depth_weight_) return false;
  return std::equal(
      pauli_gadgets_.begin(), pauli_gadgets_.end(),
      other.pauli_gadgets_.begin(), other.pauli_gadgets_.end(),
      [](const SymPauliTensor &a, const SymPauliTensor &b) {
        return a.equiv_mod(b, 4);
      });
}

Transforms::PauliSynthStrat TermSequenceBox::get_synth_strategy() const {
  return synth_strategy_;
}

PauliPartitionStrat TermSequenceBox::get_partition_strategy() const {
  return partition_strategy_;
}

GraphColourMethod TermSequenceBox::get_graph_colouring() const {
  return graph_colouring_;
}

std::vector<SymPauliTensor> TermSequenceBox::get_pauli_gadgets() const {
  return pauli_gadgets_;
}

CXConfigType TermSequenceBox::get_cx_config() const {
  return cx_configuration_;
}

double TermSequenceBox::get_depth_weight() const { return depth_weight_; }

nlohmann::json TermSequenceBox::to_json(const Op_ptr &op) {
  const auto &box = static_cast<const TermSequenceBox &>(*op);
  nlohmann::json j = core_box_json(box);
  // Encode SymPauliTensor as unlabelled pair of Pauli vector and Expr for
  // backwards compatibility from before templated PauliTensor
  std::vector<std::pair<std::vector<Pauli>, Expr>> gadget_encoding;
  for (const SymPauliTensor &g : box.get_pauli_gadgets())
    gadget_encoding.push_back({g.string, g.coeff});
  j["pauli_gadgets"] = gadget_encoding;
  j["synth_strategy"] = box.get_synth_strategy();
  j["partition_strategy"] = box.get_partition_strategy();
  j["graph_colouring"] = box.get_graph_colouring();
  j["cx_config"] = box.get_cx_config();
  j["depth_weight"] = box.get_depth_weight();

  return j;
}

Op_ptr TermSequenceBox::from_json(const nlohmann::json &j) {
  std::vector<std::pair<std::vector<Pauli>, Expr>> gadget_encoding =
      j.at("pauli_gadgets")
          .get<std::vector<std::pair<std::vector<Pauli>, Expr>>>();
  std::vector<SymPauliTensor> gadgets;
  for (const std::pair<std::vector<Pauli>, Expr> &g : gadget_encoding)
    gadgets.push_back(SymPauliTensor(g.first, g.second));
  TermSequenceBox box = TermSequenceBox(
      gadgets, j.at("synth_strategy").get<Transforms::PauliSynthStrat>(),
      j.at("partition_strategy").get<PauliPartitionStrat>(),
      j.at("graph_colouring").get<GraphColourMethod>(),
      j.at("cx_config").get<CXConfigType>(),
      j.at("depth_weight").get<double>());
  return set_box_id(
      box,
      boost::lexical_cast<boost::uuids::uuid>(j.at("id").get<std::string>()));
}

DensePauliMap pad_sparse_pauli_map(
    const QubitPauliMap &sparse_string, size_t size) {
  // n.b. as a helper method this can assume that all qubits in sparse string
  // are labelled with the default qubit register and no index is larger than or
  // equal to the size
  std::vector<Pauli> dense_string(size, Pauli::I);
  for (const auto &letter : sparse_string) {
    TKET_ASSERT(letter.first.reg_name() == q_default_reg());
    std::vector<unsigned> index = letter.first.index();
    TKET_ASSERT(index.size() == 1);
    TKET_ASSERT(index[0] < size);
    dense_string[index[0]] = letter.second;
  }
  return dense_string;
}

void TermSequenceBox::generate_circuit() const {
  // constructor guarantees all gadgets are the same size
  unsigned n_qubits = pauli_gadgets_[0].size();
  Circuit circ(n_qubits);

  // Saves later handling
  if (this->pauli_gadgets_.empty()) {
    circ_ = std::make_shared<Circuit>(circ);
    return;
  }

  // First combine any PauliGadgets with the same dense string into the same
  // term
  // n.b. the map assigns an ordering to the Pauli terms, meaning it
  // is an algorithm-affecting step
  std::map<DensePauliMap, Expr> reduced_pauli_gadgets;
  for (const SymPauliTensor &pauli_gadget : this->pauli_gadgets_) {
    auto it = reduced_pauli_gadgets.find(pauli_gadget.string);
    if (it != reduced_pauli_gadgets.end()) {
      it->second += pauli_gadget.coeff;
    } else {
      reduced_pauli_gadgets.insert(
          std::make_pair(pauli_gadget.string, pauli_gadget.coeff));
    }
  }

  // Construct circuit depending on specified strategy, using other Pauli
  // Exponential Boxes
  switch (this->synth_strategy_) {
    case Transforms::PauliSynthStrat::Individual: {
      for (auto it = reduced_pauli_gadgets.begin();
           it != reduced_pauli_gadgets.end(); ++it) {
        circ.add_box(
            PauliExpBox(
                SymPauliTensor(it->first, it->second), this->cx_configuration_),
            circ.all_qubits());
      }
      break;
    }
    case Transforms::PauliSynthStrat::Pairwise: {
      auto it = reduced_pauli_gadgets.begin();
      // if there is an odd number of gadgets, add first gadget as a single
      if (reduced_pauli_gadgets.size() % 2 == 1) {
        circ.add_box(
            PauliExpBox(
                SymPauliTensor(it->first, it->second), this->cx_configuration_),
            circ.all_qubits());
        ++it;
      }
      while (it != reduced_pauli_gadgets.end()) {
        auto next_it = std::next(it);
        TKET_ASSERT(next_it != reduced_pauli_gadgets.end());
        circ.add_box(
            PauliExpPairBox(
                SymPauliTensor(it->first, it->second),
                SymPauliTensor(next_it->first, next_it->second),
                this->cx_configuration_),
            circ.all_qubits());
        std::advance(it, 2);
      }
      break;
    }
    case Transforms::PauliSynthStrat::Sets: {
      // term_sequence expects std::list<SpPauliString>, i.e. QubitPauliMap
      // instead of DensePauliMap first convert the keys of
      // reduced_pauli_gadgets to a std::list<SpPauliString> object
      std::list<SpPauliString> all_terms;
      for (const auto &gadget : reduced_pauli_gadgets) {
        all_terms.push_back(SpPauliString(gadget.first));
      }
      // then get the term sequence into commuting sets
      std::list<std::list<SpPauliString>> commuting_sets = term_sequence(
          all_terms, this->partition_strategy_, this->graph_colouring_);
      // construct each commuting set using PauliExpCommutingSetBox:
      // 1) pad each QubitPauliMap into a DensePauliMap by adding Pauli.I terms
      // for missing keys 2) retrieve the exponential coefficient from
      // reduced_pauli_gadgets 3) construct a SymPauliTensor object as the
      // argument to PauliExpCommutingSetBox
      for (const std::list<SpPauliString> &strings : commuting_sets) {
        std::vector<SymPauliTensor> commuting_gadgets;
        for (const SpPauliString &string : strings) {
          DensePauliMap padded_string =
              pad_sparse_pauli_map(string.string, n_qubits);
          commuting_gadgets.push_back(SymPauliTensor(
              padded_string, reduced_pauli_gadgets[padded_string]));
        }
        circ.add_box(
            PauliExpCommutingSetBox(commuting_gadgets, this->cx_configuration_),
            circ.all_qubits());
      }
      break;
    }
    case Transforms::PauliSynthStrat::Greedy: {
      std::vector<SymPauliTensor> unordered_gadgets;
      for (const auto &pair : reduced_pauli_gadgets) {
        unordered_gadgets.push_back(SymPauliTensor(pair.first, pair.second));
      }
      circ = Transforms::GreedyPauliSimp::greedy_pauli_set_synthesis(
          unordered_gadgets, this->depth_weight_);
      Transforms::singleq_clifford_sweep().apply(circ);
      break;
    }
    default: {
      throw std::logic_error(
          "TermSequenceBox passed an unsupported PauliSynthStrat");
    }
  }
  circ_ = std::make_shared<Circuit>(circ);
}

REGISTER_OPFACTORY(TermSequenceBox, TermSequenceBox)

void append_single_pauli_gadget_as_pauli_exp_box(
    Circuit &circ, const SpSymPauliTensor &pauli, CXConfigType cx_config) {
  std::vector<Pauli> string;
  std::vector<Qubit> mapping;
  for (const std::pair<const Qubit, Pauli> &term : pauli.string) {
    string.push_back(term.second);
    mapping.push_back(term.first);
  }
  PauliExpBox box(SymPauliTensor(string, pauli.coeff), cx_config);
  circ.add_box(box, mapping);
}

void append_pauli_gadget_pair_as_box(
    Circuit &circ, const SpSymPauliTensor &pauli0,
    const SpSymPauliTensor &pauli1, CXConfigType cx_config) {
  std::vector<Qubit> mapping;
  std::vector<Pauli> paulis0;
  std::vector<Pauli> paulis1;
  QubitPauliMap p1map = pauli1.string;
  // add paulis for qubits in pauli0_string
  for (const std::pair<const Qubit, Pauli> &term : pauli0.string) {
    mapping.push_back(term.first);
    paulis0.push_back(term.second);
    auto found = p1map.find(term.first);
    if (found == p1map.end()) {
      paulis1.push_back(Pauli::I);
    } else {
      paulis1.push_back(found->second);
      p1map.erase(found);
    }
  }
  // add paulis for qubits in pauli1_string that weren't in pauli0_string
  for (const std::pair<const Qubit, Pauli> &term : p1map) {
    mapping.push_back(term.first);
    paulis1.push_back(term.second);
    paulis0.push_back(Pauli::I);  // If pauli0_string contained qubit, would
                                  // have been handled above
  }
  PauliExpPairBox box(
      SymPauliTensor(paulis0, pauli0.coeff),
      SymPauliTensor(paulis1, pauli1.coeff), cx_config);
  circ.add_box(box, mapping);
}

void append_commuting_pauli_gadget_set_as_box(
    Circuit &circ, const std::list<SpSymPauliTensor> &gadgets,
    CXConfigType cx_config) {
  // Translate from QubitPauliTensors to vectors of Paulis of same length
  // Preserves ordering of qubits

  std::set<Qubit> all_qubits;
  for (const SpSymPauliTensor &gadget : gadgets) {
    for (const std::pair<const Qubit, Pauli> &qubit_pauli : gadget.string) {
      all_qubits.insert(qubit_pauli.first);
    }
  }

  std::vector<Qubit> mapping;
  for (const auto &qubit : all_qubits) {
    mapping.push_back(qubit);
  }

  std::vector<SymPauliTensor> pauli_gadgets;
  for (const SpSymPauliTensor &gadget : gadgets) {
    SymPauliTensor &new_gadget =
        pauli_gadgets.emplace_back(DensePauliMap{}, gadget.coeff);
    for (const Qubit &qubit : mapping) {
      new_gadget.string.push_back(gadget.get(qubit));
    }
  }

  PauliExpCommutingSetBox box(pauli_gadgets, cx_config);
  circ.add_box(box, mapping);
}

}  // namespace tket
