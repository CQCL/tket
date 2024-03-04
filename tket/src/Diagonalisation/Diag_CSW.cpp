#include "tket/Clifford/ChoiMixTableau.hpp"
#include "tket/Diagonalisation/Diagonalisation.hpp"

namespace tket {

Circuit mutual_diagonalise_CSW_CZ(
    std::list<SpSymPauliTensor> &, std::set<Qubit> qubits,
    CXConfigType) {
  Circuit cliff_circ;
  for (const Qubit &qb : qubits) {
    cliff_circ.add_qubit(qb);
  }
  // ChoiMixTableau tab = tab_from_gadgets(gadgets);
  // const MatrixXb &xmat = tab.tab_.xmat_;
  // const MatrixXb &zmat = tab.tab_.zmat_;

  // // Make xmat have full rank
  // MatrixXb echelon = xmat;

  return cliff_circ;
}

Circuit mutual_diagonalise_CSW_CX(
    std::list<SpSymPauliTensor> &, std::set<Qubit> qubits,
    CXConfigType) {
  Circuit cliff_circ;
  for (const Qubit &qb : qubits) {
    cliff_circ.add_qubit(qb);
  }
  // ChoiMixTableau tab = tab_from_gadgets(gadgets);
  // const MatrixXb &xmat = tab.tab_.xmat_;
  // const MatrixXb &zmat = tab.tab_.zmat_;

  return cliff_circ;
}

}  // namespace tket
