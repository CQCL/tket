#include "tket/Clifford/ChoiMixTableau.hpp"
#include "tket/Diagonalisation/Diagonalisation.hpp"

namespace tket {

ChoiMixTableau tab_from_gadgets(const std::list<SpSymPauliTensor> &gadgets) {
  std::list<ChoiMixTableau::row_tensor_t> rows;
  for (const SpSymPauliTensor &g : gadgets) {
    rows.push_back({{}, SpPauliStabiliser(g.string)});
  }
  // TODO:: May need to check for independence of rows
  return ChoiMixTableau(rows);
}

Circuit mutual_diagonalise_JGM(
    std::list<SpSymPauliTensor> &gadgets, std::set<Qubit> qubits,
    CXConfigType) {
  // See Appendix A of Jena, Genin, Mosca, "Pauli Partitioning with Respect to
  // Gate Sets", https://arxiv.org/pdf/1907.07859.pdf
  Circuit cliff_circ;
  for (const Qubit &qb : qubits) {
    cliff_circ.add_qubit(qb);
  }
  ChoiMixTableau tab = tab_from_gadgets(gadgets);
  const MatrixXb &xmat = tab.tab_.xmat_;
  const MatrixXb &zmat = tab.tab_.zmat_;

  // Translation from paper
  // F_2 = H
  // R_2 = S
  // SUM_2 = CX
  for (unsigned q = 0; q < xmat.cols(); ++q) {
    for (unsigned r = 0; r < xmat.rows(); ++r) {
      if (xmat(r, q)) {
        // Found a nonzero X component
        // H1: we skip the SUM operator moving this to the first qubit, instead
        // change which qubits gates act on

        // H2: eliminate the X component at that qubit on that row
        Qubit q_qb = tab.col_index_.right.at(q).first;
        if (zmat(r, q)) {
          cliff_circ.add_op<Qubit>(OpType::S, {q_qb});
          tab.apply_gate(OpType::S, {q_qb});
        }
        cliff_circ.add_op<Qubit>(OpType::H, {q_qb});
        tab.apply_gate(OpType::H, {q_qb});

        // H3: remove other Z components on row r
        for (unsigned q2 = 0; q2 < xmat.cols(); ++q2) {
          if (q2 != q && zmat(r, q2)) {
            Qubit q2_qb = tab.col_index_.right.at(q2).first;
            cliff_circ.add_op<Qubit>(OpType::CX, {q2_qb, q_qb});
            tab.apply_gate(OpType::CX, {q2_qb, q_qb});
          }
        }

        // H4/5: remove other X components on row r
        for (unsigned q2 = 0; q2 < xmat.cols(); ++q2) {
          if (q2 != q && xmat(r, q2)) {
            Qubit q2_qb = tab.col_index_.right.at(q2).first;
            cliff_circ.add_op<Qubit>(OpType::H, {q2_qb});
            tab.apply_gate(OpType::H, {q2_qb});
            cliff_circ.add_op<Qubit>(OpType::CX, {q2_qb, q_qb});
            tab.apply_gate(OpType::CX, {q2_qb, q_qb});
          }
        }
      }
    }
  }

  return cliff_circ;
}

}  // namespace tket
