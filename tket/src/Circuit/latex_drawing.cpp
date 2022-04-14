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

#include <limits>

#include "Boxes.hpp"
#include "Circuit.hpp"
#include "OpType/OpType.hpp"
namespace tket {

struct LineBufferInfo {
  std::stringstream buffer;
  unsigned depth;
  bool is_quantum;
};

struct LatexContext {
  std::map<UnitID, unsigned> line_ids;
  std::vector<LineBufferInfo> lines;
};

void add_latex_for_command(LatexContext& context, const Command& command) {
  std::map<UnitID, unsigned>& line_ids = context.line_ids;
  std::vector<LineBufferInfo>& lines = context.lines;
  unit_vector_t args = command.get_args();
  const Op_ptr op = command.get_op_ptr();
  switch (op->get_type()) {
    case OpType::CnRy:
    case OpType::CnX: {
      int target_index = line_ids.at(args.back());
      args.pop_back();
      for (const UnitID& u : args) {
        int control_index = line_ids.at(u);
        lines.at(control_index).buffer
            << "\\ctrl{" << target_index - control_index << "} & ";
        lines.at(control_index).depth++;
      }
      if (op->get_type() == OpType::CnRy) {
        lines.at(target_index).buffer
            << "\\gate{\\text{"
            << get_op_ptr(OpType::Ry, op->get_params())->get_name(true)
            << "}} & ";
      } else {
        lines.at(target_index).buffer << "\\targ{} & ";
      }
      lines.at(target_index).depth++;
      break;
    }
    case OpType::CCX: {
      int control0_index = line_ids.at(args.at(0));
      int control1_index = line_ids.at(args.at(1));
      int target_index = line_ids.at(args.at(2));
      lines.at(control0_index).buffer
          << "\\ctrl{" << target_index - control0_index << "} & ";
      lines.at(control0_index).depth++;
      lines.at(control1_index).buffer
          << "\\ctrl{" << target_index - control1_index << "} & ";
      lines.at(control1_index).depth++;
      lines.at(target_index).buffer << "\\targ{} & ";
      lines.at(target_index).depth++;
      break;
    }
    case OpType::CSWAP: {
      int control_index = line_ids.at(args.at(0));
      int target0_index = line_ids.at(args.at(1));
      int target1_index = line_ids.at(args.at(2));
      lines.at(control_index).buffer << "\\ctrl{"
                                     << target0_index - control_index << "} & ";
      lines.at(control_index).depth++;
      lines.at(target0_index).buffer << "\\swap{"
                                     << target1_index - target0_index << "} & ";
      lines.at(target0_index).depth++;
      lines.at(target1_index).buffer << "\\targX{} & ";
      lines.at(target1_index).depth++;
      break;
    }
    case OpType::CH:
    case OpType::CRz:
    case OpType::CRx:
    case OpType::CRy:
    case OpType::CU1:
    case OpType::CU3:
    case OpType::CY: {
      std::stringstream gate_name;
      switch (op->get_type()) {
        case OpType::CH: {
          gate_name << get_op_ptr(OpType::H)->get_name(true);
          break;
        }
        case OpType::CRz: {
          gate_name << get_op_ptr(OpType::Rz, op->get_params())->get_name(true);
          break;
        }
        case OpType::CRx: {
          gate_name << get_op_ptr(OpType::Rx, op->get_params())->get_name(true);
          break;
        }
        case OpType::CRy: {
          gate_name << get_op_ptr(OpType::Ry, op->get_params())->get_name(true);
          break;
        }
        case OpType::CU1: {
          gate_name << get_op_ptr(OpType::U1, op->get_params())->get_name(true);
          break;
        }
        case OpType::CU3: {
          gate_name << get_op_ptr(OpType::U3, op->get_params())->get_name(true);
          break;
        }
        case OpType::CY: {
          gate_name << get_op_ptr(OpType::Y)->get_name(true);
          break;
        }
        case OpType::CV: {
          gate_name << get_op_ptr(OpType::V)->get_name(true);
          break;
        }
        case OpType::CVdg: {
          gate_name << get_op_ptr(OpType::Vdg)->get_name(true);
          break;
        }
        case OpType::CSX: {
          gate_name << get_op_ptr(OpType::SX)->get_name(true);
          break;
        }
        case OpType::CSXdg: {
          gate_name << get_op_ptr(OpType::SXdg)->get_name(true);
          break;
        }
        default: {
        }
      }
      int control_index = line_ids.at(args.at(0));
      int target_index = line_ids.at(args.at(1));
      lines.at(control_index).buffer << "\\ctrl{"
                                     << target_index - control_index << "} & ";
      lines.at(control_index).depth++;
      lines.at(target_index).buffer << "\\gate{\\text{" << gate_name.str()
                                    << "}} & ";
      lines.at(target_index).depth++;
      break;
    }
    case OpType::CX: {
      int control_index = line_ids.at(args.at(0));
      int target_index = line_ids.at(args.at(1));
      lines.at(control_index).buffer << "\\ctrl{"
                                     << target_index - control_index << "} & ";
      lines.at(control_index).depth++;
      lines.at(target_index).buffer << "\\targ{} & ";
      lines.at(target_index).depth++;
      break;
    }
    case OpType::CZ: {
      int control_index = line_ids.at(args.at(0));
      int target_index = line_ids.at(args.at(1));
      lines.at(control_index).buffer << "\\ctrl{"
                                     << target_index - control_index << "} & ";
      lines.at(control_index).depth++;
      lines.at(target_index).buffer << "\\control{} & ";
      lines.at(target_index).depth++;
      break;
    }
    case OpType::Measure: {
      int qb_index = line_ids.at(args.at(0));
      int cb_index = line_ids.at(args.at(1));
      lines.at(qb_index).buffer << "\\meter{} \\vcw{" << cb_index - qb_index
                                << "} & ";
      lines.at(qb_index).depth++;
      lines.at(cb_index).buffer << "\\cw & ";
      lines.at(cb_index).depth++;
      break;
    }
    case OpType::Collapse: {
      int qb_index = line_ids.at(args.at(0));
      lines.at(qb_index).buffer << "\\meter{} & ";
      lines.at(qb_index).depth++;
      break;
    }
    case OpType::SWAP: {
      int qb0_index = line_ids.at(args.at(0));
      int qb1_index = line_ids.at(args.at(1));
      lines.at(qb0_index).buffer << "\\swap{" << qb1_index - qb0_index
                                 << "} & ";
      lines.at(qb0_index).depth++;
      lines.at(qb1_index).buffer << "\\targX{} & ";
      lines.at(qb1_index).depth++;
      break;
    }
    case OpType::QControlBox: {
      const QControlBox& box = static_cast<const QControlBox&>(*op);
      unsigned n_controls = box.get_n_controls();
      int max_index = 0;
      unit_vector_t target_args;
      for (unsigned i = n_controls; i < args.size(); ++i) {
        target_args.push_back(args[i]);
        int index = line_ids.at(args[i]);
        if (index > max_index) max_index = index;
      }
      add_latex_for_command(context, Command(box.get_op(), target_args));
      for (unsigned i = 0; i < n_controls; ++i) {
        int index = line_ids.at(args[i]);
        lines.at(index).buffer << "\\ctrl{" << max_index - index << "} & ";
        lines.at(index).depth++;
      }
      break;
    }
    case OpType::Conditional: {
      const Conditional& box = static_cast<const Conditional&>(*op);
      unsigned n_controls = box.get_width();
      int max_index = 0;
      unit_vector_t target_args;
      for (unsigned i = n_controls; i < args.size(); ++i) {
        target_args.push_back(args[i]);
        int index = line_ids.at(args[i]);
        if (index > max_index) max_index = index;
      }
      add_latex_for_command(context, Command(box.get_op(), target_args));
      for (unsigned i = 0; i < n_controls; ++i) {
        int index = line_ids.at(args[i]);
        lines.at(index).buffer << "\\cwbend{" << max_index - index << "} & ";
        lines.at(index).depth++;
      }
      break;
    }
    default: {
      unsigned min_index = std::numeric_limits<unsigned>::max();
      unsigned max_index = 0;
      for (const UnitID& arg : args) {
        unsigned index = line_ids.at(arg);
        if (index < min_index) min_index = index;
        if (index > max_index) max_index = index;
      }
      lines.at(min_index).buffer << "\\gate[" << (max_index + 1 - min_index)
                                 << "]{\\text{" + op->get_name(true) + "}} & ";
      lines.at(min_index).depth++;
      for (const UnitID& arg : args) {
        unsigned index = line_ids.at(arg);
        if (index != min_index) {
          lines.at(index).buffer << "& ";
          lines.at(index).depth++;
        }
      }
    }
  }
}

std::string Circuit::to_latex_str() const {
  // Initial header
  std::stringstream buffer;
  buffer << "\\documentclass[tikz]{standalone}\n";
  buffer << "\\usetikzlibrary{quantikz}\n";
  buffer << "\\begin{document}\n";
  buffer << "\\begin{quantikz}\n";

  // Wire labels
  LatexContext context;
  std::map<UnitID, unsigned>& line_ids = context.line_ids;
  std::vector<LineBufferInfo>& lines = context.lines;
  for (const Qubit& qb : this->all_qubits()) {
    unsigned n_lines = lines.size();
    line_ids.insert({qb, n_lines});
    LineBufferInfo& line = *lines.emplace(lines.end());
    line.buffer << "\\lstick{" + qb.repr() + "} & ";
    line.is_quantum = true;
  }
  for (const Bit& cb : this->all_bits()) {
    unsigned n_lines = lines.size();
    line_ids.insert({cb, n_lines});
    LineBufferInfo& line = *lines.emplace(lines.end());
    line.buffer << "\\lstick{" + cb.repr() + "} & ";
    line.is_quantum = false;
  }

  // Commands
  for (const Command& com : this->get_commands()) {
    std::set<unsigned> used_lines;
    unsigned min_index = std::numeric_limits<unsigned>::max();
    unsigned max_index = 0;
    for (const UnitID& arg : com.get_args()) {
      used_lines.insert(line_ids.at(arg));
    }
    for (unsigned index : used_lines) {
      if (index < min_index) min_index = index;
      if (index > max_index) max_index = index;
    }

    unsigned max_depth = 0;
    for (unsigned index = min_index; index <= max_index; index++) {
      if (lines.at(index).depth > max_depth) max_depth = lines.at(index).depth;
    }
    for (unsigned index : used_lines) {
      for (unsigned d = lines.at(index).depth; d < max_depth; d++) {
        if (lines.at(index).is_quantum) {
          lines.at(index).buffer << "\\qw & ";
        } else {
          lines.at(index).buffer << "\\cw & ";
        }
      }
      lines.at(index).depth = max_depth;
    }

    add_latex_for_command(context, com);

    for (unsigned index = min_index; index <= max_index; index++) {
      for (unsigned d = lines.at(index).depth; d <= max_depth; d++) {
        if (lines.at(index).is_quantum) {
          lines.at(index).buffer << "\\qw & ";
        } else {
          lines.at(index).buffer << "\\cw & ";
        }
      }
      lines.at(index).depth = max_depth + 1;
    }
  }

  // Fill out ends
  unsigned max_depth = 0;
  for (const LineBufferInfo& l : lines) {
    if (l.depth > max_depth) max_depth = l.depth;
  }
  for (LineBufferInfo& l : lines) {
    for (unsigned d = l.depth; d < max_depth; d++) {
      if (l.is_quantum) {
        l.buffer << "\\qw & ";
      } else {
        l.buffer << "\\cw & ";
      }
    }
  }
  for (LineBufferInfo& l : lines) {
    if (l.is_quantum) {
      l.buffer << "\\qw \\\\";
    } else {
      l.buffer << "\\cw \\\\";
    }
  }

  // Build up overall buffer
  for (LineBufferInfo& l : lines) {
    buffer << l.buffer.str() << "\n";
  }
  buffer << "\\end{quantikz}\n";
  buffer << "\\end{document}";

  return buffer.str();
}

void Circuit::to_latex_file(const std::string& filename) const {
  std::ofstream file(filename);
  file << to_latex_str();
}

}  // namespace tket
