// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "Circuit/Circuit.hpp"

namespace tket {

struct FlowVertProperties {
  Circuit circ;
  std::optional<Bit> branch_condition;
  std::optional<std::string> label;
};

struct FlowEdgeProperties {
  bool branch;
};

typedef boost::adjacency_list<
    boost::listS, boost::listS, boost::bidirectionalS, FlowVertProperties,
    FlowEdgeProperties>
    FlowGraph;
typedef boost::graph_traits<FlowGraph>::vertex_descriptor FGVert;
typedef boost::graph_traits<FlowGraph>::edge_descriptor FGEdge;

typedef std::vector<FGVert> FGVertVec;
typedef std::vector<FGEdge> FGEdgeVec;

typedef boost::multi_index::multi_index_container<
    UnitID, boost::multi_index::indexed_by<
                boost::multi_index::ordered_unique<
                    boost::multi_index::tag<TagID>,
                    boost::multi_index::identity<UnitID> >,
                boost::multi_index::ordered_non_unique<
                    boost::multi_index::tag<TagType>,
                    boost::multi_index::const_mem_fun<
                        UnitID, UnitType, &UnitID::type> >,
                boost::multi_index::ordered_non_unique<
                    boost::multi_index::tag<TagReg>,
                    boost::multi_index::const_mem_fun<
                        UnitID, std::string, &UnitID::reg_name> > > >
    unit_lookup_t;

class ProgramError : public std::logic_error {
 public:
  explicit ProgramError(const std::string &message)
      : std::logic_error(message) {}
};

/**
 * A control flow graph where basic blocks are Circuits.
 *
 * Each block ends with either an unconditional jump, or
 * branching based on a single bit condition.
 *
 * Each block can contain OpenQASM-style conditional gates
 * using the Conditional without the need for explicit branching.
 */
class Program {
 public:
  Program();
  explicit Program(unsigned qubits, unsigned bits = 0);
  Program(const Program &to_copy);

  /** Unit control */
  qubit_vector_t all_qubits() const;
  bit_vector_t all_bits() const;
  std::vector<UnitID> all_units() const;

  std::map<Bit, unsigned> bit_readout() const;
  std::map<Qubit, unsigned> qubit_readout() const;

  opt_reg_info_t get_reg_info(std::string reg_name) const;
  register_t get_reg(std::string reg_name) const;

  void add_qubit(const Qubit &id, bool reject_dups = true);
  void add_bit(const Bit &id, bool reject_dups = true);
  register_t add_q_register(std::string reg_name, unsigned size);
  register_t add_c_register(std::string reg_name, unsigned size);

  /** Graph accessors */
  Circuit &get_circuit_ref(const FGVert &vert);
  const Circuit &get_circuit_ref(const FGVert &vert) const;
  std::optional<Bit> get_condition(const FGVert &vert) const;
  std::optional<std::string> get_label(const FGVert &vert) const;

  bool get_branch(const FGEdge &edge) const;

  FGVertVec get_successors(const FGVert &vert) const;
  FGVertVec get_predecessors(const FGVert &vert) const;

  FGVert get_branch_successor(const FGVert &vert, bool branch = false) const;

  FGEdgeVec get_in_edges(const FGVert &vert) const;
  FGEdgeVec get_out_edges(const FGVert &vert) const;
  unsigned n_in_edges(const FGVert &vert) const;
  unsigned n_out_edges(const FGVert &vert) const;

  FGVert get_source(const FGEdge &edge) const;
  FGVert get_target(const FGEdge &edge) const;

  unsigned get_n_vertices() const;

  /** Graph manipulation */
  FGVert add_vertex(
      const Circuit &circ, std::optional<Bit> branch_condition = std::nullopt,
      const std::optional<std::string> &label = std::nullopt);
  void remove_vertex(const FGVert &vert);

  FGEdge add_edge(
      const FGVert &source, const FGVert &target, bool branch = false);
  void remove_edge(const FGEdge &edge);

  std::map<FGVert, FGVert> copy_graph(const Program &to_copy);

  /** Adding instructions to program */
  FGVert add_block(const Circuit &circ);

  void add_op(const Op_ptr &op, const std::vector<unsigned> &args);
  void add_op(const Op_ptr &op, const unit_vector_t &args);

  template <class ID>
  void add_op(OpType type, const std::vector<ID> &args) {
    add_op(get_op_ptr(type, std::vector<Expr>{}, args.size()), args);
  }
  template <class ID>
  void add_op(OpType type, const Expr &param, const std::vector<ID> &args) {
    add_op(get_op_ptr(type, param, args.size()), args);
  }
  template <class ID>
  void add_op(
      OpType type, const std::vector<Expr> &params,
      const std::vector<ID> &args) {
    add_op(get_op_ptr(type, params, args.size()), args);
  }

  void append(const Program &to_append);
  void append_if(const Bit &condition_bit, const Program &body);
  void append_if_else(
      const Bit &condition_bit, const Program &if_body,
      const Program &else_body);
  void append_while(const Bit &condition_bit, const Program &body);

  friend Program operator>>(const Program &p1, const Program &p2);

  /** Graph analysis */

  /**
   * Checks whether flow graph is in the correct format.
   * Returning false should be seen as a fatal error.
   */
  bool check_valid() const;

  void to_graphviz_file(const std::string &filename) const;
  void to_graphviz(std::ostream &out) const;

  // Insert Control Flow optimisations here

  /** Command Iteration */

  /**
   * Iterates through the vertices of the flow graph in a depth-first,
   * preorder traversal
   */
  class BlockIterator {
   public:
    BlockIterator();
    explicit BlockIterator(const Program &prog);

    const FGVert &operator*() const { return current_vert_; }
    const FGVert *operator->() const { return &current_vert_; }
    bool operator==(const BlockIterator &other) const {
      return this->current_vert_ == other.current_vert_;
    }
    bool operator!=(const BlockIterator &other) const {
      return !(*this == other);
    }

    BlockIterator operator++(int);
    BlockIterator &operator++();

    const Circuit &get_circuit_ref() const {
      return prog_->get_circuit_ref(current_vert_);
    }
    std::optional<Bit> get_condition() const {
      return prog_->get_condition(current_vert_);
    }

   private:
    const Program *prog_;
    FGVert current_vert_;
    std::list<FGVert> stack_;
    std::set<FGVert> visited_;
  };

  BlockIterator block_begin() const;
  static BlockIterator block_end();

  /**
   * Iterates through the commands of the program.
   * Commands in each block are given by the iteration order of the Circuit
   * class. Blocks are then ordered according to BlockIterator. Labels, jumps,
   * branches and exit commands are inserted at beginning and end of blocks.
   */
  class CommandIterator {
   public:
    explicit CommandIterator(const Program &prog);
    CommandIterator() : stage_(ComItStage::End) {}

    Command operator*() const { return current_command_; }
    const Command *operator->() const { return &current_command_; }
    bool operator==(const CommandIterator &other) const {
      return (current_block_ == other.current_block_) &&
             (current_com_iterator_ == other.current_com_iterator_) &&
             (stage_ == other.stage_);
    }
    bool operator!=(const CommandIterator &other) const {
      return !(*this == other);
    }
    // Postfix incrementer
    CommandIterator operator++(int);
    // Prefix incrementer
    CommandIterator &operator++();

    friend std::ostream &operator<<(
        std::ostream &out, const CommandIterator &cit) {
      out << cit->to_str();
      return out;
    };

   private:
    Command current_command_;
    BlockIterator current_block_;
    Circuit::CommandIterator current_com_iterator_;
    std::map<FGVert, std::string> labels_;
    const Program *prog_;

    /**
     * Stage of command iteration to come next
     */
    enum class ComItStage {
      Label,         // Just reached a new block; possibly add a label
      FirstCommand,  // Just output a label; start iterating through commands of
                     // block
      Command,       // Just output a command; continue
      Branch,        // Finished commands for block; possibly add a branch
      Goto,   // Reached end of block; possibly add a Goto before moving to next
              // block
      Stop,   // Just output label for exit block; output STOP
      Final,  // Just output STOP; move to end()
      End     // Reserved for end()
    };
    ComItStage stage_;   // The next stage of a block's commands to try
    FGVert prev_block_;  // Used to determine if a label is needed

    std::string get_label(const FGVert &block);
  };

  CommandIterator begin() const;
  static CommandIterator end();

 private:
  FlowGraph flow_;
  FGVert entry_;
  FGVert exit_;

  unit_lookup_t units_;
};

}  // namespace tket
