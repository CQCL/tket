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

#include "Ops/FlowOp.hpp"
#include "Program.hpp"

namespace tket {

Program::BlockIterator::BlockIterator()
    : prog_(nullptr),
      current_vert_(boost::graph_traits<FlowGraph>::null_vertex()) {}

Program::BlockIterator::BlockIterator(const Program &p) {
  FGVert first = p.get_successors(p.entry_).front();
  if (first == p.exit_) {
    current_vert_ = boost::graph_traits<FlowGraph>::null_vertex();
    return;
  }
  prog_ = &p;
  current_vert_ = first;
  stack_.push_back(first);
  visited_.insert(first);
}

Program::BlockIterator Program::BlockIterator::operator++(int) {
  BlockIterator it = *this;
  ++*this;
  return it;
}

Program::BlockIterator &Program::BlockIterator::operator++() {
  while (!stack_.empty()) {
    FGVert last_branch = stack_.back();
    FGVertVec succs = prog_->get_successors(last_branch);
    for (const FGVert &s : succs) {
      if (visited_.find(s) == visited_.end() && s != prog_->exit_) {
        current_vert_ = s;
        stack_.push_back(s);
        visited_.insert(s);
        return *this;
      }
    }
    stack_.pop_back();
  }
  *this = BlockIterator();
  return *this;
}

Program::BlockIterator Program::block_begin() const {
  return BlockIterator(*this);
}

Program::BlockIterator Program::block_end() { return BlockIterator(); }

Program::CommandIterator::CommandIterator(const Program &prog) {
  prog_ = &prog;
  current_block_ = prog_->block_begin();
  stage_ = ComItStage::Label;
  prev_block_ = prog_->entry_;
  ++*this;
}

Program::CommandIterator Program::CommandIterator::operator++(int) {
  CommandIterator it = *this;
  ++*this;
  return it;
}

Program::CommandIterator &Program::CommandIterator::operator++() {
  /**
   * Approximately follows routine for one-pass code generation from Aho, Lam,
   * Sethi, Ullman, section 6.7.1 (No need for backpatching when labels with
   * arbitrary names allowed)
   *
   * for (Block b in prog)
   *      if b.preds.size > 1 yield LABEL
   *      for (Command c in b.circ)
   *          yield c
   *      if b.condition yield BRANCH
   *      if visited(b.successor(false)) yield GOTO
   * yield STOP
   */
  if (stage_ == ComItStage::Final) {
    // Just output STOP; move to end()
    *this = CommandIterator();
    return *this;
  }
  while (current_block_ != prog_->block_end()) {
    switch (stage_) {
      case ComItStage::Label: {
        // Just reached a new block; possibly add a label
        stage_ = ComItStage::FirstCommand;
        FGVert block = *current_block_;
        FGEdgeVec ins = prog_->get_in_edges(block);
        if (ins.size() != 1 || prog_->get_source(ins.front()) != prev_block_ ||
            prog_->get_branch(ins.front())) {
          // Add LABEL
          std::string label = get_label(block);
          Op_ptr op = std::make_shared<FlowOp>(OpType::Label, label);
          current_command_ = Command(op, {});
          return *this;
        } else
          continue;
      }
      case ComItStage::FirstCommand: {
        // Just output a label; start iterating through commands of block
        current_com_iterator_ = current_block_.get_circuit_ref().begin();
        if (current_com_iterator_ == current_block_.get_circuit_ref().end()) {
          stage_ = ComItStage::Branch;
          continue;
        } else {
          stage_ = ComItStage::Command;
          current_command_ = *current_com_iterator_;
          return *this;
        }
      }
      case ComItStage::Command: {
        // Just output a command; continue
        ++current_com_iterator_;
        if (current_com_iterator_ == current_block_.get_circuit_ref().end()) {
          stage_ = ComItStage::Branch;
          continue;
        } else {
          current_command_ = *current_com_iterator_;
          return *this;
        }
      }
      case ComItStage::Branch: {
        // Finished commands for block; possibly add a branch
        FGVert block = *current_block_;
        std::optional<Bit> condition = prog_->get_condition(block);
        stage_ = ComItStage::Goto;
        if (condition) {
          // Add BRANCH
          FGVert target = prog_->get_branch_successor(block, true);
          std::string label = get_label(target);
          Op_ptr op = std::make_shared<FlowOp>(OpType::Branch, label);
          current_command_ = Command(op, {*condition});
          return *this;
        } else
          continue;
      }
      case ComItStage::Goto: {
        // Reached end of block; possibly add a Goto before moving to next block
        prev_block_ = *current_block_;
        ++current_block_;
        stage_ = ComItStage::Label;
        FGVert target = prog_->get_branch_successor(prev_block_, false);
        if (current_block_ == prog_->block_end()) {
          if (target == prog_->exit_) continue;
        } else if (target == *current_block_)
          continue;
        // Add Goto
        std::string label = get_label(target);
        Op_ptr op = std::make_shared<FlowOp>(OpType::Goto, label);
        current_command_ = Command(op, {});
        return *this;
      }
      default: {
        // None of Stop, Final, or End should be hit before
        // the block iterator reaches end
        throw ProgramError(
            "Error in command iteration: hit final stages before "
            "reaching exit block");
      }
    }
  }
  // Label to STOP
  if (stage_ == ComItStage::Label) {
    std::map<FGVert, std::string>::iterator found = labels_.find(prog_->exit_);
    if (found != labels_.end()) {
      Op_ptr op = std::make_shared<FlowOp>(OpType::Label, found->second);
      current_command_ = Command(op, {});
      stage_ = ComItStage::Stop;
      return *this;
    }
  }
  // STOP command
  Op_ptr op = std::make_shared<FlowOp>(OpType::Stop);
  current_command_ = Command(op, {});
  stage_ = ComItStage::Final;
  return *this;
}

std::string Program::CommandIterator::get_label(const FGVert &block) {
  std::map<FGVert, std::string>::iterator found_label = labels_.find(block);
  if (found_label == labels_.end()) {
    std::optional<std::string> label = prog_->get_label(block);
    if (!label) label = "lab_" + std::to_string(labels_.size());
    labels_.insert({block, *label});
    return *label;
  } else {
    return found_label->second;
  }
}

Program::CommandIterator Program::begin() const {
  return CommandIterator(*this);
}

Program::CommandIterator Program::end() { return CommandIterator(); }

}  // namespace tket
