#include "ZX/ZXGenerator.hpp"

#include <sstream>

#include "Utils/Assert.hpp"
#include "ZX/ZXDiagram.hpp"

namespace tket {

namespace zx {

/**
 * ZXType CLASS SETS
 */

bool find_in_set(const ZXType& val, const ZXTypeSet& set) {
  return set.find(val) != set.end();
}

bool is_boundary_type(ZXType type) {
  static const ZXTypeSet boundaries = {
      ZXType::Input, ZXType::Output, ZXType::Open};
  return find_in_set(type, boundaries);
}

bool is_basic_gen_type(ZXType type) {
  static const ZXTypeSet basics = {
      ZXType::ZSpider, ZXType::XSpider, ZXType::Hbox};
  return find_in_set(type, basics);
}

bool is_spider_type(ZXType type) {
  static const ZXTypeSet spiders = {ZXType::ZSpider, ZXType::XSpider};
  return find_in_set(type, spiders);
}

bool is_directed_type(ZXType type) {
  static const ZXTypeSet directed = {ZXType::Triangle, ZXType::ZXBox};
  return find_in_set(type, directed);
}

/**
 * ZXGen (Base class) implementation
 */

ZXGen::ZXGen(ZXType type) : type_(type) {}

ZXType ZXGen::get_type() const { return type_; }

bool ZXGen::operator==(const ZXGen& other) const {
  return (this->type_ == other.type_);
}

ZXGen::~ZXGen() {}

ZXGen_ptr ZXGen::create_gen(ZXType type, QuantumType qtype) {
  ZXGen_ptr op;
  switch (type) {
    case ZXType::Input:
    case ZXType::Output:
    case ZXType::Open: {
      op = std::make_shared<const BoundaryGen>(type, qtype);
      break;
    }
    case ZXType::ZSpider:
    case ZXType::XSpider: {
      op = std::make_shared<const BasicGen>(type, 0., qtype);
      break;
    }
    case ZXType::Hbox: {
      op = std::make_shared<const BasicGen>(type, -1., qtype);
      break;
    }
    case ZXType::Triangle: {
      op = std::make_shared<const DirectedGen>(type, qtype);
      break;
    }
    default:
      throw ZXError("Cannot instantiate a ZXGen of the required type");
  }
  return op;
}

ZXGen_ptr ZXGen::create_gen(ZXType type, const Expr& param, QuantumType qtype) {
  ZXGen_ptr op;
  switch (type) {
    case ZXType::ZSpider:
    case ZXType::XSpider: {
      op = std::make_shared<const BasicGen>(type, param, qtype);
      break;
    }
    case ZXType::Hbox: {
      op = std::make_shared<const BasicGen>(type, param, qtype);
      break;
    }
    default:
      throw ZXError(
          "Cannot instantiate a parameterised ZXGen of the required "
          "type");
  }
  return op;
}

/**
 * BoundaryGen implementation
 */

BoundaryGen::BoundaryGen(ZXType type, QuantumType qtype)
    : ZXGen(type), qtype_(qtype) {
  if (!is_boundary_type(type)) {
    throw ZXError("Unsupported ZXType for BoundaryGen");
  }
}

std::optional<QuantumType> BoundaryGen::get_qtype() const { return qtype_; }

bool BoundaryGen::valid_edge(
    std::optional<unsigned> port, QuantumType qtype) const {
  return !port && (qtype == this->qtype_);
}

SymSet BoundaryGen::free_symbols() const { return {}; }

ZXGen_ptr BoundaryGen::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return ZXGen_ptr();
}

std::string BoundaryGen::get_name(bool) const {
  std::stringstream st;
  if (qtype_ == QuantumType::Quantum) {
    st << "Q-";
  } else {
    st << "C-";
  }
  switch (type_) {
    case ZXType::Input:
      st << "Input";
      break;
    case ZXType::Output:
      st << "Output";
      break;
    case ZXType::Open:
      st << "Open";
      break;
    default:
      throw ZXError("BoundaryGen with invalid ZXType");
  }
  return st.str();
}

bool BoundaryGen::operator==(const ZXGen& other) const {
  if (!ZXGen::operator==(other)) return false;
  const BoundaryGen& other_bgen = static_cast<const BoundaryGen&>(other);
  return (this->qtype_ == other_bgen.qtype_);
}

/**
 * BasicGen implementation
 */

BasicGen::BasicGen(ZXType type, const Expr& param, QuantumType qtype)
    : ZXGen(type), qtype_(qtype), param_(param) {
  if (!is_basic_gen_type(type)) {
    throw ZXError("Unsupported ZXType for BasicGen");
  }
}

std::optional<QuantumType> BasicGen::get_qtype() const { return qtype_; }

bool BasicGen::valid_edge(
    std::optional<unsigned> port, QuantumType qtype) const {
  return !port && (qtype == QuantumType::Quantum ||
                   this->qtype_ == QuantumType::Classical);
}

Expr BasicGen::get_param() const { return param_; }

SymSet BasicGen::free_symbols() const { return expr_free_symbols(param_); }

ZXGen_ptr BasicGen::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  return std::make_shared<const BasicGen>(type_, param_.subs(sub_map), qtype_);
}

std::string BasicGen::get_name(bool) const {
  std::stringstream st;
  if (qtype_ == QuantumType::Quantum) {
    st << "Q-";
  } else {
    st << "C-";
  }
  switch (type_) {
    case ZXType::ZSpider:
      st << "Z";
      break;
    case ZXType::XSpider:
      st << "X";
      break;
    case ZXType::Hbox:
      st << "H";
      break;
    default:
      throw ZXError("BasicGen with invalid ZXType");
  }
  st << "(" << param_ << ")";
  return st.str();
}

bool BasicGen::operator==(const ZXGen& other) const {
  if (!ZXGen::operator==(other)) return false;
  const BasicGen& other_basic = static_cast<const BasicGen&>(other);
  return (
      this->qtype_ == other_basic.qtype_ && this->param_ == other_basic.param_);
}

/**
 * ZXDirected (abstract intermediate class) implementation
 */

ZXDirected::ZXDirected(ZXType type) : ZXGen(type) {
  if (!is_directed_type(type)) {
    throw ZXError("Unsupported ZXType for ZXDirected");
  }
}

/**
 * DirectedGen implementation
 */

DirectedGen::DirectedGen(ZXType type, QuantumType qtype)
    : ZXDirected(type), qtype_(qtype) {
  if (type != ZXType::Triangle) {
    throw ZXError("Unsupported ZXType for DirectedGen");
  }
}

std::optional<QuantumType> DirectedGen::get_qtype() const { return qtype_; }

bool DirectedGen::valid_edge(
    std::optional<unsigned> port, QuantumType qtype) const {
  return port && (*port < this->n_ports()) && (qtype == this->qtype_);
}

SymSet DirectedGen::free_symbols() const { return {}; }

ZXGen_ptr DirectedGen::symbol_substitution(
    const SymEngine::map_basic_basic&) const {
  return ZXGen_ptr();
}

std::string DirectedGen::get_name(bool) const {
  if (qtype_ == QuantumType::Quantum) {
    return "Q-Tri";
  } else {
    return "C-Tri";
  }
}

bool DirectedGen::operator==(const ZXGen& other) const {
  if (!ZXGen::operator==(other)) return false;

  const DirectedGen& other_dir = static_cast<const DirectedGen&>(other);
  return (this->qtype_ == other_dir.qtype_);
}

unsigned DirectedGen::n_ports() const { return 2; }

std::vector<QuantumType> DirectedGen::get_signature() const {
  return std::vector<QuantumType>(2, qtype_);
}

/**
 * ZXBox implementation
 */

ZXBox::ZXBox(const ZXDiagram& diag)
    : ZXDirected(ZXType::ZXBox),
      diag_(std::make_shared<const ZXDiagram>(diag)) {}

std::shared_ptr<const ZXDiagram> ZXBox::get_diagram() const { return diag_; }

std::optional<QuantumType> ZXBox::get_qtype() const { return std::nullopt; }

bool ZXBox::valid_edge(std::optional<unsigned> port, QuantumType qtype) const {
  if (!port) return false;
  ZXVertVec boundary = diag_->get_boundary();
  if (*port >= boundary.size()) return false;
  return (qtype == diag_->get_qtype(boundary.at(*port)));
}

SymSet ZXBox::free_symbols() const { return diag_->free_symbols(); }

ZXGen_ptr ZXBox::symbol_substitution(
    const SymEngine::map_basic_basic& sub_map) const {
  ZXDiagram new_diag(*diag_);
  new_diag.symbol_substitution(sub_map);
  return std::make_shared<const ZXBox>(new_diag);
}

std::string ZXBox::get_name(bool) const { return "Box"; }

bool ZXBox::operator==(const ZXGen&) const {
  // Checking for a proper graph isomorphism is difficult. Safest to just assume
  // all boxes are unique.
  return false;
}

unsigned ZXBox::n_ports() const { return diag_->get_boundary().size(); }

std::vector<QuantumType> ZXBox::get_signature() const {
  ZXVertVec boundary = diag_->get_boundary();
  std::vector<QuantumType> sig;
  for (const ZXVert& b : boundary) {
    std::optional<QuantumType> qt = diag_->get_qtype(b);
    TKET_ASSERT(qt.has_value());
    sig.push_back(*qt);
  }
  return sig;
}

}  // namespace zx

}  // namespace tket
