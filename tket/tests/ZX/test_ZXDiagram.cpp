#include <catch2/catch.hpp>
#include <fstream>

#include "ZX/ZXDiagram.hpp"
#include "ZX/ZXGenerator.hpp"

namespace tket {
namespace zx {
namespace test_ZXDiagram {

SCENARIO("Testing generator creation") {
  BoundaryGen input(ZXType::Input, QuantumType::Quantum);
  CHECK(input.get_name() == "Q-Input");
  CHECK(input.get_type() == ZXType::Input);
  CHECK(input.get_qtype() == QuantumType::Quantum);
  CHECK(input.valid_edge(std::nullopt, QuantumType::Quantum));
  CHECK_FALSE(input.valid_edge(0, QuantumType::Quantum));
  CHECK_FALSE(input.valid_edge(std::nullopt, QuantumType::Classical));

  BasicGen zSpider(ZXType::ZSpider, 0.3, QuantumType::Classical);
  CHECK(zSpider.get_name() == "C-Z(0.3)");
  CHECK(zSpider.get_type() == ZXType::ZSpider);
  CHECK(zSpider.get_qtype() == QuantumType::Classical);
  CHECK(zSpider.free_symbols().empty());
  CHECK(zSpider.valid_edge(std::nullopt, QuantumType::Quantum));
  CHECK(zSpider.valid_edge(std::nullopt, QuantumType::Classical));
  CHECK_FALSE(zSpider.valid_edge(0, QuantumType::Quantum));

  BasicGen xSpider(ZXType::XSpider, Expr("2*a"), QuantumType::Quantum);
  CHECK(xSpider.get_name() == "Q-X(2*a)");
  CHECK(xSpider.get_type() == ZXType::XSpider);
  CHECK(xSpider.get_qtype() == QuantumType::Quantum);
  CHECK(xSpider.free_symbols().size() == 1);
  CHECK(xSpider.valid_edge(std::nullopt, QuantumType::Quantum));
  CHECK_FALSE(xSpider.valid_edge(std::nullopt, QuantumType::Classical));
  SymEngine::map_basic_basic sub_map;
  Sym a = SymEngine::symbol("a");
  sub_map[a] = Expr(0.8);
  CHECK(xSpider.symbol_substitution(sub_map)->get_name() == "Q-X(1.6)");

  // Should throw an error: type Triangle is not a BasicGen type
  REQUIRE_THROWS_AS(BasicGen(ZXType::Triangle, 0.3), ZXError);

  DirectedGen tri(ZXType::Triangle, QuantumType::Classical);
  CHECK(tri.get_name() == "C-Tri");
  CHECK_FALSE(tri.valid_edge(std::nullopt, QuantumType::Classical));
  CHECK_FALSE(tri.valid_edge(2, QuantumType::Classical));
  CHECK(tri.valid_edge(0, QuantumType::Classical));
  CHECK_FALSE(tri.valid_edge(1, QuantumType::Quantum));
}

SCENARIO("Testing diagram creation & vertex/edge additions") {
  ZXDiagram diag(1, 1, 0, 0);
  CHECK(diag.get_scalar() == Expr(1.));
  CHECK_FALSE(diag.is_symbolic());
  diag.multiply_scalar(0.4);
  diag.multiply_scalar(Expr("2*a"));
  CHECK(diag.get_scalar() == Expr("0.8*a"));
  CHECK(diag.free_symbols().size() == 1);

  ZXVert zSpid_v = diag.add_vertex(ZXType::ZSpider, 0.1);
  ZXVert xSpid_v = diag.add_vertex(ZXType::XSpider, 3.4);
  ZXVert zSpid2_v =
      diag.add_vertex(ZXType::ZSpider, 6.7, QuantumType::Classical);

  REQUIRE_THROWS_AS(diag.add_vertex(ZXType::ZXBox, 3.), ZXError);

  REQUIRE_THROWS_WITH(
      diag.check_validity(), "Boundary vertex does not have degree 1");

  diag.add_wire(diag.get_boundary().at(0), zSpid_v);
  diag.add_wire(diag.get_boundary().at(1), xSpid_v);
  diag.add_wire(zSpid_v, xSpid_v);
  diag.add_wire(xSpid_v, zSpid_v, ZXWireType::H);
  Wire extra = diag.add_wire(diag.get_boundary().at(1), zSpid_v);
  REQUIRE_THROWS_WITH(
      diag.check_validity(), "Boundary vertex does not have degree 1");

  diag.remove_wire(extra);
  REQUIRE_NOTHROW(diag.check_validity());

  Wire wrong_port = diag.add_wire(
      zSpid2_v, xSpid_v, ZXWireType::Basic, QuantumType::Quantum, 0);
  REQUIRE_THROWS_WITH(
      diag.check_validity(), "Wire at a named port of an undirected vertex");
  diag.remove_wire(wrong_port);

  ZXVert tri_v = diag.add_vertex(ZXType::Triangle);
  diag.add_wire(tri_v, zSpid_v, ZXWireType::Basic, QuantumType::Quantum, 0);
  REQUIRE_THROWS_WITH(
      diag.check_validity(),
      "Not all ports of a directed vertex have wires connected");

  Wire no_port = diag.add_wire(zSpid_v, tri_v);
  REQUIRE_THROWS_WITH(
      diag.check_validity(), "Wire at an unnamed port of a directed vertex");

  diag.remove_wire(no_port);
  diag.add_wire(
      zSpid_v, tri_v, ZXWireType::Basic, QuantumType::Quantum, std::nullopt, 1);
  REQUIRE_NOTHROW(diag.check_validity());

  Wire extra_port =
      diag.add_wire(tri_v, zSpid_v, ZXWireType::Basic, QuantumType::Quantum, 1);
  REQUIRE_THROWS_WITH(
      diag.check_validity(), "Multiple wires on the same port of a vertex");
  diag.remove_wire(extra_port);

  ZXDiagram inner(1, 2, 1, 0);
  ZXVert inner_spid =
      inner.add_vertex(ZXType::ZSpider, 0.6, QuantumType::Classical);
  inner.add_wire(inner_spid, inner.get_boundary().at(0));
  inner.add_wire(inner_spid, inner.get_boundary().at(1));
  inner.add_wire(inner_spid, inner.get_boundary().at(2), ZXWireType::H);
  inner.add_wire(
      inner_spid, inner.get_boundary().at(3), ZXWireType::Basic,
      QuantumType::Classical);
  ZXGen_ptr box = std::make_shared<const ZXBox>(inner);

  ZXVert box_v = diag.add_vertex(box);
  diag.add_wire(box_v, zSpid2_v, ZXWireType::Basic, QuantumType::Quantum, 0);
  diag.add_wire(box_v, zSpid2_v, ZXWireType::Basic, QuantumType::Quantum, 1);
  diag.add_wire(box_v, xSpid_v, ZXWireType::Basic, QuantumType::Quantum, 2);
  Wire wrong_qtype = diag.add_wire(
      box_v, zSpid2_v, ZXWireType::Basic, QuantumType::Quantum, 3);
  REQUIRE_THROWS_WITH(
      diag.check_validity(),
      "QuantumType of wire is incompatible with the given port");

  diag.set_wire_qtype(wrong_qtype, QuantumType::Classical);
  REQUIRE_NOTHROW(diag.check_validity());

  THEN("Print diagram to file") {
    std::ofstream dot_file("zxdiag.dot");
    dot_file << diag.to_graphviz_str();
    dot_file.close();
    remove("zxdiag.dot");
  }
}

}  // namespace test_ZXDiagram
}  // namespace zx
}  // namespace tket
