#pragma once
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

namespace tket {

namespace zx {

/**
 * Error messages
 **/
class ZXError : public std::logic_error {
 public:
  explicit ZXError(const std::string& message) : std::logic_error(message) {}
};

/**
 * Type for an edge in a ZX diagram.
 *  `Basic`: normal (identity) edge
 *  `H`: Hamadard edge;
 **/
enum class ZXWireType { Basic, H };

/**
 * Provides distinction between classical and quantum diagrams,
 * wires, and operators, relating to the doubling construction for CP-maps.
 **/
enum class QuantumType { Quantum, Classical };

/**
 * Ports are used for non-commutative generators, such as
 * directed generators, to distinguish the incident wires.
 * When needed, ports on the source and target vertex are stored
 * in the wire properties.
 * Since the graphs are semantically undirected, it is not given
 * whether a vertex is a source or target of each incident edge at
 * compile time, so WireEnd variables express which it is at runtime.
 **/
enum class WireEnd { Source, Target };

}  // namespace zx

}  // namespace tket
