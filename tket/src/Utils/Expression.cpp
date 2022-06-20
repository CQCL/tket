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

#include "Expression.hpp"

#include <optional>

#include "Constants.hpp"
#include "Symbols.hpp"
#include "symengine/symengine_exception.h"

namespace tket {

bool approx_0(const Expr& e, double tol) {
  std::optional<double> v = eval_expr(e);
  return v && (std::abs(v.value()) < tol);
}

double fmodn(double x, unsigned n) {
  x /= n;
  x -= floor(x);
  return n * x;
}

bool approx_eq(double x, double y, unsigned mod, double tol) {
  double r = fmodn(x - y, mod);
  return r < tol || r > mod - tol;
}

SymSet expr_free_symbols(const Expr& e) {
  SymSet symbols;
  for (auto x : SymEngine::free_symbols(e)) {
    symbols.insert(SymEngine::rcp_static_cast<const SymEngine::Symbol>(x));
  }
  return symbols;
}

SymSet expr_free_symbols(const std::vector<Expr>& es) {
  SymSet symbols;
  for (auto e : es) {
    for (auto x : SymEngine::free_symbols(e)) {
      symbols.insert(SymEngine::rcp_static_cast<const SymEngine::Symbol>(x));
    }
  }
  return symbols;
}

std::optional<double> eval_expr(const Expr& e) {
  if (!SymEngine::free_symbols(e).empty()) {
    return std::nullopt;
  } else {
    try {
      return SymEngine::eval_double(e);
    } catch (SymEngine::NotImplementedError&) {
      return std::nullopt;
    }
  }
}

std::optional<Complex> eval_expr_c(const Expr& e) {
  if (!SymEngine::free_symbols(e).empty()) {
    return std::nullopt;
  } else {
    return SymEngine::eval_complex_double(e);
  }
}

std::optional<double> eval_expr_mod(const Expr& e, unsigned n) {
  std::optional<double> reduced_val = eval_expr(e);
  if (!reduced_val) return std::nullopt;
  double val = reduced_val.value();
  double val4 = 4 * val;
  long nearest_val4 = std::lrint(val4);
  if (std::abs(val4 - nearest_val4) < 4 * EPS) {
    val = nearest_val4 * 0.25;
  }
  return fmodn(val, n);
}

// Evaluate cos(pi x / 12). If x is close to a multiple of pi/12, it is clamped
// to that exact multiple and the return value is exact.
// Is is assumed that 0 <= x < 24.
static Expr cos_pi_by_12_times(double x) {
  static const double pi_by_12 = PI / 12;
  static const Expr pi_by_12_sym =
      SymEngine::div(SymEngine::pi, SymEngine::integer(12));
  int n(x + 0.5);  // nearest integer to x (assuming x >= 0)
  if (std::abs(x - n) < EPS) {
    Expr z = SymEngine::cos(n * pi_by_12_sym);
    return z;
  } else {
    return Expr(cos(pi_by_12 * x));
  }
}

Expr cos_halfpi_times(const Expr& e) {
  std::optional<double> x = eval_expr_mod(e / 2);  // >= 0
  if (x) {
    return cos_pi_by_12_times(12 * x.value());
  } else {
    return SymEngine::cos(SymEngine::expand(e * SymEngine::pi / 2));
  }
}

Expr sin_halfpi_times(const Expr& e) {
  return cos_halfpi_times(SymEngine::expand(SymEngine::integer(1) - e));
}

Expr minus_times(const Expr& e) {
  Expr e1 = -e;
  Expr e2 = SymEngine::expand(e1);

  unsigned e1_size = e1.get_basic()->dumps().size();
  unsigned e2_size = e2.get_basic()->dumps().size();

  return e2_size < e1_size ? e2 : e1;
}

bool equiv_expr(const Expr& e0, const Expr& e1, unsigned n, double tol) {
  std::optional<double> eval0 = eval_expr(e0);
  std::optional<double> eval1 = eval_expr(e1);
  if (!eval0 || !eval1) return e0 == e1;
  return approx_eq(eval0.value(), eval1.value(), n, tol);
}

bool equiv_val(const Expr& e, double x, unsigned n, double tol) {
  std::optional<double> eval = eval_expr(e);
  if (!eval) return false;
  return approx_eq(eval.value(), x, n, tol);
}

bool equiv_0(const Expr& e, unsigned n, double tol) {
  return equiv_val(e, 0., n, tol);
}

std::optional<unsigned> equiv_Clifford(const Expr& e, unsigned n, double tol) {
  std::optional<double> eval = eval_expr_mod(e, n);
  if (!eval) return std::nullopt;
  double v_mod_n = eval.value();
  unsigned nearest = lround(v_mod_n * 2);
  if (abs(v_mod_n - (nearest * 0.5)) < tol)
    return nearest;
  else
    return std::nullopt;
}

}  // namespace tket
