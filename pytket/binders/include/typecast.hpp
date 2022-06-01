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

#pragma once

#include <pybind11/detail/typeid.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Utils/Exceptions.hpp"
#include "Utils/Expression.hpp"
#include "Utils/Symbols.hpp"
#include "unit_downcast.hpp"

namespace pybind11 {

namespace detail {

template <>
struct type_caster<SymEngine::Expression> {
 public:
  PYBIND11_TYPE_CASTER(SymEngine::Expression, _("Expression"));
  static void assert_tuple_length(tuple t, unsigned len) {
    if (t.size() != len)
      throw std::logic_error("Sympy expression is not well-formed");
  };
  static tuple get_checked_args(handle py_expr, unsigned expected_len) {
    tuple arg_tuple = py_expr.attr("args");
    if (arg_tuple.size() != expected_len) {
      std::stringstream err;
      err << "Expected " << repr(py_expr);
      err << " to have " << expected_len;
      err << " arguments, but it had " << arg_tuple.size();
      throw std::invalid_argument(err.str());
    }
    return arg_tuple;
  }
  static tket::Expr sympy_to_expr(handle py_expr) {
    pybind11::module sympy = pybind11::module::import("sympy");
    handle numbers = sympy.attr("core").attr("numbers");

    if (isinstance(py_expr, sympy.attr("Symbol"))) {
      handle name = py_expr.attr("name");
      tket::Sym sym = SymEngine::symbol(name.cast<std::string>());
      return tket::Expr(sym);
    } else if (isinstance(py_expr, sympy.attr("Mul"))) {
      tuple arg_tuple = py_expr.attr("args");
      tket::Expr res(1);
      for (handle elem : arg_tuple) {
        res *= sympy_to_expr(elem);
      }
      return res;
    } else if (isinstance(py_expr, sympy.attr("Add"))) {
      tuple arg_tuple = py_expr.attr("args");
      tket::Expr res(0);
      for (handle elem : arg_tuple) {
        res += sympy_to_expr(elem);
      }
      return res;
    } else if (isinstance(py_expr, sympy.attr("Pow"))) {
      tuple arg_tuple = get_checked_args(py_expr, 2);
      return SymEngine::pow(
          sympy_to_expr(arg_tuple[0]), sympy_to_expr(arg_tuple[1]));
    } else if (isinstance(py_expr, sympy.attr("Integer"))) {
      return tket::Expr(py_expr.attr("p").cast<long>());
    } else if (isinstance(py_expr, sympy.attr("Rational"))) {
      return tket::Expr(py_expr.attr("p").cast<long>()) /
             tket::Expr(py_expr.attr("q").cast<long>());
    } else if (isinstance(py_expr, sympy.attr("Float"))) {
      return tket::Expr(repr(py_expr));
    } else if (isinstance(py_expr, numbers.attr("ImaginaryUnit"))) {
      return tket::Expr(SymEngine::I);
    } else if (isinstance(py_expr, numbers.attr("Exp1"))) {
      return tket::Expr(SymEngine::E);
    } else if (isinstance(py_expr, numbers.attr("Pi"))) {
      return tket::Expr(SymEngine::pi);
    } else if (isinstance(py_expr, numbers.attr("NegativeInfinity"))) {
      return tket::Expr(SymEngine::NegInf);
    } else if (isinstance(py_expr, numbers.attr("Infinity"))) {
      return tket::Expr(SymEngine::Inf);
    } else if (isinstance(py_expr, numbers.attr("ComplexInfinity"))) {
      return tket::Expr(SymEngine::ComplexInf);
    } else if (isinstance(py_expr, numbers.attr("NaN"))) {
      return tket::Expr(SymEngine::Nan);
    }
// functions with a single argument that have an equivalent in symengine
#define SECONVERT(engmeth, pyclass)                     \
  else if (isinstance(py_expr, sympy.attr(#pyclass))) { \
    tuple arg_tuple = get_checked_args(py_expr, 1);     \
    tket::Expr the_arg = sympy_to_expr(arg_tuple[0]);   \
    tket::ExprPtr res = SymEngine::engmeth(the_arg);    \
    return tket::Expr(res);                             \
  }
    SECONVERT(log, log)
    SECONVERT(conjugate, conjugate)
    SECONVERT(sin, sin)
    SECONVERT(cos, cos)
    SECONVERT(tan, tan)
    SECONVERT(cot, cot)
    SECONVERT(csc, csc)
    SECONVERT(sec, sec)
    SECONVERT(asin, asin)
    SECONVERT(acos, acos)
    SECONVERT(asec, asec)
    SECONVERT(acsc, acsc)
    SECONVERT(atan, atan)
    SECONVERT(acot, acot)
    SECONVERT(sinh, sinh)
    SECONVERT(csch, csch)
    SECONVERT(cosh, cosh)
    SECONVERT(sech, sech)
    SECONVERT(tanh, tanh)
    SECONVERT(coth, coth)
    SECONVERT(asinh, asinh)
    SECONVERT(acsch, acsch)
    SECONVERT(acosh, acosh)
    SECONVERT(atanh, atanh)
    SECONVERT(acoth, acoth)
    SECONVERT(asech, asech)
    SECONVERT(erf, erf)
    SECONVERT(erfc, erfc)
    SECONVERT(abs, Abs)
#undef SECONVERT
    else if (isinstance(py_expr, sympy.attr("atan2"))) {
      tuple arg_tuple = get_checked_args(py_expr, 2);
      return SymEngine::atan2(
          sympy_to_expr(arg_tuple[0]), sympy_to_expr(arg_tuple[1]));
    }
    else {
      std::stringstream err;
      err << "Unable to convert sympy expression " << repr(py_expr);
      throw tket::NotImplemented(err.str());
    }
  }
  bool load(handle src, bool) {
    pybind11::module sympy = pybind11::module::import("sympy");
    if (isinstance(src, sympy.attr("Expr"))) {
      value = sympy_to_expr(src);
      return true;
    } else {
      double v = PyFloat_AsDouble(src.ptr());
      if (!PyErr_Occurred()) {
        value = SymEngine::Expression(v);
        return true;
      }
      PyErr_Clear();
    }
    return false;
  }
  static object basic_to_sympy(const tket::ExprPtr& e_) {
    pybind11::module sympy = pybind11::module::import("sympy");
    switch (e_->get_type_code()) {
      case SymEngine::TypeID::SYMENGINE_SYMBOL: {
        const SymEngine::Symbol* s =
            dynamic_cast<const SymEngine::Symbol*>(e_.get());
        return sympy.attr("Symbol")(s->get_name());
      }
      case SymEngine::TypeID::SYMENGINE_MUL: {
        const SymEngine::Mul* m = dynamic_cast<const SymEngine::Mul*>(e_.get());
        object res = basic_to_sympy(m->get_coef());
        const SymEngine::map_basic_basic d = m->get_dict();
        for (SymEngine::map_basic_basic::const_iterator i = d.begin();
             i != d.end(); i++) {
          res = res * sympy.attr("Pow")(
                          basic_to_sympy(i->first), basic_to_sympy(i->second));
        }
        return res;
      }
      case SymEngine::TypeID::SYMENGINE_ADD: {
        const SymEngine::Add* a = dynamic_cast<const SymEngine::Add*>(e_.get());
        object res = basic_to_sympy(a->get_coef());
        const SymEngine::umap_basic_num d = a->get_dict();
        for (SymEngine::umap_basic_num::const_iterator i = d.begin();
             i != d.end(); i++) {
          res = res + basic_to_sympy(i->first) * basic_to_sympy(i->second);
        }
        return res;
      }
      case SymEngine::TypeID::SYMENGINE_POW: {
        const SymEngine::Pow* p = dynamic_cast<const SymEngine::Pow*>(e_.get());
        return sympy.attr("Pow")(
            basic_to_sympy(p->get_base()), basic_to_sympy(p->get_exp()));
      }
      case SymEngine::TypeID::SYMENGINE_INTEGER: {
        const SymEngine::Integer* i =
            dynamic_cast<const SymEngine::Integer*>(e_.get());
        return reinterpret_borrow<object>(PyLong_FromLong(i->as_int()));
      }
      case SymEngine::TypeID::SYMENGINE_RATIONAL: {
        const SymEngine::Rational* r =
            dynamic_cast<const SymEngine::Rational*>(e_.get());
        return sympy.attr("Rational")(
            basic_to_sympy(r->get_num()), basic_to_sympy(r->get_den()));
      }
      case SymEngine::TypeID::SYMENGINE_REAL_DOUBLE: {
        const SymEngine::RealDouble* d =
            dynamic_cast<const SymEngine::RealDouble*>(e_.get());
        return reinterpret_borrow<object>(PyFloat_FromDouble(d->as_double()));
      }
      case SymEngine::TypeID::SYMENGINE_COMPLEX:
      case SymEngine::TypeID::SYMENGINE_COMPLEX_DOUBLE: {
        const SymEngine::ComplexBase* c =
            dynamic_cast<const SymEngine::ComplexBase*>(e_.get());
        return basic_to_sympy(c->real_part()) +
               ((object)sympy.attr("I")) * basic_to_sympy(c->imaginary_part());
      }
      case SymEngine::TypeID::SYMENGINE_CONSTANT: {
        const SymEngine::Constant* c =
            dynamic_cast<const SymEngine::Constant*>(e_.get());
        std::string name = c->get_name();
        if (name == "E") {
          return sympy.attr("E");
        } else if (name == "pi") {
          return sympy.attr("pi");
        } else {
          throw tket::NotImplemented(
              "Unable to convert SymEngine constant " + name);
        }
      }
      case SymEngine::TypeID::SYMENGINE_INFTY: {
        const SymEngine::Infty* i =
            dynamic_cast<const SymEngine::Infty*>(e_.get());
        if (i->is_positive()) {
          return sympy.attr("oo");
        } else if (i->is_negative()) {
          return -sympy.attr("oo");
        } else {
          return sympy.attr("zoo");
        }
      }
      case SymEngine::TypeID::SYMENGINE_NOT_A_NUMBER: {
        return sympy.attr("nan");
      }
// functions with a single argument that have an equivalent in sympy
#define SYMCONVERT(symid, symclass, sympyname)                   \
  case SymEngine::TypeID::SYMENGINE_##symid: {                   \
    const SymEngine::symclass* x =                               \
        dynamic_cast<const SymEngine::symclass*>(e_.get());      \
    return sympy.attr(#sympyname)(basic_to_sympy(x->get_arg())); \
  }
        SYMCONVERT(LOG, Log, log)
        SYMCONVERT(CONJUGATE, Conjugate, conjugate)
        SYMCONVERT(SIN, Sin, sin)
        SYMCONVERT(COS, Cos, cos)
        SYMCONVERT(TAN, Tan, tan)
        SYMCONVERT(COT, Cot, cot)
        SYMCONVERT(CSC, Csc, csc)
        SYMCONVERT(SEC, Sec, sec)
        SYMCONVERT(ASIN, ASin, asin)
        SYMCONVERT(ACOS, ACos, acos)
        SYMCONVERT(ASEC, ASec, asec)
        SYMCONVERT(ACSC, ACsc, acsc)
        SYMCONVERT(ATAN, ATan, atan)
        SYMCONVERT(ACOT, ACot, acot)
        SYMCONVERT(SINH, Sinh, sinh)
        SYMCONVERT(CSCH, Csch, csch)
        SYMCONVERT(COSH, Cosh, cosh)
        SYMCONVERT(SECH, Sech, sech)
        SYMCONVERT(TANH, Tanh, tanh)
        SYMCONVERT(COTH, Coth, coth)
        SYMCONVERT(ASINH, ASinh, asinh)
        SYMCONVERT(ACSCH, ACsch, acsch)
        SYMCONVERT(ACOSH, ACosh, acosh)
        SYMCONVERT(ATANH, ATanh, atanh)
        SYMCONVERT(ACOTH, ACoth, acoth)
        SYMCONVERT(ASECH, ASech, asech)
        SYMCONVERT(ERF, Erf, erf)
        SYMCONVERT(ERFC, Erfc, erfc)
        SYMCONVERT(ABS, Abs, Abs)
#undef SYMCONVERT
      case SymEngine::TypeID::SYMENGINE_ATAN2: {
        const SymEngine::ATan2* a =
            dynamic_cast<const SymEngine::ATan2*>(e_.get());
        return sympy.attr("atan2")(
            basic_to_sympy(a->get_num()), basic_to_sympy(a->get_den()));
      }
      default: {
        std::stringstream err;
        err << "Unable to convert SymEngine expression " << e_->__str__();
        throw tket::NotImplemented(err.str());
      }
    }
  }
  static handle cast(
      SymEngine::Expression src, return_value_policy /* policy */,
      handle /* parent */) {
    std::optional<double> eval = tket::eval_expr(src);
    if (!eval)
      return basic_to_sympy(src).release();
    else {
      return PyFloat_FromDouble(eval.value());
    }
  }
};

template <>
struct type_caster<SymEngine::RCP<const SymEngine::Symbol>> {
 public:
  PYBIND11_TYPE_CASTER(SymEngine::RCP<const SymEngine::Symbol>, _("Symbol"));
  bool load(handle src, bool) {
    pybind11::module sympy = pybind11::module::import("sympy");
    if (!isinstance(src, sympy.attr("Symbol"))) return false;
    value = SymEngine::symbol(repr(src));
    return true;
  }
  static handle cast(
      SymEngine::RCP<const SymEngine::Symbol> src,
      return_value_policy /* policy */, handle /* parent */) {
    pybind11::module sympy = pybind11::module::import("sympy");
    return sympy.attr("Symbol")(src->get_name()).release();
  }
};
}  // namespace detail
}  // namespace pybind11
