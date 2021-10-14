#include "ZX/Rewrite.hpp"

namespace tket {

namespace zx {

Rewrite::Rewrite(const RewriteFun &fun) : apply(fun) {}

Rewrite Rewrite::sequence(const std::vector<Rewrite> &rvec) {
  return Rewrite([=](ZXDiagram &diag) {
    bool success = false;
    for (std::vector<Rewrite>::const_iterator it = rvec.begin();
         it != rvec.end(); ++it) {
      success = it->apply(diag) || success;
    }
    return success;
  });
}

Rewrite Rewrite::repeat(const Rewrite &rw) {
  return Rewrite([=](ZXDiagram &diag) {
    bool success = false;
    while (rw.apply(diag)) success = true;
    return success;
  });
}

Rewrite Rewrite::repeat_with_metric(
    const Rewrite &rw, const Rewrite::Metric &eval) {
  return Rewrite([=](ZXDiagram &diag) {
    bool success = false;
    unsigned currentVal = eval(diag);
    ZXDiagram newDiag = diag;
    rw.apply(newDiag);
    unsigned newVal = eval(newDiag);
    while (newVal < currentVal) {
      currentVal = newVal;
      success = true;
      rw.apply(newDiag);
      newVal = eval(newDiag);
    }
    if (success) {
      diag = newDiag;
    }
    return success;
  });
}

Rewrite Rewrite::repeat_while(const Rewrite &cond, const Rewrite &body) {
  return Rewrite([=](ZXDiagram &diag) {
    bool success = false;
    while (cond.apply(diag)) {
      success = true;
      body.apply(diag);
    }
    return success;
  });
}

}  // namespace zx

}  // namespace tket
