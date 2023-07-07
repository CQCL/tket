self: super: {
  symengine = super.symengine.overrideAttrs (old: {
    patches = [ ./symengine.patch ];
    cmakeFlags = old.cmakeFlags ++ [
      "-DBUILD_TESTS=OFF"
      "-DBUILD_BENCHMARKS=OFF"
      "-DBUILD_SHARED_LIBS=ON"
      "-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=yes"
    ];
    propagatedBuildInputs = [
      super.flint
      super.gmp
      super.libmpc
      super.mpfr
    ];
  });
}
