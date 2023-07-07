self: super: {
  symengine-static = super.symengine.overrideAttrs (old: {
    patches = [ ./symengine.patch ];
    cmakeFlags = old.cmakeFlags ++ ["-DBUILD_TESTS=OFF" "-DBUILD_BENCHMARKS=OFF"];
    propagatedBuildInputs = [
      super.flint
      super.gmp
      super.libmpc
      super.mpfr
    ];
  });
  symengine-shared = self.symengine-static.overrideAttrs (old: {
    cmakeFlags = old.cmakeFlags ++ [
      "-DBUILD_SHARED_LIBS=ON"
      "-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=yes"
    ];
  });
}
