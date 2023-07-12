self: super: {
  # symengine in nixpkgs is built as a static library, and doesn't expose
  # propagated inputs in its pkgconfig file. This is a workaround to build
  # it as a shared library, and to expose the propagated inputs.
  symengine = super.symengine.overrideAttrs (old: {
    patches = [ ./symengine.patch ];
    cmakeFlags = old.cmakeFlags ++ [
      "-DBUILD_TESTS=OFF"
      "-DBUILD_BENCHMARKS=OFF"
      "-DBUILD_SHARED_LIBS=ON"
      "-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=yes"
    ];
    propagatedBuildInputs = [ super.flint super.gmp super.libmpc super.mpfr ];
  });
}
