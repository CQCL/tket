{
  pkgs ? import <nixpkgs> {},
  static ? true
}:
let
  symengine-static = pkgs.symengine.overrideAttrs (old: {
    # symengine incorrectly marks the following
    # as buildInputs, but they are actually propagatedBuildInputs,
    # as symengine's headers include them directly. Thus when
    # including from symengine, they must be available.
    patches = [
      ./symengine.patch
    ];
    cmakeFlags = old.cmakeFlags ++ [
      "-DBUILD_TESTS=OFF"
      "-DBUILD_BENCHMARKS=OFF"
    ];
    propagatedBuildInputs = [
      pkgs.flint
      pkgs.gmp
      pkgs.libmpc
      pkgs.mpfr
    ];
  });
  symengine-shared = symengine-static.overrideAttrs (old: {
    cmakeFlags = old.cmakeFlags ++ [
      "-DBUILD_SHARED_LIBS=ON"
      "-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=yes"
    ];
  });
in
  builtins.trace "static: ${if static then "true" else "false"}" (
  if static
    then symengine-static
    else symengine-shared)
 
