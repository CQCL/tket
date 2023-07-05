{
  pkgs ? import <nixpkgs> {},
  static ? true,
  libs ? import ../libs { inherit pkgs static; },
  symengine ? import ../nix-helpers/symengine.nix { inherit pkgs static; }
}:
let
  src = builtins.filterSource(p: _: baseNameOf p != "default.nix") ./.;
  inputs = [ 
      pkgs.boost
      symengine
      pkgs.eigen
      pkgs.nlohmann_json
      libs.tklog
      libs.tkassert
      libs.tkrng
      libs.tktokenswap
      libs.tkwsm
  ] ++ symengine.buildInputs;
  rpath = builtins.concatStringsSep ":" (
    map (x: "${x}/lib")
    [libs.tklog libs.tkassert libs.tkrng libs.tktokenswap libs.tkwsm
    pkgs.gcc.libc pkgs.gcc.cc.lib.lib symengine]
  );
  tket-static = pkgs.stdenv.mkDerivation{
    name = "tket";
    inherit src;
    nativeBuildInputs = [pkgs.cmake];
    buildInputs = inputs;
    postFixup = ''
      # fix bogus include paths
      # trick found here: https://github.com/NixOS/nixpkgs/blob/master/pkgs/development/libraries/crc32c/default.nix
      for f in $(find $out/lib/cmake -name '*.cmake'); do
        substituteInPlace "$f" --replace "\''${_IMPORT_PREFIX}/$out/include" "\''${_IMPORT_PREFIX}/include"
      done
    '';
  };
  tket-shared = pkgs.stdenv.mkDerivation{
    name = "tket";
    inherit src;
    nativeBuildInputs = [pkgs.cmake];
    buildInputs = inputs;
    propagatedBuildInputs = [symengine];
    cmakeFlags = ["-DBUILD_SHARED_LIBS=ON"];
    postFixup = ''
      # fix bogus include paths
      # trick found here:
      for f in $(find $out/lib/cmake -name '*.cmake'); do
        substituteInPlace "$f" --replace "\''${_IMPORT_PREFIX}/$out/include" "\''${_IMPORT_PREFIX}/include"
      done
      patchelf --add-needed "libsymengine.so" $out/lib/libtket.so;
      patchelf --set-rpath "''${ORIGIN}:${rpath}" $out/lib/libtket.so;
    '';
  };
in
  if static then tket-static else tket-shared
