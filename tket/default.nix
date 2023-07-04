{
  pkgs ? import <nixpkgs> {},
  static ? true,
  libs ? import ../libs { inherit pkgs static; }
}:
let
  src = builtins.filterSource(p: _: baseNameOf p != "default.nix") ./.;
in
  pkgs.stdenv.mkDerivation{
    name = "tket";
    inherit src;
    nativeBuildInputs = [pkgs.cmake];
    buildInputs = [
      pkgs.boost
      pkgs.symengine
      pkgs.eigen
      pkgs.nlohmann_json
      libs.tklog
      libs.tkassert
      libs.tkrng
      libs.tktokenswap
      libs.tkwsm
    ] ++ pkgs.symengine.buildInputs;
    cmakeFlags = if static then [ "-DBUILD_SHARED_LIBS=OFF" ] else [ "-DBUILD_SHARED_LIBS=ON" ];
    postFixup = ''
      # fix bogus include paths
      # trick found here: https://github.com/NixOS/nixpkgs/blob/master/pkgs/development/libraries/crc32c/default.nix
      for f in $(find $out/lib/cmake -name '*.cmake'); do
        substituteInPlace "$f" --replace "\''${_IMPORT_PREFIX}/$out/include" "\''${_IMPORT_PREFIX}/include"
      done
    '';
  }
