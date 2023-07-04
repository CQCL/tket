{
  pkgs ? import <nixpkgs> {},
  static ? true,
  tklog ? import ../tklog { inherit pkgs static; }
}:
let
  src = builtins.filterSource(p: _: baseNameOf p != "default.nix") ./.;
in
  pkgs.stdenv.mkDerivation{
    name = "tkassert";
    inherit src;
    buildInputs = [ pkgs.cmake ];
    nativeBuildInputs = [ tklog ];
    cmakeFlags = if static then [ "-DBUILD_SHARED_LIBS=OFF" ] else [ "-DBUILD_SHARED_LIBS=ON" ];
    postFixup = ''
      # fix bogus include paths
      # trick found here: https://github.com/NixOS/nixpkgs/blob/master/pkgs/development/libraries/crc32c/default.nix
      for f in $(find $out/lib/cmake -name '*.cmake'); do
        substituteInPlace "$f" --replace "\''${_IMPORT_PREFIX}/$out/include" "\''${_IMPORT_PREFIX}/include"
      done
    '';
  }
