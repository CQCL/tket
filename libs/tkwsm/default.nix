{
  pkgs ? import <nixpkgs> {},
  tklog ? import ../tklog { inherit pkgs; },
  tkrng ? import ../tkrng { inherit pkgs; },
  tkassert ? import ../tkassert { inherit pkgs; inherit tklog; },
}:
let
  src = builtins.filterSource(p: _: baseNameOf p != "default.nix") ./.;
in
  pkgs.stdenv.mkDerivation{
    name = "tkwsm";
    inherit src;
    buildInputs = [ pkgs.cmake pkgs.boost tklog tkassert tkrng ];
    postFixup = ''
      # fix bogus include paths
      # trick found here: https://github.com/NixOS/nixpkgs/blob/master/pkgs/development/libraries/crc32c/default.nix
      for f in $(find $out/lib/cmake -name '*.cmake'); do
        substituteInPlace "$f" --replace "\''${_IMPORT_PREFIX}/$out/include" "\''${_IMPORT_PREFIX}/include"
      done
    '';
  }
