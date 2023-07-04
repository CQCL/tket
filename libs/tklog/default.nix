{ 
  pkgs ? import <nixpkgs> {},
  static ? true
}:
let
  src = builtins.filterSource(p: _: baseNameOf p != "default.nix") ./.;
in
  pkgs.stdenv.mkDerivation{
    name = "tklog";
    inherit src;
    buildInputs = [ pkgs.cmake ];
    cmakeFlags = if static then [ "-DBUILD_SHARED_LIBS=OFF" ] else [ "-DBUILD_SHARED_LIBS=ON" ];
    postFixup = ''
      # fix bogus include paths
      # trick found here: https://github.com/NixOS/nixpkgs/blob/master/pkgs/development/libraries/crc32c/default.nix
      # this is needed because when Nix builds the cmake project, it doubly prefixes the include paths:
      # e.g.:
      #  set(_IMPORT_PREFIX "/nix/store/wbvzlgvc7984xmpjd6lf3cg2490xb81b-tklog")
      #  ...
      #  target_sources(tklog::tklog
      #    INTERFACE
      #      FILE_SET "HEADERS"
      #      TYPE "HEADERS"
      #      BASE_DIRS "\''${_IMPORT_PREFIX}//nix/store/wbvzlgvc7984xmpjd6lf3cg2490xb81b-tklog/include"
      #      FILES "\''${_IMPORT_PREFIX}//nix/store/wbvzlgvc7984xmpjd6lf3cg2490xb81b-tklog/include/tklog/TketLog.hpp"
      #  )
      for f in $(find $out/lib/cmake -name '*.cmake'); do
        substituteInPlace "$f" --replace "\''${_IMPORT_PREFIX}/$out/include" "\''${_IMPORT_PREFIX}/include"
      done
    '';
  }
