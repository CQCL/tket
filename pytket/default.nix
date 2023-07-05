{
  pkgs ? import <nixpkgs> {},
  python_version ? pkgs.python3,
  libs ? import ../libs { inherit pkgs; static=false; },
  tket ? import ../tket { inherit pkgs; static=false; inherit libs; },
  third-party ? import ./../nix-helpers/third-party-python-packages.nix { inherit pkgs python_version; },
  symengine ? import ../../nix-helpers/symengine.nix { inherit pkgs; static=false; }
}:
let
  binders = pkgs.stdenv.mkDerivation{
    name = "binders";
    nativeBuildInputs = [
      pkgs.cmake
      pkgs.pkgconfig
      python_version.pkgs.pybind11
      third-party.pybind11_json
    ];
    cmakeFlags = [
      "-DBUILD_SHARED_LIBS=ON"
    ];
    buildInputs = [
      pkgs.boost
      pkgs.eigen
      symengine
      pkgs.nlohmann_json
      tket
      libs.tklog
      libs.tkassert
      libs.tkrng
      libs.tktokenswap
      libs.tkwsm
    ] ++ symengine.buildInputs;
    unpackPhase = ''
      cp -r ${./binders} binders;
      cp ${./CMakeLists.txt} CMakeLists.txt;
    '';
  };
  pytket = python_version.pkgs.buildPythonPackage{
    name = "pytket";
    propagatedBuildInputs = with python_version.pkgs; [
      binders
      graphviz
      networkx
      jinja2
      sympy
      scipy
      numpy
      typing-extensions
      third-party.lark-parser
      third-party.types-pkg_resources
      third-party.qwasm
    ] ++ binders.buildInputs;
    
    unpackPhase = ''
      cp -r ${./pytket} pytket;
      cp -r ${./setup.py} setup.py;
      cp -r ${./package.md} package.md;
    '';
    preBuild = ''
      export USE_NIX=1;
    '';
    postFixup = ''
      tketlib=$out/lib/python3.10/site-packages/pytket/_tket/libtket.so;
      patchelf --add-needed "libsymengine.so" $tketlib;
      patchelf --set-rpath "''${ORIGIN}:${rpath}" $tketlib;
      echo '__version__ = "0.0.0"' > $out/lib/python3.10/site-packages/pytket/_version.py;
    '';
    doCheck = false;
  };

  rpath = builtins.concatStringsSep ":" (
    map (x: "${x}/lib")
    [libs.tklog libs.tkassert libs.tkrng libs.tktokenswap libs.tkwsm
    pkgs.gcc.libc pkgs.gcc.cc.lib.lib symengine]
  );
in
  pytket

