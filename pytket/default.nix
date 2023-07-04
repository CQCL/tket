{
  pkgs ? import <nixpkgs> {},
  python_version ? pkgs.python3,
  libs ? import ../libs { inherit pkgs; },
  tket ? import ../tket { inherit pkgs; inherit libs; },
  third-party ? import ./third-party.nix { inherit pkgs; inherit python_version; }
}:
let
  binders = pkgs.stdenv.mkDerivation{
    name = "binders";
    nativeBuildInputs = [ pkgs.cmake pkgs.pkgconfig ];
    propagatedBuildInputs = [
      python_version.pkgs.pybind11
      third-party.pybind11_json
      pkgs.boost
      pkgs.eigen
      tket
    ] ++ tket.buildInputs;
    unpackPhase = ''
      cp -r ${./binders} binders;
      cp ${./CMakeLists.txt} CMakeLists.txt;
    '';
  };
  pytket = python_version.pkgs.buildPythonPackage{
    name = "pytket";
    nativeBuildInputs = [ ];
    buildInputs = [
      binders
      third-party.qwasm
    ] ++ (with python_version.pkgs; [
      graphviz
      networkx
      jinja2
      sympy
      scipy
      numpy
      typing-extensions
      third-party.lark-parser
      third-party.types-pkg_resources
    ]);
    
    unpackPhase = ''
      cp -r ${./pytket} pytket;
      cp -r ${./setup.py} setup.py;
      cp -r ${./package.md} package.md;
    '';
    preBuild = ''
      export USE_NIX=1;
    '';
    doCheck = false;
  };
  pytket-tests = 
in {
  pytket = pytket;
  pytket-tests = pytket-tests;
}

