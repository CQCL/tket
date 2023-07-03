{
  pkgs ? import <nixpkgs> {},
  python_version ? pkgs.python3,
  libs ? import ../libs { inherit pkgs; },
  tket ? import ../tket { inherit pkgs; inherit libs; }
}:
let
  pybind11_json = pkgs.stdenv.mkDerivation{
    name = "pybind11_json";
    src = pkgs.fetchFromGitHub {
      owner = "pybind";
      repo = "pybind11_json";
      rev = "0.2.13";
      sha256 = sha256:Kl/QflV2bBoH72/LW03K8JDlhBF+DYYXL47A5s1nmTw=;
    };
    nativeBuildInputs = [ pkgs.cmake ];
    buildInputs = [
      python_version.pkgs.pybind11
      pkgs.nlohmann_json
    ];
  };
  binders = pkgs.stdenv.mkDerivation{
    name = "binders";
    nativeBuildInputs = [ pkgs.cmake pkgs.pkgconfig ];
    propagatedBuildInputs = [
      python_version.pkgs.pybind11
      pybind11_json
      pkgs.boost
      pkgs.eigen
      tket
    ] ++ tket.buildInputs;
    unpackPhase = ''
      cp -r ${./binders} binders;
      cp ${./CMakeLists.txt} CMakeLists.txt;
    '';
  };
  qwasm = python_version.pkgs.buildPythonPackage{
    name = "qwasm";
    src = pkgs.fetchFromGitHub {
      owner = "CQCL";
      repo = "qwasm";
      rev = "35ebe1e2551449d97b9948a600f8d2e4d7474df6";
      sha256 = sha256:g/QA5CpAR3exRDgVQMnXGIH8bEGtwGFBjjSblbdXRkU=;
    };
  };
  lark-parser = python_version.pkgs.buildPythonPackage {
      pname = "lark-parser";
      version = "0.7.8";
      src = pkgs.fetchFromGitHub{
        owner = "lark-parser";
        repo = "lark";
        rev = "refs/tags/0.7.8";
        hash = sha256:XwEcQO8u9UnLCKjnA3eAdTbktCyLQ8oQHDn6y/RgpT0=;
      };
      doCheck = false;
    };
in
python_version.pkgs.buildPythonPackage{
  name = "pytket";
  nativeBuildInputs = [ ];
  buildInputs = [
    binders
    qwasm
  ] ++ (with python_version.pkgs; [
    graphviz
    networkx
    jinja2
    sympy
    scipy
    numpy
    typing-extensions
    lark-parser
  ]);
  
  unpackPhase = ''
    cp -r ${./pytket} pytket;
    cp -r ${./setup.py} setup.py;
    cp -r ${./package.md} package.md;
  '';
  preBuild = ''
    export USE_NIX=1;
  '';
}

