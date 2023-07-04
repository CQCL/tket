{
  pkgs ? import <nixpkgs> {},
  python_version ? pkgs.python3
}:
{
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
  types-pkg_resources = let
    pname = "types-pkg_resources";
    version = "0.1.3"; 
  in python_version.pkgs.buildPythonPackage {
    inherit pname version;
    src = pkgs.fetchPypi {
      inherit pname version;
      sha256 = sha256:g0qbjT2+o0NWL9mdXTNZpyb2v503M7zNK08wlvurna4=;
    };
    doCheck = false;
  };
}
