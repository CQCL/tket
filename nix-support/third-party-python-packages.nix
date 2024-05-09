self: super: {
  pybind11_json = super.stdenv.mkDerivation {
    name = "pybind11_json";
    src = super.fetchFromGitHub {
      owner = "pybind";
      repo = "pybind11_json";
      rev = "0.2.13";
      sha256 = "sha256:Kl/QflV2bBoH72/LW03K8JDlhBF+DYYXL47A5s1nmTw=";
    };
    nativeBuildInputs = [ super.cmake ];
    buildInputs = [ super.python3Packages.pybind11 super.nlohmann_json ];
  };
  qwasm = super.python3.pkgs.buildPythonPackage {
    name = "qwasm";
    src = super.fetchFromGitHub {
      owner = "CQCL";
      repo = "qwasm";
      rev = "35ebe1e2551449d97b9948a600f8d2e4d7474df6";
      sha256 = "sha256:g/QA5CpAR3exRDgVQMnXGIH8bEGtwGFBjjSblbdXRkU=";
    };
  };
  lark = super.python3.pkgs.buildPythonPackage {
    pname = "lark";
    version = "1.1.9";
    format = "pyproject";
    src = super.fetchFromGitHub {
      owner = "lark-parser";
      repo = "lark";
      rev = "refs/tags/1.1.9";
      hash = "sha256:pWLKjELy10VNumpBHjBYCO2TltKsZx1GhQcGMHsYJNk=";
    };
    nativeBuildInputs = with super.python3Packages; [ setuptools ];
    doCheck = false;
  };
  types-pkg_resources = let
    pname = "types-pkg_resources";
    version = "0.1.3";
  in super.python3.pkgs.buildPythonPackage {
    inherit pname version;
    src = super.fetchPypi {
      inherit pname version;
      sha256 = "sha256:g0qbjT2+o0NWL9mdXTNZpyb2v503M7zNK08wlvurna4=";
    };
    doCheck = false;
  };
}
