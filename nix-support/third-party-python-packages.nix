self: super: {
  pybind11_json = super.stdenv.mkDerivation {
    name = "pybind11_json";
    src = super.fetchFromGitHub {
      owner = "pybind";
      repo = "pybind11_json";
      rev = "0.2.14";
      sha256 = "sha256-6L675DsfafzRv0mRR3b0eUFFjUpll3jCPoBAAffk7U0=";
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
  sympy' = super.python3.pkgs.buildPythonPackage rec{
    # version bump - nixpkgs' version is 1.12 at the time of writing
    pname = "sympy";
    version = "1.13.0";
    format = "setuptools";
    src = super.python3Packages.fetchPypi {
      inherit pname version;
      sha256 = "sha256:O2r49NAIuaGmpCaLM1uYSyODXybR1gsFJuvHHUiiX1c=";
    };
    nativeCheckInputs = [ super.glibcLocales ];
    propagatedBuildInputs = [ super.python3Packages.mpmath ];
    # tests take ~1h
    doCheck = false;
    pythonImportsCheck = [ "sympy" ];
    preCheck = ''
      export LANG="en_US.UTF-8"
    '';
  };
}
