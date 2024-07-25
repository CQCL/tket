self: super:
let
  config_contents = builtins.readFile ../pytket/docs/conf.py;
  versions =
    builtins.match ''.*release *= *["']([^"']+)["'].*'' config_contents;
  version = if builtins.length versions > 0 then
    builtins.elemAt versions 0
  else
    builtins.trace "Warning: Unable to find version. Defaulting to 0.0.0" "0.0.0";

  jsonschema-4180 = super.python3Packages.jsonschema.overrideAttrs (_: rec {
    version = "4.18.0";
    src = super.fetchPypi {
      pname = "jsonschema";
      version = "4.18.0";
      hash = sha256:jK9bV6mQqY6bOYMu88s1wXb+MxQUJStuGyb9WGb4kaQ=;
    };
  });
in {
  pytket = super.python3.pkgs.buildPythonPackage {
    pname = "pytket";
    inherit version;
    src = ../pytket;
    nativeBuildInputs = with super.python3.pkgs; [
      setuptools
      super.cmake
      super.pkg-config
    ];
    propagatedBuildInputs = [
      super.tket
      super.pybind11_json
    ];
    dependencies = with super.python3.pkgs; [
      super.lark
      super.qwasm
      pybind11
      graphviz
      networkx
      jinja2
      super.sympy'
      scipy
      numpy
      typing-extensions
    ];
    configurePhase = "true"; # skip configure
    #                          (by default it runs cmake, which we
    #                           want python to manage instead)
    preBuild = ''
      # explicitly provide the pytket version to setup.py
      # and to _version.py, as we can't rely on git version
      # information in nix.

      cat ${../pytket/setup.py} \
        | sed 's/setup(/setup(version="${version}",/' \
        > setup.py;
      echo '__version__ = "${version}"' > pytket/_version.py;
      
      # instruct python to build with cmake instead of conan,
      # and to build in a temporary directory.
      export NO_CONAN=1;
      export INSTALL_DIR=$(mktemp -d);
      export BUILD_DIR=$(mktemp -d);
    '';
    checkInputs = with super.python3.pkgs; [
      mypy
      pytest
      pytest-cov
      pytest-benchmark
      py
      hypothesis
      docker
      opt-einsum
    ] ++ [jsonschema-4180];
    checkPhase = ''
      export HOME=$TMPDIR;

      # run mypy
      python -m mypy --config-file=mypy.ini --no-incremental -p pytket -p test_root.tests;

      # run tests
      chmod 700 $TMPDIR/test_root/tests/qasm_test_files;
      cd test_root/tests;
      python -m pytest -s .
    '';
    doCheck = false; #TODO revert to true
  };
}
