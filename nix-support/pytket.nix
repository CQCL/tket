self: super:
let
  config_contents = builtins.readFile ../pytket/docs/conf.py;
  versions =
    builtins.match ''.*release *= *["']([^"']+)["'].*'' config_contents;
  version = if builtins.length versions > 0 then
    builtins.elemAt versions 0
  else
    "0.0.0";

  jsonschema-4180 = super.python3Packages.jsonschema.overrideAttrs (_: rec {
    version = "4.18.0";
    src = super.fetchPypi {
      pname = "jsonschema";
      version = "4.18.0";
      hash = sha256:jK9bV6mQqY6bOYMu88s1wXb+MxQUJStuGyb9WGb4kaQ=;
    };
  });
in {
  binders = super.stdenv.mkDerivation {
    name = "binders";
    nativeBuildInputs = [
      super.cmake
      super.pkg-config
      super.python3Packages.pybind11
      super.pybind11_json
    ];
    cmakeFlags = [ "-DBUILD_SHARED_LIBS=ON" ];
    propagatedBuildInputs = [
      super.tket
    ];
    unpackPhase = ''
      cp -r ${../pytket/binders} binders;
      cp ${../pytket/CMakeLists.txt} CMakeLists.txt;
    '';
  };
  pytket = super.python3.pkgs.buildPythonPackage {
    name = "pytket";
    inherit version;
    propagatedBuildInputs = with super.python3.pkgs; [
      self.binders
      super.lark-parser
      super.types-pkg_resources
      super.qwasm
      graphviz
      networkx
      jinja2
      sympy
      scipy
      numpy
      typing-extensions
    ];

    unpackPhase = ''
      cp -r ${../pytket/pytket} pytket;
      cp -r ${../pytket/setup.py} setup.py;
      cp -r ${../pytket/package.md} package.md;
      cp -r ${../schemas} schemas;
      mkdir test_root;
      cp -r ${../pytket/tests} test_root/tests;
    '';
    preBuild = ''
      export USE_NIX=1;
    '';
    postFixup = ''
      # hardcode the version extracted from docs/conf.py.
      echo '__version__ = "${version}"' > $out/lib/python3.10/site-packages/pytket/_version.py;

      # these directories aren't copied by setup.py, so we do it manually
      cp -r ${
        ../pytket/pytket/circuit/display/js
      } $out/lib/python3.10/site-packages/pytket/circuit/display/js;
      cp -r ${
        ../pytket/pytket/circuit/display/static
      } $out/lib/python3.10/site-packages/pytket/circuit/display/static;
    '';
    checkInputs = with super.python3.pkgs; [
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
      chmod 700 $TMPDIR/test_root/tests/qasm_test_files;
      cd test_root/tests;
      python -m pytest -s .
    '';
    doCheck = true;
  };
}
