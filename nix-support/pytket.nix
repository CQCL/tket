self: super: 
let
  config_contents = builtins.readFile ../pytket/docs/conf.py;
  versions = builtins.match ''.*''\nrelease *= *["']([^"']+)["'].*'' config_contents;
  version = if builtins.length versions > 0 then builtins.elemAt versions 0 else "0.0.0";
in
{
  binders = super.stdenv.mkDerivation{
    name = "binders";
    nativeBuildInputs = [
      super.cmake
      super.pkgconfig
      super.python3.pkgs.pybind11
      super.pybind11_json
    ];
    cmakeFlags = ["-DBUILD_SHARED_LIBS=ON" ];
    buildInputs = [
      super.boost
      super.eigen
      super.symengine-shared
      super.nlohmann_json
      super.tket-shared
    ] ++ super.tklibs-shared ++ super.symengine-shared.buildInputs;
    unpackPhase = ''
      cp -r ${../pytket/binders} binders;
      cp ${../pytket/CMakeLists.txt} CMakeLists.txt;
    '';
  };
  pytket = super.python3.pkgs.buildPythonPackage{
    name = "pytket";
    inherit version;
    propagatedBuildInputs = with super.python3.pkgs; [
      self.binders
      graphviz
      networkx
      jinja2
      sympy
      scipy
      numpy
      typing-extensions
      super.lark-parser
      super.types-pkg_resources
      super.qwasm
    ] ++ self.binders.buildInputs;
    
    unpackPhase = ''
      cp -r ${../pytket/pytket} pytket;
      cp -r ${../pytket/setup.py} setup.py;
      cp -r ${../pytket/package.md} package.md;
    '';
    preBuild = ''
      export USE_NIX=1;
    '';
    postFixup = ''
      tketlib=$out/lib/python3.10/site-packages/pytket/_tket/libtket.so;
      echo '__version__ = "${version}"' > $out/lib/python3.10/site-packages/pytket/_version.py;
    '';
    doCheck = false;
  };
}
