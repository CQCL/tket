self: super:
let
  postFixup = import ./includes-fixup.nix;
  # only import necessary directories so that changes to unrelated
  # files inside of ../tket won't trigger a rebuild
  tket-without-tests = super.stdenv.mkDerivation {
    name = "tket-sources";
    phases = [ "installPhase" ];
    installPhase = ''
      mkdir -p $out;
      cp -r ${../tket/cmake} $out/cmake;
      cp -r ${../tket/include} $out/include;
      cp -r ${../tket/src} $out/src;
      cp -r ${../tket/CMakeLists.txt} $out/CMakeLists.txt;
    '';
  };
in {
  tket = super.stdenv.mkDerivation {
    name = "tket";
    src = tket-without-tests;
    nativeBuildInputs = [ super.cmake ];
    propagatedBuildInputs = super.tklibs
      ++ [ super.boost super.symengine super.eigen super.nlohmann_json ];
    cmakeFlags = [ "-DBUILD_SHARED_LIBS=ON" ];
    inherit postFixup;
  };

  tket-tests = super.stdenv.mkDerivation {
    name = "tket-tests";
    src = ../tket/test;
    nativeBuildInputs = [ super.cmake super.pkg-config ];
    buildInputs = [ self.tket super.catch2_3 ];
  };
  run-tket-tests = super.stdenv.mkDerivation {
    name = "run-tket-tests";
    stages = [ "build" ];
    buildCommand = ''
      pushd ${self.tket-tests}/bin;
      mkdir -p $out;
      ./test-tket > $out/test_result.txt;
      popd;
    '';
  };
}
