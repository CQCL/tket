self: super:
let
  # only import necessary directories so that changes to unrelated
  # files inside of ../tket won't trigger a rebuild
  src = super.stdenv.mkDerivation {
    name = "tket-sources";
    phases = [ "installPhase" ];
    installPhase = ''
      mkdir -p $out;
      cp -r ${../tket/cmake} $out/cmake;
      cp -r ${../tket/include} $out/include;
      cp -r ${../tket/src} $out/src;
      cp -r ${../tket/test} $out/test;
      cp -r ${../tket/proptest} $out/proptest;
      cp ${../tket/CMakeLists.txt} $out/CMakeLists.txt;
    '';
  };
in {
  tket = super.stdenv.mkDerivation {
    name = "tket";
    inherit src;
    nativeBuildInputs = [ super.cmake ];
    propagatedBuildInputs = super.tklibs
      ++ (with super; [ boost symengine eigen nlohmann_json ]);
    # TODO: add rapidcheck once nixpkgs packaging is correctly implemented.
    # At current the package fails because the .pc files incorrectly reference
    # the library dir rather than the dev dir.
    # See https://github.com/NixOS/nixpkgs/issues/296348
    buildInputs = with super; [ catch2_3 ];
    cmakeFlags = [
      "-DBUILD_SHARED_LIBS=ON"
      "-DINSTALL_NAME_DIR=OFF"
      "-DBUILD_TKET_TEST=ON"
      "-DBUILD_TKET_PROPTEST=OFF"
    ];
    doCheck = true;
  };
}
