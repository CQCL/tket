self: super:
let
  inputs-base = [ 
      super.cmake
      super.boost
      super.eigen
      super.nlohmann_json
  ];
  postFixup = import ./includes-fixup.nix;
in
{
  tket-static = super.stdenv.mkDerivation{
    name = "tket";
    src = ../tket;
    buildInputs = inputs-base ++ super.tklibs-static ++ [super.symengine-static];
    cmakeFlags = ["-DBUILD_SHARED_LIBS=OFF"];
    inherit postFixup;
  };
  tket-shared = super.stdenv.mkDerivation{
    name = "tket";
    src = ../tket;
    buildInputs = inputs-base ++ super.tklibs-shared ++ [super.symengine-shared];
    cmakeFlags = ["-DBUILD_SHARED_LIBS=ON"];
    inherit postFixup;
  };
  tket = self.tket-static;
}
