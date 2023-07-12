self: super:
let
  shared-flags = [ "-DBUILD_SHARED_LIBS=ON" ];
  postFixup = import ./includes-fixup.nix;
in {
  tklog = super.stdenv.mkDerivation {
    name = "tklog";
    src = ../libs/tklog;
    nativeBuildInputs = [ super.cmake ];
    cmakeFlags = shared-flags;
    inherit postFixup;
  };
  tkrng = super.stdenv.mkDerivation {
    name = "tkrng";
    src = ../libs/tkrng;
    nativeBuildInputs = [ super.cmake ];
    cmakeFlags = shared-flags;
    inherit postFixup;
  };
  tkassert = super.stdenv.mkDerivation {
    name = "tkassert";
    src = ../libs/tkassert;
    nativeBuildInputs = [ super.cmake ];
    buildInputs = [ self.tklog ];
    cmakeFlags = shared-flags;
    inherit postFixup;
  };
  tktokenswap = super.stdenv.mkDerivation {
    name = "tktokenswap";
    src = ../libs/tktokenswap;
    nativeBuildInputs = [ super.cmake super.boost ];
    buildInputs = [ self.tklog self.tkassert self.tkrng ];
    cmakeFlags = shared-flags;
    inherit postFixup;
  };
  tkwsm = super.stdenv.mkDerivation {
    name = "tkwsm";
    src = ../libs/tkwsm;
    nativeBuildInputs = [ super.cmake super.boost ];
    buildInputs = [ self.tklog self.tkassert self.tkrng ];
    cmakeFlags = shared-flags;
    inherit postFixup;
  };
  tklibs = [ self.tklog self.tkrng self.tkassert self.tktokenswap self.tkwsm ];
}
