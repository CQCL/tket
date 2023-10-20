self: super:
let
  default-flags = [ "-DBUILD_SHARED_LIBS=ON" "-DINSTALL_NAME_DIR=OFF" ];
  postFixup = import ./includes-fixup.nix;
in {
  tklog = super.stdenv.mkDerivation {
    name = "tklog";
    src = ../libs/tklog;
    nativeBuildInputs = [ super.cmake ];
    cmakeFlags = default-flags;
    inherit postFixup;
  };
  tkrng = super.stdenv.mkDerivation {
    name = "tkrng";
    src = ../libs/tkrng;
    nativeBuildInputs = [ super.cmake ];
    cmakeFlags = default-flags;
    inherit postFixup;
  };
  tkassert = super.stdenv.mkDerivation {
    name = "tkassert";
    src = ../libs/tkassert;
    nativeBuildInputs = [ super.cmake ];
    buildInputs = [ self.tklog ];
    cmakeFlags = default-flags;
    inherit postFixup;
  };
  tktokenswap = super.stdenv.mkDerivation {
    name = "tktokenswap";
    src = ../libs/tktokenswap;
    nativeBuildInputs = [ super.cmake super.boost ];
    buildInputs = [ self.tklog self.tkassert self.tkrng ];
    cmakeFlags = default-flags;
    inherit postFixup;
  };
  tkwsm = super.stdenv.mkDerivation {
    name = "tkwsm";
    src = ../libs/tkwsm;
    nativeBuildInputs = [ super.cmake super.boost ];
    buildInputs = [ self.tklog self.tkassert self.tkrng ];
    cmakeFlags = default-flags;
    inherit postFixup;
  };
  tklibs = [ self.tklog self.tkrng self.tkassert self.tktokenswap self.tkwsm ];
}
