self: super:
let
  default-flags = [ "-DBUILD_SHARED_LIBS=ON" "-DINSTALL_NAME_DIR=OFF" ];
in {
  tklog = super.stdenv.mkDerivation {
    name = "tklog";
    src = ../libs/tklog;
    nativeBuildInputs = [ super.cmake ];
    cmakeFlags = default-flags;
  };
  tkrng = super.stdenv.mkDerivation {
    name = "tkrng";
    src = ../libs/tkrng;
    nativeBuildInputs = [ super.cmake ];
    cmakeFlags = default-flags;
  };
  tkassert = super.stdenv.mkDerivation {
    name = "tkassert";
    src = ../libs/tkassert;
    nativeBuildInputs = [ super.cmake ];
    buildInputs = [ self.tklog ];
    cmakeFlags = default-flags;
  };
  tktokenswap = super.stdenv.mkDerivation {
    name = "tktokenswap";
    src = ../libs/tktokenswap;
    nativeBuildInputs = [ super.cmake super.boost ];
    buildInputs = [ self.tklog self.tkassert self.tkrng ];
    cmakeFlags = default-flags;
  };
  tkwsm = super.stdenv.mkDerivation {
    name = "tkwsm";
    src = ../libs/tkwsm;
    nativeBuildInputs = [ super.cmake super.boost ];
    buildInputs = [ self.tklog self.tkassert self.tkrng ];
    cmakeFlags = default-flags;
  };
  tklibs = [ self.tklog self.tkrng self.tkassert self.tktokenswap self.tkwsm ];
}
