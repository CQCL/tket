self: super:
let
  static-shared-flags = static: if static then [ "-DBUILD_SHARED_LIBS=OFF" ] else [ "-DBUILD_SHARED_LIBS=ON" ];
  postFixup = import ./includes-fixup.nix;

  tklog-base = {build-static}: super.stdenv.mkDerivation{
      name = "tklog";
      src = ../libs/tklog;
      nativeBuildInputs = [ super.cmake ];
      cmakeFlags = static-shared-flags build-static;
      inherit postFixup;
    };
  tkrng-base = {build-static}: super.stdenv.mkDerivation{
      name = "tkrng";
      src = ../libs/tkrng;
      nativeBuildInputs = [ super.cmake ];
      cmakeFlags = static-shared-flags build-static;
      inherit postFixup;
    };
  tkassert-base = {build-static, tklog}: super.stdenv.mkDerivation{
      name = "tkassert";
      src = ../libs/tkassert;
      nativeBuildInputs = [ super.cmake ];
      buildInputs = [ tklog ];
      cmakeFlags = static-shared-flags build-static;
      inherit postFixup;
    };
  tktokenswap-base = {build-static, tklog, tkassert, tkrng}:
    super.stdenv.mkDerivation{
      name = "tktokenswap";
      src = ../libs/tktokenswap;
      nativeBuildInputs = [ super.cmake super.boost ];
      buildInputs = [ tklog tkassert tkrng ];
      cmakeFlags = static-shared-flags build-static;
      inherit postFixup;
    };
  tkwsm-base = {build-static, tklog, tkassert, tkrng}:
    super.stdenv.mkDerivation {
      name = "tkrng";
      src = ../libs/tkwsm;
      nativeBuildInputs = [ super.cmake super.boost ];
      buildInputs = [ tklog tkassert tkrng ];
      cmakeFlags = static-shared-flags build-static;
      inherit postFixup;
    };
in
{
  tklog-static = tklog-base {build-static = true;};
  tklog-shared = tklog-base {build-static = false;};
  tklog = self.tklog-static;

  tkrng-static = tkrng-base {build-static = true;};
  tkrng-shared = tkrng-base {build-static = false;};
  tkrng = self.tkrng-static;

  tkassert-static = tkassert-base {build-static = true; tklog = self.tklog-static;};
  tkassert-shared = tkassert-base {build-static = false; tklog = self.tklog-shared;};
  tkassert = self.tkassert-static;

  tktokenswap-static = tktokenswap-base {build-static = true; tklog = self.tklog-static; tkassert = self.tkassert-static; tkrng = self.tkrng-static;};
  tktokenswap-shared = tktokenswap-base {build-static = false; tklog = self.tklog-shared; tkassert = self.tkassert-shared; tkrng = self.tkrng-shared;};
  tktokenswap = self.tktokenswap-static;

  tkwsm-static = tkwsm-base {build-static = true; tklog = self.tklog-static; tkassert = self.tkassert-static; tkrng = self.tkrng-static;};
  tkwsm-shared = tkwsm-base {build-static = false; tklog = self.tklog-shared; tkassert = self.tkassert-shared; tkrng = self.tkrng-shared;};
  tkwsm = self.tkwsm-static;
  
  tklibs-static = [
    self.tklog-static
    self.tkrng-static
    self.tkassert-static
    self.tktokenswap-static
    self.tkwsm-static
  ];

  tklibs-shared = [
    self.tklog-shared
    self.tkrng-shared
    self.tkassert-shared
    self.tktokenswap-shared
    self.tkwsm-shared
  ];
}
