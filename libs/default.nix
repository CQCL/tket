{
  pkgs ? import <nixpkgs> {},
  static ? true
}:
let
  tklog = import ./tklog { inherit pkgs static; };
  tkrng = import ./tkrng { inherit pkgs static; };
  tkassert = import ./tkassert { inherit pkgs static tklog; };
  tktokenswap = import ./tktokenswap { inherit pkgs static tklog tkassert tkrng; };
  tkwsm = import ./tkwsm { inherit pkgs static tklog tkassert tkrng; };
in
  {
    inherit tklog;
    inherit tkassert;
    inherit tkrng;
    inherit tktokenswap;
    inherit tkwsm;
  }
