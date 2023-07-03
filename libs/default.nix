{ pkgs ? import <nixpkgs> {} }:
let
  tklog = import ./tklog { inherit pkgs; };
  tkrng = import ./tkrng { inherit pkgs; };
  tkassert = import ./tkassert { inherit pkgs; inherit tklog; };
  tktokenswap = import ./tktokenswap { inherit pkgs; inherit tklog; inherit tkassert; inherit tkrng; };
  tkwsm = import ./tkwsm { inherit pkgs; inherit tklog; inherit tkassert; inherit tkrng; };
in
  {
    inherit tklog;
    inherit tkassert;
    inherit tkrng;
    inherit tktokenswap;
    inherit tkwsm;
  }
