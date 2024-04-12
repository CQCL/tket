{
  description = "Tket Quantum SDK";
  nixConfig.extra-substituters = "https://tket.cachix.org";
  nixConfig.trusted-public-keys = "tket.cachix.org-1:ACdm5Zg19qPL0PpvUwTPPiIx8SEUy+D/uqa9vKJFwh0=";
  inputs.nixpkgs.url = "github:nixos/nixpkgs";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [
            (import ./nix-support/libs.nix)
            (import ./nix-support/symengine.nix)
            (import ./nix-support/tket.nix)
            (import ./nix-support/third-party-python-packages.nix)
            (import ./nix-support/pytket.nix)
          ];
        };
      in {
        packages = {
          tket = pkgs.tket;
          pytket = pkgs.pytket;
        };
        devShells = {
          default = pkgs.mkShell { buildInputs = [ pkgs.tket pkgs.pytket ]; };
        };
        checks = {
          tket-tests = pkgs.run-tket-tests;
          pytket-tests = pkgs.pytket;
        };
      });
}
