{
  description = "Tket Quantum SDK";
  nixConfig.extra-substituters = "https://tket.cachix.org https://cache.nixos.org";
  nixConfig.trusted-public-keys = ''
    tket.cachix.org-1:ACdm5Zg19qPL0PpvUwTPPiIx8SEUy+D/uqa9vKJFwh0=
    cache.nixos.org-1:6NCHdD59X431o0gWypbMrAURkbJ16ZPMQFGspcDShjY=
  '';
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
          conan = pkgs.mkShell { buildInputs = with pkgs; [ conan cmake ninja gcc python3 ]; };
        };
        checks = {
          tket-tests = pkgs.tket;
          pytket-tests = pkgs.pytket;
        };
      });
}
