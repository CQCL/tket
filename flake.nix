{
  description = "Tket Quantum SDK";
  inputs.nixpkgs.url = github:nixos/nixpkgs/nixos-23.05;
  inputs.flake-utils.url = github:numtide/flake-utils;
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
      in
        {
          defaultPackage = pkgs.tket-static;
          packages = {
            tket = pkgs.tket;
            tket-static = pkgs.tket-static;
            tket-shared = pkgs.tket-shared;
            pytket = pkgs.pytket;
          };
          devShell = with pkgs; mkShell {
            buildInputs = [
              pkgs.tket-static
              pkgs.pytket
            ];
          };
        }
    );
}
