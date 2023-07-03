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
            (self: super: {
              tket = import ./tket.nix { pkgs=super; };
            })
          ];
        };
      in
        {
          defaultPackage = pkgs.tket;
        }
    );
}
