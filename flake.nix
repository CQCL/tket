{
  description = "Tket Quantum SDK";
  inputs.nixpkgs.url = github:nixos/nixpkgs/nixos-23.05;
  inputs.flake-utils.url = github:numtide/flake-utils;
  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };
        libs = pkgs.callPackage ./libs { };
        tket = pkgs.callPackage ./tket { inherit libs; };
        pytket = pkgs.callPackage ./pytket { inherit tket; };
      in
        {
          defaultPackage = tket;
          packages = {
            tket = tket;
            pytket = pytket;
          };
          devShell = with pkgs; mkShell {
            buildInputs = [
              tket
              pytket
            ];
          };
        }
    );
}
