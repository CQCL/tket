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
        libs-static = pkgs.callPackage ./libs { static=true; };
        libs-shared = pkgs.callPackage ./libs { static=false; };
        tket-static = pkgs.callPackage ./tket { libs=libs-static; static=true; };
        tket-shared = pkgs.callPackage ./tket { libs=libs-shared; static=false; };
        pytket = pkgs.callPackage ./pytket { tket=tket-shared; };
      in
        {
          defaultPackage = tket-static;
          packages = {
            tket = tket-static;
            tket-shared = tket-shared;
            tket-static = tket-static;
            pytket = pytket;
          };
          devShell = with pkgs; mkShell {
            buildInputs = [
              tket-static
              pytket
            ];
          };
        }
    );
}
