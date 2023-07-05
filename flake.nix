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

        # static libs
        libs-static = pkgs.callPackage ./libs { static=true; };
        tket-static = pkgs.callPackage ./tket { libs=libs-static; static=true; };
        symengine-static = pkgs.callPackage ./nix-helpers/symengine.nix { static=true; };
        # shared libs
        libs-shared = pkgs.callPackage ./libs { static=false; };
        symengine-shared = pkgs.callPackage ./nix-helpers/symengine.nix { static=false; };
        tket-shared = pkgs.callPackage ./tket { libs=libs-shared; static=false; symengine=symengine-shared; };
        # pytket
        pytket = pkgs.callPackage ./pytket {
          libs=libs-shared;
          tket=tket-shared;
          symengine=symengine-shared;
        };
      in
        {
          defaultPackage = tket-static;
          packages = {
            tket = tket-static;
            tket-static = tket-static;
            tket-shared = tket-shared;
            pytket = pytket;
            symengine-static = symengine-static;
            symengine-shared = symengine-shared;
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
