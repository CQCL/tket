{
  inputs.nixpkgs.url = "github:nixos/nixpkgs/nixos-23.05";
  inputs.flake-utils.url = "github:numtide/flake-utils";
  # fetch tket from github
  inputs.tket.url = "github:CQCL/tket";

  outputs = { self, nixpkgs, flake-utils, tket }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          overlays = [
            (self: super: {
              # add tket and pytket to pkgs for later use
              tket = tket.packages.${system}.tket;
              pytket = tket.packages.${system}.pytket;
            })
          ];
        };
        examples = pkgs.python3Packages.buildPythonApplication {
          pname = "examples";
          version = "0.0.1";
          format = "pyproject";
          # copy pyproject.toml and src to the build directory
          unpackPhase = ''
            cp ${./pyproject.toml} pyproject.toml
            cp -r ${./src} src
            chmod 700 src;
          '';
          # provide pytket as a dependency
          propagatedBuildInputs = [ pkgs.pytket ];
        };
      in {
        apps = {
          entanglement = {
            type = "app";
            program = "${examples}/bin/entanglement";
          };
          teleport = {
            type = "app";
            program = "${examples}/bin/teleport";
          };
        };
      });
}
