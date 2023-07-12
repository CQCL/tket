# Example flake project

In this example, we create a [flake.nix](flake.nix) file that tells Nix
where to fetch nixpkgs, flake-utils and tket. It also builds an example
application, using [pyproject.toml](pyproject.toml) and the
[src](src) directory. This application [exposes](src/examples/basic_circuits)
two scripts: `entanglement` and `teleport`, that set up trivial circuits
and display them in your browser.

To run these scripts, simply run:

```
nix run .#entanglement
```
and
```
nix run .#teleport
```
