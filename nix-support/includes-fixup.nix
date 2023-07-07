''
    # fix bogus include paths
    # trick found here: https://github.com/NixOS/nixpkgs/blob/master/pkgs/development/libraries/crc32c/default.nix
    for f in $(find $out/lib/cmake -name '*.cmake'); do
      substituteInPlace "$f" --replace "\''${_IMPORT_PREFIX}/$out/include" "\''${_IMPORT_PREFIX}/include"
    done
''
