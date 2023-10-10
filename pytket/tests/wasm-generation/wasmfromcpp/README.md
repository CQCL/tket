Ways to compile wasm from C/C++ are:

```
emcc <filename>.c -o <filename>.html
```

```
clang --target=wasm32 --no-standard-libraries -Wl,--export-all -Wl,--no-entry -o <filename>.wasm <filename>.c
```

```
clang --target=wasm32 -mmultivalue -Xclang -target-abi -Xclang experimental-mv --no-standard-libraries -Wl,--export-all -Wl,--no-entry -o <filename>.wasm <filename>.c
```


You can then run:

```
wasm2wat <filename>.wasm -o <filename>.wast
```

to convert the WASM to human-readable text format.
