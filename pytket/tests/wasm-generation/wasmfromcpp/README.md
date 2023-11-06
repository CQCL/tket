Ways to compile wasm from C/C++ are:

Using emcc:
```
emcc <filename>.c -o <filename>.html
```
This will generate an wasm file and a html file. You can ignore the html.

Using clang, this should work in all clang vesions.
```
clang --target=wasm32 --no-standard-libraries -Wl,--export-all -Wl,--no-entry -o <filename>.wasm <filename>.c
```

If you want to compile wasm files that contains function which return more than one int you need to give some more parameters, this only works with > clang 16

```
clang --target=wasm32 -mmultivalue -Xclang -target-abi -Xclang experimental-mv --no-standard-libraries -Wl,--export-all -Wl,--no-entry -o <filename>.wasm <filename>.c
```


After generating the wasm file you can then run:

```
wasm2wat <filename>.wasm -o <filename>.wast
```

to convert the WASM to human-readable text format.
