#include <emscripten/emscripten.h>
#include <stdio.h>

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN
#endif

int global_v = 0;
int global_c = 0;
int global_b = 0;

EXTERN EMSCRIPTEN_KEEPALIVE void init() {
  global_v = 0;
  global_c = 0;
  global_b = 8;
}

EXTERN EMSCRIPTEN_KEEPALIVE void myFunction(int argc) {
  printf("Hello World\n");
}

EXTERN EMSCRIPTEN_KEEPALIVE int myFunction2(
    int value, int value2, int value3, int value4, int value5, int value6,
    int value7) {
  return value + value2;
}

EXTERN EMSCRIPTEN_KEEPALIVE int myFunction3(
    int value,  //
    int value2, int value3, int value4, int value5, int value6, int value7,
    int value8, int value9, int value10, int value11) {
  return value + value2;
}

EXTERN EMSCRIPTEN_KEEPALIVE int mightloop_returns1() {
  if (global_c == 0) {
    return 1;
  } else {
    while (global_c > 0) {
      global_c = global_c * 1;
    }
    return 2;
  }
}

EXTERN EMSCRIPTEN_KEEPALIVE int longrun() {
  for (int i = 0; i < 100000; ++i) {
    global_b = global_b * 2;
    global_b = global_b * 2;
    global_b = global_b / 8;
    global_b = global_b * 2;

    global_v = global_v + 1;
  }
  return global_v;
}

EXTERN EMSCRIPTEN_KEEPALIVE int verylongrun() {
  for (int i = 0; i < 100000000; ++i) {
    global_b = global_b * 2;
    global_b = global_b * 2;
    global_b = global_b / 8;
    global_b = global_b * 2;

    global_v = global_v + 1;
  }
  return global_v;
}

EXTERN EMSCRIPTEN_KEEPALIVE int get_v() { return global_v; }

EXTERN EMSCRIPTEN_KEEPALIVE int get_c() { return global_c; }

EXTERN EMSCRIPTEN_KEEPALIVE int get_b() { return global_b; }
