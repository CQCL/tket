int global_v = 0;
int global_c = 0;
int global_b = 0;

void init() {
  global_v = 0;
  global_c = 0;
  global_b = 8;
}

void myFunction(int argc) { global_v = 0; }

int myFunction2(
    int value,   //
    int value2,  //
    int value3,  //
    int value4,  //
    int value5,  //
    int value6,  //
    int value7) {
  return value + value2;
}

int myFunction3(
    int value,    //
    int value2,   //
    int value3,   //
    int value4,   //
    int value5,   //
    int value6,   //
    int value7,   //
    int value8,   //
    int value9,   //
    int value10,  //
    int value11) {
  return value + value2;
}

int mightloop_returns1() {
  if (global_c == 0) {
    return 1;
  } else {
    while (global_c > 0) {
      global_c = global_c * 1;
    }
    return 2;
  }
}

int longrun() {
  for (int i = 0; i < 100000; ++i) {
    global_b = global_b * 2;
    global_b = global_b * 2;
    global_b = global_b / 8;
    global_b = global_b * 2;

    global_v = global_v + 1;
  }
  return global_v;
}

int verylongrun() {
  for (int i = 0; i < 100000000; ++i) {
    global_b = global_b * 2;
    global_b = global_b * 2;
    global_b = global_b / 8;
    global_b = global_b * 2;

    global_v = global_v + 1;
  }
  return global_v;
}

int get_v() { return global_v; }

int get_c() { return global_c; }

int get_b() { return global_b; }
