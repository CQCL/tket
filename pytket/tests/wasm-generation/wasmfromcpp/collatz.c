void init() {}

unsigned collatz(unsigned n) {
    if (n == 0) {
        return 0;
    }
    unsigned m = 0;
    while (n != 1) {
        if (n & 1) {
            n = (3 * n + 1) / 2;
        }
        else {
            n /= 2;
        }
        m++;
    }
    return m;
}
