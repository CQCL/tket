void init() {}

typedef struct _twoints
{
    unsigned x;
    unsigned y;
} twoints;

twoints divmod(unsigned x, unsigned y)
{
    unsigned q, r;
    if (y == 0) {
        q = 0;
        r = 0;
    } else {
        q = x / y;
        r = x % y;
    }
    twoints result = {q, r};
    return result;
}
