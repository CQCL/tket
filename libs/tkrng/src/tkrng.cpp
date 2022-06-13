#include <iostream>
#include "tkrng.h"

void tkrng(){
    #ifdef NDEBUG
    std::cout << "tkrng/0.1.0: Hello World Release!" <<std::endl;
    #else
    std::cout << "tkrng/0.1.0: Hello World Debug!" <<std::endl;
    #endif
}
