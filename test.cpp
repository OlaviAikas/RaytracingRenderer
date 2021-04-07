#include "headers/Vect.hpp"
#include "headers/Ball.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <limits>

void print(Vect v) {
    std::cout << "(" << v.get_x() << " " << v.get_y() << " " << v.get_z() << ")" << std::endl;
}

int main() {
    for (int i = 0; i < 100; i++) {
        std::cout << rand() << std::endl;
    }
    return 0;
}