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
    Vect b(1,2,3);
    Vect bn = b.normalize();
    print(bn);
    print(b.cross(bn));
    std::cout << bn.norm() << std::endl;
    return 0;
}