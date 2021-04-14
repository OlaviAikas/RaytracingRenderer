#include "headers/Vect.hpp"
#include "headers/Ball.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <limits>
#include "headers/MonteCarloIntegral.hpp"
using namespace MonteCarloIntegral;

void print(Vect v) {
    std::cout << "(" << v.get_x() << " " << v.get_y() << " " << v.get_z() << ")" << std::endl;
}

double f(double v) {
    return v;
}

int main() {
    std::cout << integrate_gauss(f, 0, 1, 100) << std::endl;
    return 0;
}