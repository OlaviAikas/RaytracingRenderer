#include "headers/Vect.hpp"
#include "headers/Ball.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <limits>
#include <random>

void print(Vect v) {
    std::cout << "(" << v.get_x() << " " << v.get_y() << " " << v.get_z() << ")" << std::endl;
}

static std::default_random_engine engine(time(0));
static std::uniform_real_distribution<double> uniform(0, 1);

void box_muller(double& x, double& y) {
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    x = sqrt(-2 * log(r1))*cos(2*M_PI*r2);
    y = sqrt(-2 * log(r1))*sin(2*M_PI*r2);
}

double gpdf(const double& x, const double& y, const double& z) {
    return pow(1/(sqrt(2*M_PI)), 3)*exp(-(x*x + y*y + z*z)/2);
}

int main() {
    double sum = 0;
    for (unsigned s = 0; s < 50000; s++) {
        double x1, x2, y1, y2, z1, z2;
        box_muller(x1, x2);
        box_muller(y1, y2);
        box_muller(z1, z2);
        double t1 = 1;
        double t2 = 1;
        if (x1 < -M_PI_2 || y1 < -M_PI_2 || z1 < -M_PI_2 || x1 > M_PI_2 || y1 > M_PI_2 || z1 > M_PI_2) {
            t1 = 0;
        }
        if (x2 < -M_PI_2 || y2 < -M_PI_2 || z2 < -M_PI_2 || x2 > M_PI_2 || y2 > M_PI_2 || z2 > M_PI_2) {
            t2 = 0;
        }
        sum += t1*cos(x1*y1*z1)/gpdf(x1, y1, z1);
        sum += t2*cos(x2*y2*z2)/gpdf(x2, y2, z2);
    }
    std::cout << sum/100000 << std::endl;
    return 0;
}