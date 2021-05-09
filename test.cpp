#include "headers/Vect.hpp"
#include "headers/Ball.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <limits>
#include <random>
#include "headers/Tmesh.hpp"

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
    Vect normal = Vect(1, 2, 0).normalize();
    double r1 = uniform(engine);
    double r2 = uniform(engine);
    double sqrt1mr2 = sqrt(1-r2);
    double twopir1 = 2*M_PI*r1;
    double x = cos(twopir1)*sqrt1mr2;
    double y = sin(twopir1)*sqrt1mr2;
    double z = sqrt(r2);
    Vect T1(0,0,0);
    Vect T2(0,0,0);
    if (abs(x) <= abs(y), abs(x) <= abs(z)) {
        T1 = Vect(0, -z, y).normalize();
    } else if (abs(y) <= abs(x), abs(y) <= abs(z)) {
        T1 = Vect(-z, 0, x).normalize();
    } else {
        T1 = Vect(-y, x, 0).normalize();
    }
    T2 = T1.cross(normal);
    print(T1);
    print(T2);
    print(normal);
    std::cout << normal.dot(T2) << std::endl;
    std::cout << T1.dot(T2) << std::endl;
    return 0;
}