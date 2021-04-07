#include "../headers/Vect.hpp"
#include <iostream>
#include <math.h>

Vect::Vect() {
    coors[0] = 0;
    coors[1] = 0;
    coors[2] = 0;
}

Vect::Vect(double x, double y, double z) {
    coors[0] = x;
    coors[1] = y;
    coors[2] = z;
}

double Vect::get_x() const {
    return coors[0];
}

double Vect::get_y() const {
    return coors[1];
}

double Vect::get_z() const {
    return coors[2];
}

Vect Vect::operator+(const Vect other) {
    return Vect(coors[0] + other.get_x(), coors[1] + other.get_y(), coors[2] + other.get_z());
}

Vect Vect::operator-(const Vect other) {
    return Vect(coors[0] - other.get_x(), coors[1] - other.get_y(), coors[2] - other.get_z());
}

Vect& Vect::operator+=(const Vect other) {
    coors[0] += other.get_x();
    coors[1] += other.get_y();
    coors[2] += other.get_z();
    return *this;
}

Vect Vect::operator/(const double scalar) {
    return Vect(coors[0]/scalar, coors[1]/scalar, coors[2]/scalar);
}

Vect Vect::operator*(const double scalar) {
    return Vect(coors[0]*scalar, coors[1]*scalar, coors[2]*scalar);
}

double Vect::dot(const Vect other) {
    return coors[0]*other.get_x() + coors[1]*other.get_y() + coors[2]*other.get_z();
}

Vect Vect::cross(const Vect other) {
    return Vect(coors[1]*other.get_z() - coors[2]*other.get_y(), coors[2]*other.get_x() - coors[0]*other.get_z(), coors[0]*other.get_y() - coors[1]*other.get_x());
}

double Vect::norm() const {
    return std::sqrt(coors[0]*coors[0] + coors[1]*coors[1] + coors[2]*coors[2]);
}

Vect Vect::normalize() const {
    double n = this->norm();
    return(Vect(coors[0]/n, coors[1]/n, coors[2]/n));
}