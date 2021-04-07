#pragma once

class Vect {
    public:
        Vect();
        Vect(double x, double y, double z);

        Vect operator+(const Vect other);
        Vect operator-(const Vect other);
        Vect operator/(const double scalar);
        Vect operator*(const double scalar);
        Vect& operator+=(const Vect other);
        double dot(const Vect other);
        Vect cross(const Vect other);
        double norm() const;
        Vect normalize() const;
        double get_x() const;
        double get_y() const;
        double get_z() const;
        void print() const;

    private:
        double coors[3];
};