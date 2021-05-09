#pragma once

class Vect {
    public:
        Vect();
        Vect(double x, double y, double z);

        Vect operator+(const Vect other) const;
        Vect operator-(const Vect other) const;
        Vect operator/(const double scalar) const;
        Vect operator*(const double scalar) const;
        Vect operator*(const Vect other) const;
        Vect& operator+=(const Vect other);
        double dot(const Vect other) const;
        Vect cross(const Vect other) const;
        double norm() const;
        Vect normalize() const;
        double get_x() const;
        double get_y() const;
        double get_z() const;
        void print() const;
        double operator[](int i) const;
        double& operator[](int i);
        bool operator==(const Vect other) const;

    private:
        double coors[3];
};