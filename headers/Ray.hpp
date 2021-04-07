#pragma once
#include "Vect.hpp"

class Ray {
    private:
        Vect start;
        Vect dir;
        Vect colour;

    public:
        Ray(Vect start, Vect dir);

        Vect get_start() const;
        Vect get_dir() const;
        Vect get_colour() const;
};