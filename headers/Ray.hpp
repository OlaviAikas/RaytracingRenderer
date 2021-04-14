#pragma once
#include "Vect.hpp"

class Ray {
    private:
        const Vect start;
        const Vect dir;

    public:
        Ray(const Vect& start, const Vect& dir);

        Vect get_start() const;
        Vect get_dir() const;
        Vect get_colour() const;
};