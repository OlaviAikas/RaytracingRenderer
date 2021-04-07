#pragma once
#include "Vect.hpp"

class Light {
    private:
        Vect pos;
        double intensity;

    public:
        Light(Vect pos, double intensity);

        Vect get_pos() const;
        double get_intensity() const;
};