#include "../headers/Light.hpp"
#include "../headers/Vect.hpp"

Light::Light(Vect pos, double intensity)
    : pos(pos), intensity(intensity) { }

Vect Light::get_pos() const {
    return pos;
}

double Light::get_intensity() const {
    return intensity;
}