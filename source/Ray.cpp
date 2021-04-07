#include "../headers/Ray.hpp"
#include "../headers/Vect.hpp"


Ray::Ray(Vect start, Vect dir) : start(start), 
                                 dir(dir.normalize()), 
                                 colour(Vect(255, 255, 255)) { }

Vect Ray::get_start() const {
    return start;
}

Vect Ray::get_dir() const {
    return dir;
}

Vect Ray::get_colour() const {
    return colour;
}
