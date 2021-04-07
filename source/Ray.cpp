#include "../headers/Ray.hpp"
#include "../headers/Vect.hpp"


Ray::Ray(Vect start, Vect dir) : start(start), 
                                 dir(dir.normalize()) { }

Vect Ray::get_start() const {
    return start;
}

Vect Ray::get_dir() const {
    return dir;
}
