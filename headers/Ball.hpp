#pragma once
#include "Vect.hpp"
#include "Ray.hpp"

struct Ilist {
    unsigned char n; //How many intersections
    Vect i1; //Coordinates of intersections (if any).
    Vect i2; //i1 is closer to the ray origin than i2, or of the intersections are behind the camera, then the other way around
    Ilist(unsigned char n, Vect i1, Vect i2) : n(n), i1(i1), i2(i2) { };
};


class Ball {
    private:
        Vect pos;
        double r;
        Vect colour;
        double alb; //Albedo
        double refl; //Reflectivity. 1 is a mirror, 0 doesn't reflect
        double ri; //Index of refraction

    public:
        Ball(Vect pos, double r);
        Ball(Vect pos, double r, Vect colour);
        Ball(Vect pos, double r, Vect colour, double alb);
        Ball(Vect pos, double r, double refl);
        Ball(Vect pos, double r, double refl, double ri);
        Ball(Vect pos, double r, Vect colour, double alb, double refl, double ri);

        Vect get_pos() const;
        double get_r() const;
        Vect get_colour() const;
        double get_alb() const;
        double get_refl() const;
        double get_ri() const;
        //Get the VISIBLE intersections of the ball by the ray, e.g not behind the ray origin
        Ilist intersections(Ray ray);
};