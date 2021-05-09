#pragma once
#include "Vect.hpp"
#include "Ray.hpp"

struct Intersection {
    bool going_in;
    Vect pos; //Coordinates of the intersection
    Vect normal; //Normal vector (0 if no intersection)
    Intersection(Vect pos, Vect normal, bool going_in) : pos(pos), normal(normal), going_in(going_in) { };
};

class Geometry {
    public:
        Geometry(const Vect& colour);
        const Vect colour; //Albedo
        virtual Intersection intersection(const Ray& ray) const = 0;
        virtual bool is_mirror() const = 0;
        virtual bool is_transparent() const = 0;
        virtual double get_ri_in() const = 0;
        virtual double get_ri_out() const = 0;
        virtual double get_k0_in() const = 0;
        virtual double get_k0_out() const = 0;
};