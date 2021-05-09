#pragma once
#include "Vect.hpp"
#include "Ray.hpp"
#include "Geometry.hpp"

class Ball : public Geometry {
    private:
        const Vect pos;
        const double r;
        const bool mirror; //Reflectivity. 1 is a mirror, 0 doesn't reflect
        const double ri_in; //Index of refraction inside
        const double ri_out; //Index of refraction outside
        const double k0_in; //reflection coefficient for transparent balls inside
        const double k0_out; //reflection coefficient for transparent balls outside

    public:
        Ball(const Vect& pos, const double& r);
        Ball(const Vect& pos, const double& r, const Vect& colour);
        Ball(const Vect& pos, const double& r, const double& ri_in);
        Ball(const Vect& pos, const double& r, const double& ri_in, const double& ri_out);
        Ball(const Vect& pos, const double& r, const bool& mirror);

        Vect get_pos() const;
        double get_r() const;
        Vect get_colour() const;
        bool is_mirror() const;
        bool is_transparent() const;
        double get_ri_in() const;
        double get_ri_out() const;
        double get_k0_in() const;
        double get_k0_out() const;
        //Get the nearest VISIBLE intersections of the ball by the ray, 
        //e.g not behind the ray origin and indicate if the ray is going in or out
        Intersection intersection(const Ray& ray) const;
};