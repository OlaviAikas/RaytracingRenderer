#include "../headers/Ball.hpp"
#include "../headers/Vect.hpp"
#include "../headers/Ray.hpp"
#include <math.h>

Ball::Ball(const Vect& pos, const double& r) 
    : pos(pos), r(r), 
      colour(Vect(255, 255, 255)), 
      alb(0.9), mirror(false), ri_in(-1), ri_out(1), k0_in(0), k0_out(0) { }

Ball::Ball(const Vect& pos, const double& r, const Vect& colour)
    : pos(pos), r(r),
      colour(colour),
      alb(0.9), mirror(false), ri_in(-1), ri_out(1), k0_in(0), k0_out(0) { }
    
Ball::Ball(const Vect& pos, const double& r, const Vect& colour, const double& alb)
    : pos(pos), r(r),
      colour(colour),
      alb(alb), mirror(false), ri_in(-1), ri_out(1), k0_in(0), k0_out(0) { }

Ball::Ball(const Vect& pos, const double& r, const double& ri_in)
    : pos(pos), r(r),
      colour(Vect(255,255,255)),
      alb(0.9), mirror(false), ri_in(ri_in), ri_out(1),
      k0_in(((ri_in - ri_out)*(ri_in - ri_out))/((ri_in + ri_out)*(ri_in + ri_out))), 
      k0_out(((ri_out - ri_in)*(ri_out - ri_in))/((ri_out + ri_in)*(ri_out + ri_in))) { }

Ball::Ball(const Vect& pos, const double& r, const double& ri_in, const double& ri_out)
    : pos(pos), r(r),
      colour(Vect(255,255,255)),
      alb(0.9), mirror(false), ri_in(ri_in), ri_out(ri_out), 
      k0_in(((ri_in - ri_out)*(ri_in - ri_out))/((ri_in + ri_out)*(ri_in + ri_out))), 
      k0_out(((ri_out - ri_in)*(ri_out - ri_in))/((ri_out + ri_in)*(ri_out + ri_in))) { }

Ball::Ball(const Vect& pos, const double& r, const bool& mirror)
    : pos(pos), r(r),
      colour(Vect(255,255,255)),
      alb(0.9), mirror(mirror), ri_in(-1), ri_out(1), k0_in(0), k0_out(0) { }

Vect Ball::get_pos() const {
    return pos;
}

double Ball::get_r() const {
    return r;
}

Vect Ball::get_colour() const {
    return colour;
}

double Ball::get_alb() const {
    return alb;
}

bool Ball::get_mirror() const {
    return mirror;
}

double Ball::get_ri_in() const {
    return ri_in;
}

double Ball::get_ri_out() const {
    return ri_out;
}

double Ball::get_k0_out() const {
    return k0_out;
}

double Ball::get_k0_in() const {
    return k0_in;
}

Ilist Ball::intersection(const Ray& ray) const {
    double p = ray.get_dir().dot(ray.get_start() - pos);
    double ocn = (ray.get_start()-pos).norm();
    double discriminant = p*p - (ocn*ocn - r*r);
    if (discriminant < 0) {
        return Ilist(0, Vect(0,0,0));
    }
    if (discriminant == 0) {
        double sol = -p - sqrt(discriminant);
        if (sol >= 0) {
            return Ilist(1, ray.get_start() + ray.get_dir()*sol);
        }
        return Ilist(0, Vect(0,0,0));
    }
    double sol1 = -p - sqrt(discriminant);
    double sol2 = -p + sqrt(discriminant);
    if (sol1 >= 0) { //Can see both
        return Ilist(1, ray.get_start() + ray.get_dir()*sol1);
    }
    if (sol2 >= 0) { //Can only see sol2, inside the ball
        return Ilist(2, ray.get_start() + ray.get_dir()*sol2);
    }
    //Can't see either
    return Ilist(0, Vect(0,0,0));
}