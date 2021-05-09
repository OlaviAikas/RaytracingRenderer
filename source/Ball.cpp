#include "../headers/Ball.hpp"
#include "../headers/Vect.hpp"
#include "../headers/Ray.hpp"
#include <math.h>

Ball::Ball(const Vect& pos, const double& r) 
    : Geometry(Vect(1,1,1)), pos(pos), r(r),
      mirror(false), ri_in(-1), ri_out(1), k0_in(0), k0_out(0) { }

Ball::Ball(const Vect& pos, const double& r, const Vect& colour)
    : Geometry(colour), pos(pos), r(r),
      mirror(false), ri_in(-1), ri_out(1), k0_in(0), k0_out(0) { }

Ball::Ball(const Vect& pos, const double& r, const double& ri_in)
    : Geometry(Vect(1,1,1)), pos(pos), r(r),
      mirror(false), ri_in(ri_in), ri_out(1),
      k0_in(((ri_in - ri_out)*(ri_in - ri_out))/((ri_in + ri_out)*(ri_in + ri_out))), 
      k0_out(((ri_out - ri_in)*(ri_out - ri_in))/((ri_out + ri_in)*(ri_out + ri_in))) { }

Ball::Ball(const Vect& pos, const double& r, const double& ri_in, const double& ri_out)
    : Geometry(Vect(1,1,1)), pos(pos), r(r),
      mirror(false), ri_in(ri_in), ri_out(ri_out), 
      k0_in(((ri_in - ri_out)*(ri_in - ri_out))/((ri_in + ri_out)*(ri_in + ri_out))), 
      k0_out(((ri_out - ri_in)*(ri_out - ri_in))/((ri_out + ri_in)*(ri_out + ri_in))) { }

Ball::Ball(const Vect& pos, const double& r, const bool& mirror)
    : Geometry(Vect(1,1,1)), pos(pos), r(r),
      mirror(mirror), ri_in(-1), ri_out(1), k0_in(0), k0_out(0) { }

Vect Ball::get_pos() const {
    return pos;
}

double Ball::get_r() const {
    return r;
}

Vect Ball::get_colour() const {
    return colour;
}

bool Ball::is_mirror() const {
    return mirror;
}

bool Ball::is_transparent() const {
    if (ri_in == -1) {
        return false;
    }
    return true;
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

Intersection Ball::intersection(const Ray& ray) const {
    double p = ray.get_dir().dot(ray.get_start() - pos);
    double ocn = (ray.get_start()-pos).norm();
    double discriminant = p*p - (ocn*ocn - r*r);
    if (discriminant < 0) {
        return Intersection(Vect(0,0,0), Vect(0,0,0), false);
    }
    if (discriminant == 0) {
        double sol = -p;
        if (sol >= 0) {
            Vect int_pos = ray.get_start() + ray.get_dir()*sol;
            return Intersection(int_pos, (int_pos - pos).normalize(), true);
        }
        return Intersection(Vect(0,0,0), Vect(0,0,0), false);
    }
    double sol1 = -p - sqrt(discriminant);
    double sol2 = -p + sqrt(discriminant);
    if (sol1 >= 0) { //Can see both
        Vect int_pos = ray.get_start() + ray.get_dir()*sol1;
        return Intersection(int_pos, (int_pos - pos).normalize(), true);
    }
    if (sol2 >= 0) { //Can only see sol2, inside the ball
        Vect int_pos = ray.get_start() + ray.get_dir()*sol2;
        return Intersection(int_pos, (int_pos - pos).normalize()*-1, false);
    }
    //Can't see either
    return Intersection(Vect(0,0,0), Vect(0,0,0), false);
}