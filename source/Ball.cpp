#include "../headers/Ball.hpp"
#include "../headers/Vect.hpp"
#include "../headers/Ray.hpp"
#include <math.h>

Ball::Ball(Vect pos, double r) 
    : pos(pos), r(r), 
      colour(Vect(255, 255, 255)), 
      alb(1), refl(0), ri(1) { }

Ball::Ball(Vect pos, double r, Vect colour)
    : pos(pos), r(r),
      colour(colour),
      alb(0.9), refl(0), ri(1) { }
    
Ball::Ball(Vect pos, double r, Vect colour, double alb)
    : pos(pos), r(r),
      colour(colour),
      alb(alb), refl(0), ri(1) { }

Ball::Ball(Vect pos, double r, double refl)
    : pos(pos), r(r),
      colour(Vect(255,255,255)),
      alb(0.9), refl(refl), ri(1) { }

Ball::Ball(Vect pos, double r, double refl, double ri)
    : pos(pos), r(r),
      colour(Vect(255,255,255)),
      alb(0.9), refl(refl), ri(ri) { }

Ball::Ball(Vect pos, double r, Vect colour, double alb, double refl, double ri)
    : pos(pos), r(r),
      colour(colour),
      alb(alb), refl(refl), ri(ri) { }


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

double Ball::get_refl() const {
    return refl;
}

double Ball::get_ri() const {
    return ri;
}

Ilist Ball::intersections(const Ray ray) {
    double p =  ray.get_dir().dot(ray.get_start() - pos);
    double ocn = (ray.get_start()-pos).norm();
    double discriminant = p*p - (ocn*ocn - r*r);
    if (discriminant < 0) {
        return Ilist(0, Vect(0,0,0), Vect(0,0,0));
    }
    if (discriminant == 0) {
        double sol = -p - sqrt(discriminant);
        if (sol >= 0) {
            return Ilist(1, ray.get_start() + ray.get_dir()*sol, Vect(0,0,0));
        }
        return Ilist(0, Vect(0,0,0), Vect(0,0,0));
    }
    double sol1 = -p - sqrt(discriminant);
    double sol2 = -p + sqrt(discriminant);
    if (sol1 <= sol2) {
        if (sol1 >= 0) { //Can see both
            return Ilist(2, ray.get_start() + ray.get_dir()*sol1, ray.get_start() + ray.get_dir()*sol2);
        }
        if (sol2 >= 0) { //Can only see sol2
            return Ilist(1, ray.get_start() + ray.get_dir()*sol2, Vect(0,0,0));
        }
        //Can't see either
        return Ilist(0, Vect(0,0,0), Vect(0,0,0));
    }
    if (sol2 <= sol1) {
        if (sol2 >= 0) { //Can see both
            return Ilist(2, ray.get_start() + ray.get_dir()*sol2, ray.get_start() + ray.get_dir()*sol1);
        }
        if (sol1 >= 0) { //Can only see sol1
            return Ilist(1, ray.get_start() + ray.get_dir()*sol1, Vect(0,0,0));
        }
        //Can't see either
        return Ilist(0, Vect(0,0,0), Vect(0,0,0));
    }
    return Ilist(0, Vect(0,0,0), Vect(0,0,0));
}