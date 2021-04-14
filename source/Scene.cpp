#include "../headers/Scene.hpp"
#include "../headers/Ball.hpp"
#include "../headers/Light.hpp"
#include "../headers/Ray.hpp"
#include <random>
#include <chrono>
#include <vector>
#include <math.h>
#include <limits>
#include <iostream>

#define EPS 0.000000001
#define GAMMA 2.2

Scene::Scene() { }

void Scene::add_ball(Ball ball) {
    balls.push_back(ball);
}

void Scene::add_light(Light light) {
    lights.push_back(light);
}

int Scene::numballs() const {
    return balls.size();
}

int Scene::numlights() const {
    return lights.size();
}

Vect Scene::colour(const Ray& ray, unsigned int depth, bool& break_early) const {
    if (depth <= 0) {
        return Vect(0,0,0);
    }
    int nballs = numballs();
    int nlights = numlights();
    Vect closest_hit(INFINITY, INFINITY, INFINITY);
    double closest_hit_d = INFINITY;
    bool going_in = true;
    int hitbn = -1;
    for (int i = 0; i < nballs; i++) {
        Ilist l = balls[i].intersection(ray);
        if (l.n == 0) { continue; }
        if ((l.i - ray.get_start()).norm() < closest_hit_d) {
            closest_hit = l.i;
            closest_hit_d = (closest_hit - ray.get_start()).norm();
            if (l.n == 1) {
                going_in = true;
            } else {
                going_in = false;
            }
            hitbn = i;
        }
    }
    if (hitbn == -1) {
        return Vect(0, 0, 0);
    }
    Vect normal = (closest_hit - balls[hitbn].get_pos()).normalize();
    bool be = false;
    if (!going_in) {
        normal = normal*-1;
    }
    bool tot_int_refl = false;
    //Refraction
    if (balls[hitbn].get_ri_in() != -1) {
        double k0 = 0;
        if (going_in) {
            k0 = balls[hitbn].get_k0_out();
        } else {
            k0 = balls[hitbn].get_k0_in();
        }
        double rc = k0 + (1 - k0)*pow((1-abs(normal.dot(ray.get_dir()))), 5);
        if (rand()/((double) RAND_MAX) < rc) {
            Vect wr = ray.get_dir() - normal*(2*ray.get_dir().dot(normal));
            closest_hit += normal*EPS;
            return colour(Ray(closest_hit, wr), depth - 1, be);
        }
        double n1on2 = 0;
        if (going_in) {
            n1on2 = balls[hitbn].get_ri_out()/balls[hitbn].get_ri_in();
        } else {
            n1on2 = balls[hitbn].get_ri_in()/balls[hitbn].get_ri_out();
        }
        double rdn = ray.get_dir().dot(normal);
        double rt = 1 - (n1on2*n1on2)*(1 - (rdn*rdn));
        if (rt < 0) {
            tot_int_refl = true;
        } else {
            Vect wtT = (ray.get_dir() - (normal*rdn))*n1on2;
            Vect wtN = normal*-1*sqrt(rt);
            closest_hit = closest_hit - normal*EPS;
            return colour(Ray(closest_hit, wtT + wtN), depth - 1, be);
        }
    }
    if (tot_int_refl || balls[hitbn].get_mirror() == true) {
        //std::cout << "reflecting" << std::endl;
        Vect wr = ray.get_dir() - normal*(2*ray.get_dir().dot(normal));
        closest_hit += normal*EPS;
        return colour(Ray(closest_hit, wr), depth - 1, be);
    }
    break_early = true;
    closest_hit += normal*EPS;
    double tot_intensity = 0;
    for (int i = 0; i < nlights; i++) {
        Vect line_of_sight = lights[i].get_pos() - closest_hit;
        double sd = line_of_sight.norm();
        bool visible = true;
        for (int b = 1; b < nballs; b++) {
            Ilist l = balls[b].intersection(Ray(closest_hit, line_of_sight));
            if (l.n == 0) { continue; }
            if ((l.i - closest_hit).norm() < sd) {
                visible = false;
                break;
            }
        }
        if (visible == false) {
            //std::cout << "This hit can't see the light" << std::endl;
            continue;
        }
        double ct = (lights[i].get_intensity()*balls[hitbn].get_alb())/(4*M_PI*M_PI*sd*sd);
        double p = (closest_hit - balls[hitbn].get_pos()).normalize().dot((line_of_sight/sd));
        tot_intensity += ct*p;
    }
    return balls[hitbn].get_colour()*tot_intensity;
}