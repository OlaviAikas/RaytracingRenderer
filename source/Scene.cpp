#include "../headers/Scene.hpp"
#include "../headers/Geometry.hpp"
#include "../headers/Light.hpp"
#include "../headers/Ray.hpp"
#include <random>
#include <chrono>
#include <vector>
#include <math.h>
#include <limits>
#include <iostream>

#define EPS 0.00000001

Scene::Scene() { }

void Scene::free_geos() {
    for (int i = 0; i < numgeos(); i++) {
        free(geos[i]);
    }
}

void Scene::add_geo(Geometry* geo) {
    geos.push_back(geo);
}

void Scene::add_light(Light light) {
    lights.push_back(light);
}

int Scene::numgeos() const {
    return geos.size();
}

int Scene::numlights() const {
    return lights.size();
}

Vect Scene::colour(const Ray& ray, unsigned int depth, std::default_random_engine& engine, std::uniform_real_distribution<double>& distr) const {
    if (depth <= 0) {
        return Vect(0,0,0);
    }
    int ngeos = numgeos();
    int nlights = numlights();
    Vect closest_hit(INFINITY, INFINITY, INFINITY);
    double closest_hit_d = INFINITY;
    bool going_in = true;
    int hitn = -1;
    Vect normal;
    for (int i = 0; i < ngeos; i++) {
        Intersection l = geos[i]->intersection(ray);
        if (l.normal == Vect(0,0,0)) { continue; }
        if ((l.pos - ray.get_start()).norm() < closest_hit_d) {
            closest_hit = l.pos;
            normal = l.normal;
            going_in = l.going_in;
            closest_hit_d = (closest_hit - ray.get_start()).norm();
            hitn = i;
        }
    }
    if (hitn == -1) {
        return Vect(0, 0, 0);
    }
    bool tot_int_refl = false;
    //Refraction
    if (geos[hitn]->is_transparent()) {
        double k0 = 0;
        if (going_in) {
            k0 = geos[hitn]->get_k0_out();
        } else {
            k0 = geos[hitn]->get_k0_in();
        }
        double rc = k0 + (1 - k0)*pow((1-abs(normal.dot(ray.get_dir()))), 5);
        if (distr(engine) < rc) {
            Vect wr = ray.get_dir() - normal*(2*ray.get_dir().dot(normal));
            closest_hit += normal*EPS;
            return colour(Ray(closest_hit, wr), depth - 1, engine, distr);
        }
        double n1on2 = 0;
        if (going_in) {
            n1on2 = geos[hitn]->get_ri_out()/geos[hitn]->get_ri_in();
        } else {
            n1on2 = geos[hitn]->get_ri_in()/geos[hitn]->get_ri_out();
        }
        double rdn = ray.get_dir().dot(normal);
        double rt = 1 - (n1on2*n1on2)*(1 - (rdn*rdn));
        if (rt < 0) {
            tot_int_refl = true;
        } else {
            Vect wtT = (ray.get_dir() - (normal*rdn))*n1on2;
            Vect wtN = normal*-1*sqrt(rt);
            closest_hit = closest_hit - normal*EPS;
            return colour(Ray(closest_hit, wtT + wtN), depth - 1, engine, distr);
        }
    }
    //Reflection
    if (tot_int_refl || geos[hitn]->is_mirror()) {
        //std::cout << "reflecting" << std::endl;
        Vect wr = ray.get_dir() - normal*(2*ray.get_dir().dot(normal));
        closest_hit += normal*EPS;
        return colour(Ray(closest_hit, wr), depth - 1, engine, distr);
    }
    //Direct lighting
    closest_hit += normal*EPS;
    Vect L0(0,0,0);
    for (int i = 0; i < nlights; i++) {
        Vect line_of_sight = lights[i].get_pos() - closest_hit;
        double sd = line_of_sight.norm();
        bool visible = true;
        for (int b = 1; b < ngeos; b++) {
            Intersection l = geos[b]->intersection(Ray(closest_hit, line_of_sight));
            if (l.normal == Vect(0,0,0)) { continue; }
            if ((l.pos - closest_hit).norm() < sd) {
                visible = false;
                break;
            }
        }
        if (visible == false) {
            //std::cout << "This hit can't see the light" << std::endl;
            continue;
        }
        Vect ct = (geos[hitn]->colour*lights[i].get_intensity())/(4*M_PI*M_PI*sd*sd);
        double p = std::max(normal.dot((line_of_sight/sd)), 0.);
        L0 += ct*p;
    }
    //Indirect lighting random sample
    double r1 = distr(engine);
    double r2 = distr(engine);
    double sqrt1mr2 = sqrt(1-r2);
    double twopir1 = 2*M_PI*r1;
    double x = cos(twopir1)*sqrt1mr2;
    double y = sin(twopir1)*sqrt1mr2;
    double z = sqrt(r2);
    Vect T1(0,0,0);
    Vect T2(0,0,0);
    if (abs(x) <= abs(y), abs(x) <= abs(z)) {
        T1 = Vect(0, -z, y).normalize();
    } else if (abs(y) <= abs(x), abs(y) <= abs(z)) {
        T1 = Vect(-z, 0, x).normalize();
    } else {
        T1 = Vect(-y, x, 0).normalize();
    }
    T2 = T1.cross(normal);
    Vect wi = T1*x + T2*y + normal*z;
    L0 += geos[hitn]->colour*colour(Ray(closest_hit, wi), depth - 1, engine, distr);
    return L0;
}