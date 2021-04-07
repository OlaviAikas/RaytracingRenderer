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

Scene::Scene() { }

void Scene::add_ball(Ball ball) {
    balls.push_back(ball);
}

void Scene::add_light(Light light) {
    lights.push_back(light);
}

Ball* Scene::get_balls() {
    return balls.data();
}

Light* Scene::get_lights() {
    return lights.data();
}

int Scene::numballs() {
    return balls.size();
}

int Scene::numlights() {
    return lights.size();
}

Vect Scene::colour(Ray ray, unsigned int depth) {
    int nballs = numballs();
    int nlights = numlights();
    Ball* balls = get_balls();
    Light* lights = get_lights();
    Vect closest_hit(INFINITY, INFINITY, INFINITY);
    int hitbn = -1;
    for (int i = 0; i < nballs; i++) {
        Ilist l = balls[i].intersections(ray);
        if (l.n == 0) { continue; }
        if (l.i1.norm() < closest_hit.norm()) {
            closest_hit = l.i1;
            hitbn = i;
        }
    }
    if (hitbn == -1) {
        return Vect(0, 0, 0);
    }
    Vect normal = (closest_hit - balls[hitbn].get_pos()).normalize();
    closest_hit += normal*0.0000001;
    if (balls[hitbn].get_refl() > 0 && depth > 0) {
        if (rand()/ ((double) RAND_MAX) < balls[hitbn].get_refl()) {
            Vect wr = ray.get_dir() - normal*(2*ray.get_dir().dot(normal));
            return colour(Ray(closest_hit, wr), depth - 1);
        }
    }
    double tot_intensity = 0;
    for (int i = 0; i < nlights; i++) {
        Vect line_of_sight = lights[i].get_pos() - closest_hit;
        double sd = line_of_sight.norm();
        bool visible = true;
        for (int b = 1; b < nballs; b++) {
            Ilist l = balls[b].intersections(Ray(closest_hit, line_of_sight));
            if (l.n == 0) { continue; }
            if ((l.i1 - closest_hit).norm() < sd) {
                //std::cout << "entered here" << std::endl;
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
    return balls[hitbn].get_colour()*fmin(1, tot_intensity);
}