#pragma once
#include "Light.hpp"
#include "Geometry.hpp"
#include "Ray.hpp"
#include <vector>
#include <random>

class Scene {
    private:
        std::vector<Geometry*> geos;
        std::vector<Light> lights;

    public:
        Scene();
        void free_geos();
        void add_geo(Geometry* geo);
        void add_light(Light light);

        int numgeos() const;
        int numlights() const;

        Vect colour(const Ray& ray, unsigned int depth, std::default_random_engine& engine, std::uniform_real_distribution<double>& distr) const;
};
