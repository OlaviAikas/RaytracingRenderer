#pragma once
#include "Light.hpp"
#include "Ball.hpp"
#include "Ray.hpp"
#include <vector>
#include <random>

class Scene {
    private:
        std::vector<Ball> balls;
        std::vector<Light> lights;

    public:
        Scene();
        void add_ball(Ball ball);
        void add_light(Light light);

        int numballs() const;
        int numlights() const;

        Vect colour(const Ray& ray, unsigned int depth, std::default_random_engine& engine, std::uniform_real_distribution<double>& distr) const;
};
