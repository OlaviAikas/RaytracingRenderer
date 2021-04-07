#pragma once
#include "Light.hpp"
#include "Ball.hpp"
#include "Ray.hpp"
#include <vector>

class Scene {
    private:
        std::vector<Ball> balls;
        std::vector<Light> lights;

    public:
        Scene();
        void add_ball(Ball ball);
        void add_light(Light light);

        Ball* get_balls();
        Light* get_lights();

        int numballs();
        int numlights();

        Vect colour(Ray ray);
};
