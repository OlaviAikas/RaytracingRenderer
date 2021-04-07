#include "headers/Vect.hpp"
#include "headers/Ray.hpp"
#include "headers/Ball.hpp"
#include "headers/Scene.hpp"
#include <iostream>
#include <vector>
#include <math.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define W 1024 //Width of camera sensor in pixels
#define H 1024 //Height of camera sensor in pixels
#define FOV 1.0 //In radians
#define NRAYS 3 //How many rays are checked per pixel
#define NJUMPS 12 //How many times a ray can bounce
#define GAMMA 2.2

void print(Vect v) {
    std::cout << "(" << v.get_x() << " " << v.get_y() << " " << v.get_z() << ")" << std::endl;
}

int main() {
    Ball b = Ball(Vect(0,1,-10), 1);
    Scene scene;
    srand(time(0));
    //Background balls
    scene.add_ball(Ball(Vect(0,0,-1000), 940, Vect(0,1,0)));
    scene.add_ball(Ball(Vect(0,0,1000), 940, 1));
    scene.add_ball(Ball(Vect(0,-1000,0), 940, Vect(0, 1, 1)));
    scene.add_ball(Ball(Vect(0,1000,0), 940, Vect(1, 1, 0)));
    scene.add_ball(Ball(Vect(-1000,0,0), 940, Vect(1, 0, 1)));
    scene.add_ball(Ball(Vect(950,0,0), 940, Vect(0, 0, 1)));
    //Target ball(s)
    scene.add_ball(Ball(Vect(0,10,0), 10, Vect(1,0,0), 1, 0.2, 1));
    scene.add_ball(Ball(Vect(0,-10,0), 10, Vect(0,0,1), 1, 0.2, 1));
    //scene.add_ball(Ball(Vect(0,0,-10), 1));
    //Light it up!
    scene.add_light(Light(Vect(-25, 0, 40), 38000));
    double f = (W/2.f)/tan(FOV/2.f);
    //Render the scene from the point of view of a camera at (0,0,55)
    //With the conventions in the course notes
    std::vector<unsigned char> image(W*H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vect dir = Vect(i + 0.5 - W/2.f, j + 0.5 - H/2.f, 55 - f) - Vect(0,0,55);
            double r_i = 0;
            double g_i = 0;
            double b_i = 0;
            for (int k = 0; k < NRAYS; k++) {
                Vect c = scene.colour(Ray(Vect(0,0,55), dir), NJUMPS);
                r_i += c.get_x();
                g_i += c.get_y();
                b_i += c.get_z();
            }
            //print(r_color);
            image[(i*W + j) * 3 + 0] = (unsigned char) 255*pow(r_i/NRAYS, 1/GAMMA); //R
            image[(i*W + j) * 3 + 1] = (unsigned char) 255*pow(g_i/NRAYS, 1/GAMMA); //G
            image[(i*W + j) * 3 + 2] = (unsigned char) 255*pow(b_i/NRAYS, 1/GAMMA); //B
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    return 0;
}