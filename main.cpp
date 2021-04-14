#include "headers/Vect.hpp"
#include "headers/Ray.hpp"
#include "headers/Ball.hpp"
#include "headers/Scene.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <thread>
#include <algorithm>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define W 1024 //Width of camera sensor in pixels
#define H 1024 //Height of camera sensor in pixels
#define FOV 1.1 //In radians
#define NRAYS 100 //How many rays are checked per pixel
#define NJUMPS 5 //How many times a ray can bounce
#define GAMMA 2.2
#define NTHREADS 8

void print(Vect v) {
    std::cout << "(" << v.get_x() << " " << v.get_y() << " " << v.get_z() << ")" << std::endl;
}

void render_rows(const Scene& s, const double& f, const unsigned& thread_n, std::vector<unsigned char>& image) {
    for (unsigned i = thread_n; i < H; i += NTHREADS) {
        for (unsigned j = 0; j < W; j++) {
            Vect dir = Vect(i + 0.5 - W/2.f, j + 0.5 - H/2.f, 55 - f) - Vect(0,0,55);
            double r_i = 0;
            double g_i = 0;
            double b_i = 0;
            bool break_early = false;
            for (int k = 0; k < NRAYS; k++) {
                Vect c = s.colour(Ray(Vect(0,0,55), dir), NJUMPS, break_early);
                r_i += c.get_x();
                g_i += c.get_y();
                b_i += c.get_z();
                if (break_early) {
                    break;
                }
            }
            double sent_rays = NRAYS;
            if (break_early) {
                sent_rays = 1;
            }
            //print(r_color);
            image[(i*W + j) * 3 + 0] = (unsigned char) std::min(255, (int) pow((r_i/sent_rays), 1/GAMMA)); //R
            image[(i*W + j) * 3 + 1] = (unsigned char) std::min(255, (int) pow((g_i/sent_rays), 1/GAMMA)); //G
            image[(i*W + j) * 3 + 2] = (unsigned char) std::min(255, (int) pow((b_i/sent_rays), 1/GAMMA)); //B
        }
    }
}

int main() {
    Scene scene;
    srand(time(0));
    //Background balls
    scene.add_ball(Ball(Vect(0,0,-1000), 940, Vect(0,1,0)));
    scene.add_ball(Ball(Vect(0,0,1000), 940, Vect(1,0,1)));
    scene.add_ball(Ball(Vect(0,-1000,0), 940, Vect(0, 1, 1)));
    scene.add_ball(Ball(Vect(0,1000,0), 940, Vect(1, 1, 0)));
    scene.add_ball(Ball(Vect(-1000,0,0), 940, Vect(1, 0, 0)));
    scene.add_ball(Ball(Vect(950,0,0), 940, Vect(0, 0, 1)));
    //Target ball(s)
    scene.add_ball(Ball(Vect(0, -20, 0), 10, true));
    scene.add_ball(Ball(Vect(0,0,0), 10, 1.5));
    scene.add_ball(Ball(Vect(0,20,0), 10, 1.5));
    scene.add_ball(Ball(Vect(0,20,0), 9.2, 1., 1.5));
    //Light it up!
    scene.add_light(Light(Vect(-20, -10, 40),20000000000));
    double f = (W/2.f)/tan(FOV/2.f);
    //Render the scene from the point of view of a camera at (0,0,55)
    //With the conventions in the course notes
    std::cout << "Starting render" << std::endl;
    time_t t0 = time(0);
    std::vector<unsigned char> image(W*H * 3, 0);
    std::thread threads[NTHREADS-1];
    for (unsigned k = 0; k < NTHREADS - 1; k++) {
        threads[k] = std::thread(render_rows, scene, f, k, std::ref(image));
    }
    render_rows(scene, f, NTHREADS-1, std::ref(image));
    for (unsigned k = 0; k < NTHREADS - 1; k++) {
        threads[k].join();
    }
    unsigned remaining_rows = H % NTHREADS;
    for (unsigned i = H-remaining_rows; i < H; i++) {
        for (unsigned j = 0; j < W; j++) {
            Vect dir = Vect(i + 0.5 - W/2.f, j + 0.5 - H/2.f, 55 - f) - Vect(0,0,55);
            double r_i = 0;
            double g_i = 0;
            double b_i = 0;
            bool break_early = false;
            for (int k = 0; k < NRAYS; k++) {
                Vect c = scene.colour(Ray(Vect(0,0,55), dir), NJUMPS, break_early);
                r_i += c.get_x();
                g_i += c.get_y();
                b_i += c.get_z();
                if (break_early) {
                    break;
                }
            }
            double sent_rays = NRAYS;
            if (break_early) {
                sent_rays = 1;
            }
            //print(r_color);
            image[(i*W + j) * 3 + 0] = (unsigned char) std::min(255, (int) pow((r_i/sent_rays), 1/GAMMA)); //R
            image[(i*W + j) * 3 + 1] = (unsigned char) std::min(255, (int) pow((g_i/sent_rays), 1/GAMMA)); //G
            image[(i*W + j) * 3 + 2] = (unsigned char) std::min(255, (int) pow((b_i/sent_rays), 1/GAMMA)); //B
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    time_t t1 = time(0) - t0;
    std::cout << "Finished in " << t1 << " seconds" << std::endl;
    return 0;
}