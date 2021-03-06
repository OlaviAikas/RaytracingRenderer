#include "headers/Vect.hpp"
#include "headers/Ray.hpp"
#include "headers/Ball.hpp"
#include "headers/Geometry.hpp"
#include "headers/Tmesh.hpp"
#include "headers/Scene.hpp"
#include <iostream>
#include <vector>
#include <math.h>
#include <thread>
#include <algorithm>
#include <random> 

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define W 1024 //Width of camera sensor in pixels
#define H 1024 //Height of camera sensor in pixels
#define FOV 1.1 //In radians
#define NRAYS 1 //How many rays are checked per pixel
#define NJUMPS 5 //How many times a ray can bounce
#define GAMMA 0.45
#define NTHREADS 8

void print(Vect v) {
    std::cout << "(" << v.get_x() << " " << v.get_y() << " " << v.get_z() << ")" << std::endl;
}

void render_rows(const Scene& s, const double& f, const unsigned& thread_n, std::vector<unsigned char>& image) {
    std::default_random_engine engine(time(0)+thread_n);
    std::uniform_real_distribution<double> distr(0, 1);
    for (unsigned i = thread_n; i < H; i += NTHREADS) {
        for (unsigned j = 0; j < W; j++) {
            double r_i = 0;
            double g_i = 0;
            double b_i = 0;
            unsigned sent_rays = 0;
            while (sent_rays < NRAYS) {
                double r1 = distr(engine);
                double r2 = distr(engine);
                double bmx = sqrt(-2 * log(r1))*cos(2*M_PI*r2);
                double bmy = sqrt(-2 * log(r1))*sin(2*M_PI*r2);
                if (abs(bmx) > 0.5 || abs(bmy) > 0.5) {
                    continue;
                }
                sent_rays++;
                Vect dir = Vect(j + 0.5 + bmx - W/2.f, -(i + 0.5 + bmy - H/2.f), 55 - f) - Vect(0,0,55);
                const Ray pr = Ray(Vect(0,0,55), dir);
                bool db = false;
                Vect c = s.colour(pr, NJUMPS, engine, distr);
                r_i += c.get_x();
                g_i += c.get_y();
                b_i += c.get_z();
            }
            //print(r_color);
            image[(i*W + j) * 3 + 0] = (unsigned char) std::min(255, (int) pow((r_i/NRAYS), GAMMA)); //R
            image[(i*W + j) * 3 + 1] = (unsigned char) std::min(255, (int) pow((g_i/NRAYS), GAMMA)); //G
            image[(i*W + j) * 3 + 2] = (unsigned char) std::min(255, (int) pow((b_i/NRAYS), GAMMA)); //B
        }
    }
    std::cout << "Thread " << thread_n << " finished" << std::endl;
}

int main() {
    Scene scene;
    srand(time(0));
    //Background balls
    scene.add_geo(new Ball(Vect(0,0,-1000), 940, Vect(0,1,0)));
    scene.add_geo(new Ball(Vect(0,0,1000), 940, Vect(1,0,1)));
    scene.add_geo(new Ball(Vect(-1000,0,0), 940, Vect(0, 1, 1)));
    scene.add_geo(new Ball(Vect(1000,0,0), 940, Vect(1, 1, 0)));
    scene.add_geo(new Ball(Vect(0,1000,0), 940, Vect(1, 0, 0)));
    scene.add_geo(new Ball(Vect(0,-950,0), 940, Vect(0, 0, 1)));
    //Target ball(s)
    //scene.add_geo(new Ball(Vect(0, 20, 0), 2));
    //scene.add_geo(new Ball(Vect(0,0,0), 10, true));
    //scene.add_geo(new Ball(Vect(0,20,0), 10, 1.5));
    //scene.add_geo(new Ball(Vect(0,20,0), 9.3, 1., 1.5));
    TriangleMesh* cat = new TriangleMesh("cat.obj");
    //std::cout << cat->normals.size() << std::endl;
    //std::cout << cat->vertices.size() << std::endl;
    //std::cout << cat->vertices.size() << std::endl;
    cat->scale(0.6);
    cat->transpose(Vect(0, -10, 0));
    cat->compute_bboxes();
    //cat->print_bbox_tree();
    scene.add_geo(cat);
    //Light it up!
    scene.add_light(Light(Vect(-10, 20, 40),20000000000));
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
    std::default_random_engine engine(time(0));
    std::uniform_real_distribution<double> distr(0, 1);
    for (unsigned i = H-remaining_rows; i < H; i++) {
        for (unsigned j = 0; j < W; j++) {
            double r_i = 0;
            double g_i = 0;
            double b_i = 0;
            unsigned sent_rays = 0;
            while (sent_rays < NRAYS) {
                double r1 = distr(engine);
                double r2 = distr(engine);
                double bmx = sqrt(-2 * log(r1))*cos(2*M_PI*r2);
                double bmy = sqrt(-2 * log(r1))*sin(2*M_PI*r2);
                if (abs(bmx) > 0.5 || abs(bmy) > 0.5) {
                    continue;
                }
                sent_rays++;
                Vect dir = Vect(j + 0.5 + bmx - W/2.f, -(i + 0.5 + bmy - H/2.f), 55 - f) - Vect(0,0,55);
                const Ray pr = Ray(Vect(0,0,55), dir);
                Vect c = scene.colour(pr, NJUMPS, engine, distr);
                r_i += c.get_x();
                g_i += c.get_y();
                b_i += c.get_z();
            }
            //print(r_color);
            image[(i*W + j) * 3 + 0] = (unsigned char) std::min(255, (int) pow((r_i/NRAYS), GAMMA)); //R
            image[(i*W + j) * 3 + 1] = (unsigned char) std::min(255, (int) pow((g_i/NRAYS), GAMMA)); //G
            image[(i*W + j) * 3 + 2] = (unsigned char) std::min(255, (int) pow((b_i/NRAYS), GAMMA)); //B
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    time_t t1 = time(0) - t0;
    std::cout << "Finished in " << t1 << " seconds" << std::endl;
    scene.free_geos();
    return 0;
}