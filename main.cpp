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

#define W 512 //Width of camera sensor in pixels
#define H 512 //Height of camera sensor in pixels
#define FOV 120 //In degrees

void print(Vect v) {
    std::cout << "(" << v.get_x() << " " << v.get_y() << " " << v.get_z() << ")" << std::endl;
}

int main() {
    Ball b = Ball(Vect(0,1,-10), 1);
    Scene scene;
    //Background balls
    scene.add_ball(Ball(Vect(0,0,-1000), 900, Vect(0, 255, 255)));
    scene.add_ball(Ball(Vect(0,0,1000), 900, Vect(0, 255, 0)));
    scene.add_ball(Ball(Vect(0,-1000,0), 900, Vect(0, 0, 255)));
    scene.add_ball(Ball(Vect(0,1000,0), 900, Vect(255, 255, 0)));
    scene.add_ball(Ball(Vect(-1000,0,0), 900, Vect(255, 0, 255)));
    scene.add_ball(Ball(Vect(1000,0,0), 900, Vect(255, 0, 0)));
    //Target ball(s)
    //scene.add_ball(Ball(Vect(0,0,-10), 1));
    //Light it up!
    scene.add_light(Light(Vect(20, 0, 0), 100));
    double f = (W/2.f)/tan(FOV/2.f);
    //Render the scene from the point of view of a camera at (0,0,0)
    //With the conventions in the course notes
    std::vector<unsigned char> image(W*H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vect dir = Vect(i + 0.5 - W/2.f, j + 0.5 - H/2.f, - f) - Vect(0,0,0);
            Ray r1(Vect(0,0,0), dir);
            Vect r_color = scene.colour(r1);
            //print(r_color);
            image[(i*W + j) * 3 + 0] = r_color.get_x(); //R
            image[(i*W + j) * 3 + 1] = r_color.get_y(); //G
            image[(i*W + j) * 3 + 2] = r_color.get_z(); //B
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    return 0;
}