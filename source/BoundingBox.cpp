#include "../headers/BoundingBox.hpp"
#include "../headers/Ray.hpp"
#include "../headers/Vect.hpp"
#include "../headers/Tmesh.hpp"
#include <algorithm>

BoundingBox::BoundingBox(TriangleMesh mesh) {
    min_point = *std::min_element(mesh.vertices.begin(), mesh.vertices.end(),
                [&](Vect v1, Vect v2) -> bool { 
                    return v1.get_x() <= v2.get_x() && v1.get_y() <= v2.get_y() && v1.get_z() <= v2.get_z(); });
    max_point = *std::max_element(mesh.vertices.begin(), mesh.vertices.end(),
                [&](Vect v1, Vect v2) -> bool { 
                    return v1.get_x() <= v2.get_x() && v1.get_y() <= v2.get_y() && v1.get_z() <= v2.get_z(); });
}

bool BoundingBox::intersects(const Ray& ray) {
    const double tx0 = (min_point - ray.get_start()).dot(Vect(1,0,0))/(ray.get_dir().dot(Vect(1,0,0)));
    const double tx1 = (max_point - ray.get_start()).dot(Vect(1,0,0))/(ray.get_dir().dot(Vect(1,0,0)));
    if (tx0 < 0 && tx1 < 0) { return false; }
    const double ty0 = (min_point - ray.get_start()).dot(Vect(0,1,0))/(ray.get_dir().dot(Vect(0,1,0)));
    const double ty1 = (max_point - ray.get_start()).dot(Vect(0,1,0))/(ray.get_dir().dot(Vect(0,1,0)));
    if (ty0 < 0 && ty1 < 0) { return false; }
    const double tz0 = (min_point - ray.get_start()).dot(Vect(0,0,1))/(ray.get_dir().dot(Vect(0,0,1)));
    const double tz1 = (max_point - ray.get_start()).dot(Vect(0,0,1))/(ray.get_dir().dot(Vect(0,0,1)));
    if (tz0 < 0 && tz1 < 0) { return false; }
    if (std::max(std::max(std::min(tx0, tx1), std::min(ty0, ty1)), std::min(tz0, tz1)) <= std::min(std::min(std::max(tx0, tx1), std::max(ty0, ty1)), std::max(tz0, tz1))) {
        return true;
    }
    return false;
}
