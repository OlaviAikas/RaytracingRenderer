#include "Vect.hpp"
#include "Tmesh.hpp"
#include "Ray.hpp"

class BoundingBox {
    private:
        Vect min_point;
        Vect max_point;

    public:
        BoundingBox(TriangleMesh mesh);

        bool intersects(const Ray& ray);
};