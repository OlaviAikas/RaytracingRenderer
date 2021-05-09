#pragma once
#include "Vect.hpp"
#include "Geometry.hpp"
#include "Ray.hpp"
#include <vector>

class TriangleIndices {
public:
	TriangleIndices(int vtxi, int vtxj, int vtxk, int ni, int nj, int nk, int uvi, int uvj, int uvk, int group, bool added);
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;  // indices within the uv coordinates array
	int ni, nj, nk;  // indices within the normals array
	int group;       // face group
};

struct bbox {
    Vect b_min;
    Vect b_max;
    std::vector<TriangleIndices>::const_iterator begin;
    std::vector<TriangleIndices>::const_iterator end;
    bbox(Vect b_min, Vect b_max, std::vector<TriangleIndices>::const_iterator begin, std::vector<TriangleIndices>::const_iterator end) : b_max(b_max), b_min(b_min), begin(begin), end(end) {};
};

struct bbox_tree {
    bbox_tree* left_child;
    bbox_tree* right_child;
    bbox box;
    bbox_tree(bbox box, bbox_tree* lc, bbox_tree* rc) : left_child(lc), right_child(rc), box(box) {};
    ~bbox_tree();
};


class TriangleMesh : public Geometry {
    private:
        bbox_tree* bboxes;
        bbox compute_bbox(std::vector<TriangleIndices>::const_iterator begin, std::vector<TriangleIndices>::const_iterator end);
        bbox_tree* compute_bboxes_aux(std::vector<TriangleIndices>::iterator begin, std::vector<TriangleIndices>::iterator end);
        double bbox_intersect(const bbox& box, const Ray& ray) const;
        void print_aux(bbox_tree* node);

    public:
        ~TriangleMesh();
	    TriangleMesh();
        TriangleMesh(const Vect& ind1, const Vect& ind2, const Vect& ind3);
        TriangleMesh(const char* obj);

        //debugging
        void print_bbox_tree();
    
        std::vector<TriangleIndices> indices;
	    std::vector<Vect> vertices;
	    std::vector<Vect> normals;
	    std::vector<Vect> uvs;
	    std::vector<Vect> vertexcolors;

        Intersection intersection(const Ray& ray) const;
        bool is_mirror() const;
        bool is_transparent() const;
        double get_ri_in() const;
        double get_ri_out() const;
        double get_k0_in() const;
        double get_k0_out() const;

        void scale(const double& scalar);
        void transpose(const Vect& vector);

	    void readOBJ(const char* obj);
        void compute_bboxes();
};