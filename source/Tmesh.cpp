#include <string>
#include <string.h>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <list>
#include <math.h>
#include "../headers/Tmesh.hpp"
#include "../headers/Vect.hpp"
#include "../headers/Geometry.hpp"


TriangleIndices::TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1, bool added = false) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk), group(group) {
}

bbox_tree::~bbox_tree() {
    free(left_child);
    free(right_child);
}

TriangleMesh::TriangleMesh() : Geometry(Vect(1,1,1)) {
    bboxes = NULL;
}

TriangleMesh::TriangleMesh(const Vect& ind1, const Vect& ind2, const Vect& ind3) : Geometry(Vect(1,1,1)) {
    vertices.push_back(ind1);
    vertices.push_back(ind2);
    vertices.push_back(ind3);
    normals.push_back(Vect(0,0,0));
    normals.push_back(Vect(0,0,0));
    normals.push_back(Vect(0,0,0));
    uvs.push_back(Vect(0,0,0));
    uvs.push_back(Vect(0,0,0));
    uvs.push_back(Vect(0,0,0));
    vertexcolors.push_back(Vect(0,0,0));
    vertexcolors.push_back(Vect(0,0,0));
    vertexcolors.push_back(Vect(0,0,0));
    indices.push_back(TriangleIndices(0, 1, 2, 0, 1, 2, 0, 1, 2, 0, true));
    bboxes = NULL;
}

TriangleMesh::TriangleMesh(const char* obj) : Geometry(Vect(1,1,1)) {
    readOBJ(obj);
    compute_bboxes();
}

TriangleMesh::~TriangleMesh() {
    free(bboxes);
}

bbox TriangleMesh::compute_bbox(std::vector<TriangleIndices>::const_iterator begin, std::vector<TriangleIndices>::const_iterator end) {
    double min_x = vertices[begin->vtxi].get_x();
    double max_x = vertices[begin->vtxi].get_x();
    double min_y = vertices[begin->vtxi].get_y();
    double max_y = vertices[begin->vtxi].get_y();
    double min_z = vertices[begin->vtxi].get_z();
    double max_z = vertices[begin->vtxi].get_z();
    for (std::vector<TriangleIndices>::const_iterator it = begin; it < end; it++) {
        if (vertices[it->vtxi].get_x() < min_x) { min_x = vertices[it->vtxi].get_x(); }
        if (vertices[it->vtxj].get_x() < min_x) { min_x = vertices[it->vtxj].get_x(); }
        if (vertices[it->vtxk].get_x() < min_x) { min_x = vertices[it->vtxk].get_x(); }
        if (vertices[it->vtxi].get_x() > max_x) { max_x = vertices[it->vtxi].get_x(); }
        if (vertices[it->vtxj].get_x() > max_x) { max_x = vertices[it->vtxj].get_x(); }
        if (vertices[it->vtxk].get_x() > max_x) { max_x = vertices[it->vtxk].get_x(); }
        if (vertices[it->vtxi].get_y() < min_y) { min_y = vertices[it->vtxi].get_y(); }
        if (vertices[it->vtxj].get_y() < min_y) { min_y = vertices[it->vtxj].get_y(); }
        if (vertices[it->vtxk].get_y() < min_y) { min_y = vertices[it->vtxk].get_y(); }
        if (vertices[it->vtxi].get_y() > max_y) { max_y = vertices[it->vtxi].get_y(); }
        if (vertices[it->vtxj].get_y() > max_y) { max_y = vertices[it->vtxj].get_y(); }
        if (vertices[it->vtxk].get_y() > max_y) { max_y = vertices[it->vtxk].get_y(); }
        if (vertices[it->vtxi].get_z() < min_z) { min_z = vertices[it->vtxi].get_z(); }
        if (vertices[it->vtxj].get_z() < min_z) { min_z = vertices[it->vtxj].get_z(); }
        if (vertices[it->vtxk].get_z() < min_z) { min_z = vertices[it->vtxk].get_z(); }
        if (vertices[it->vtxi].get_z() > max_z) { max_z = vertices[it->vtxi].get_z(); }
        if (vertices[it->vtxj].get_z() > max_z) { max_z = vertices[it->vtxj].get_z(); }
        if (vertices[it->vtxk].get_z() > max_z) { max_z = vertices[it->vtxk].get_z(); }
    }
    Vect b_min = Vect(min_x, min_y, min_z);
    Vect b_max = Vect(max_x, max_y, max_z);
    return bbox(b_min, b_max, begin, end);
}

bbox_tree* TriangleMesh::compute_bboxes_aux(std::vector<TriangleIndices>::iterator begin, std::vector<TriangleIndices>::iterator end) {
    bbox_tree* node = new bbox_tree(compute_bbox(begin, end), NULL, NULL);
    Vect diag = node->box.b_max - node->box.b_min;
    Vect mid_diag = node->box.b_min + diag*0.5;
    unsigned char longest_axis;
    if (abs(diag.get_x()) >= abs(diag.get_y()) && abs(diag.get_x()) >= abs(diag.get_z())) {
        longest_axis = 0;
    } else if (abs(diag.get_y()) >= abs(diag.get_x()) && abs(diag.get_y()) >= abs(diag.get_z())) {
        longest_axis = 1;
    } else {
        longest_axis = 2;
    }
    std::sort(begin, end, [&](TriangleIndices t1, TriangleIndices t2) -> bool {
        Vect barycenter_t1 = (vertices[t1.vtxi] + vertices[t1.vtxj] + vertices[t1.vtxk])/(1.f/3);
        Vect barycenter_t2 = (vertices[t2.vtxi] + vertices[t2.vtxj] + vertices[t2.vtxk])/(1.f/3);
        return barycenter_t1[longest_axis] <= barycenter_t2[longest_axis];
    });
    //std::vector<TriangleIndices>::iterator pivot = begin;
    //for (std::vector<TriangleIndices>::iterator i = begin; i < end; i++) {
    //    Vect barycenter = (vertices[i->vtxi] + vertices[i->vtxj] + vertices[i->vtxk])/(1.f/3);
    //    if (barycenter[longest_axis] < mid_diag[longest_axis]) {
    //        std::swap(*i, *pivot);
    //        pivot++;
    //    }
    //}
    std::vector<TriangleIndices>::iterator midpoint = begin + std::distance(begin,end)/2;
    if (std::distance(begin, end) <= 500) {
        return node;
    }
    //std::cout << std::distance(begin, pivot) << std::endl;
    node->left_child = compute_bboxes_aux(begin, midpoint);
    node->right_child = compute_bboxes_aux(midpoint, end);
    return node;
}

void TriangleMesh::compute_bboxes() {
    if (bboxes != NULL) {
        free(bboxes);
    }
    bboxes = compute_bboxes_aux(indices.begin(), indices.end());
}

double TriangleMesh::bbox_intersect(const bbox& box, const Ray& ray) const {
    const double xi0 = (box.b_min - ray.get_start()).dot(Vect(1,0,0))/(ray.get_dir().dot(Vect(1,0,0)));
    const double xi1 = (box.b_max - ray.get_start()).dot(Vect(1,0,0))/(ray.get_dir().dot(Vect(1,0,0)));
    if (xi0 < 0 && xi1 < 0) { return false; }
    const double yi0 = (box.b_min - ray.get_start()).dot(Vect(0,1,0))/(ray.get_dir().dot(Vect(0,1,0)));
    const double yi1 = (box.b_max - ray.get_start()).dot(Vect(0,1,0))/(ray.get_dir().dot(Vect(0,1,0)));
    if (yi0 < 0 && yi1 < 0) { return false; }
    const double zi0 = (box.b_min - ray.get_start()).dot(Vect(0,0,1))/(ray.get_dir().dot(Vect(0,0,1)));
    const double zi1 = (box.b_max - ray.get_start()).dot(Vect(0,0,1))/(ray.get_dir().dot(Vect(0,0,1)));
    if (zi0 < 0 && zi1 < 0) { return false; }
    const double tx0 = std::min(xi0, xi1);
    const double tx1 = std::max(xi0, xi1);
    const double ty0 = std::min(yi0, yi1);
    const double ty1 = std::max(yi0, yi1);
    const double tz0 = std::min(zi0, zi1);
    const double tz1 = std::max(zi0, zi1);
    const double pos_i = std::max(std::max(tx0, ty0), tz0);
    if (pos_i <= std::min(std::min(tx1, ty1), tz1)) {
        return pos_i;
    }
    return -1;
}

double TriangleMesh::get_k0_in() const {
    return 1;
}

double TriangleMesh::get_k0_out() const {
    return 1;
}

double TriangleMesh::get_ri_in() const {
    return -1;
}

double TriangleMesh::get_ri_out() const {
    return 1;
}

bool TriangleMesh::is_mirror() const {
    return false;
}

bool TriangleMesh::is_transparent() const {
    return false;
}

void TriangleMesh::scale(const double& scalar) {
    for (unsigned i = 0; i < vertices.size(); i++) {
        vertices[i] = vertices[i]*scalar;
    }
}

void TriangleMesh::transpose(const Vect& vector) {
    for (unsigned i = 0; i < vertices.size(); i++) {
        vertices[i] += vector;
    }
}

void TriangleMesh::print_bbox_tree() {
    print_aux(bboxes);
}

void TriangleMesh::print_aux(bbox_tree* node) {
    if ((node->left_child == NULL) && (node->right_child == NULL)) {
        std::cout << std::distance(node->box.begin, node->box.end) << std::endl;
    }
    if (node->left_child != NULL) {
        print_aux(node->left_child);
    }
    if (node->right_child != NULL) {
        print_aux(node->right_child);
    }
}

Intersection TriangleMesh::intersection(const Ray& ray) const {
    if (bbox_intersect(bboxes->box, ray) == -1) {
        return Intersection(Vect(0,0,0), Vect(0,0,0), false);
    }
    std::list<bbox_tree*> nodes_to_visit = {};
    nodes_to_visit.push_front(bboxes);
    double smallest_t = INFINITY;
    std::vector<TriangleIndices>::const_iterator t_it = indices.begin();
    double t_alpha = 0;
    double t_beta = 0;
    double t_gamma = 0;
    while(!nodes_to_visit.empty()) {
        bbox_tree* cur_node = nodes_to_visit.back();
        nodes_to_visit.pop_back();
        bool hit_left = false;
        bool hit_right = false;
        if (cur_node->left_child != NULL) {
            double pos_i = bbox_intersect(cur_node->left_child->box, ray);
            if (pos_i >= 0) {
                if (pos_i < smallest_t) {
                    hit_left = true;
                    nodes_to_visit.push_back(cur_node->left_child);
                }
            }
        }
        if (cur_node->right_child != NULL) {
            double pos_i = bbox_intersect(cur_node->right_child->box, ray);
            if (pos_i >= 0) {
                if (pos_i < smallest_t) {
                    hit_right = true;
                    nodes_to_visit.push_back(cur_node->right_child);
                }
            }
        }
        if (!hit_right && !hit_left) {
            //std::cout << std::distance(cur_node->box.begin, cur_node->box.end) << std::endl;
            for (std::vector<TriangleIndices>::const_iterator i = cur_node->box.begin; i < cur_node->box.end; i++) {
                const TriangleIndices triangle = *i;
                const Vect e1 = vertices[triangle.vtxj] - vertices[triangle.vtxi];
                const Vect e2 = vertices[triangle.vtxk] - vertices[triangle.vtxi];
                const Vect normal = e1.cross(e2);
                const double udN = ray.get_dir().dot(normal);
                const Vect AmO = vertices[triangle.vtxi] - ray.get_start();
                const double t = AmO.dot(normal)/udN;
                if (t > 0 && t < smallest_t) {
                    const Vect AmOcrossDir = AmO.cross(ray.get_dir());
                    const double beta = (e2.dot(AmOcrossDir))/udN;
                    if (beta < 0) { continue; }
                    const double gamma = -(e1.dot(AmOcrossDir))/udN;
                    if (gamma < 0) { continue; }
                    const double alpha = 1 - beta - gamma;
                    if (alpha < 0) { continue; }
                    smallest_t = t;
                    t_it = i;
                    t_alpha = alpha;
                    t_beta = beta;
                    t_gamma = gamma;
                }
            }
        }
    }
    if (smallest_t == INFINITY) {
        return Intersection(Vect(0,0,0), Vect(0,0,0), false);
    }
    Vect g_normal = (normals[t_it->ni]*t_alpha + normals[t_it->nj]*t_beta + normals[t_it->nk]*t_gamma).normalize();
    return Intersection(ray.get_start() + ray.get_dir()*smallest_t, g_normal, true);
}

void TriangleMesh::readOBJ(const char* obj) {
	char matfile[255];
	char grp[255];
	FILE* f;
	f = fopen(obj, "r");
	int curGroup = -1;
	while (!feof(f)) {
		char line[255];
		if (!fgets(line, 255, f)) break;
		std::string linetrim(line);
		linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
		strcpy(line, linetrim.c_str());
		if (line[0] == 'u' && line[1] == 's') {
			sscanf(line, "usemtl %[^\n]\n", grp);
			curGroup++;
		}
		if (line[0] == 'v' && line[1] == ' ') {
			Vect vec;
			Vect col;
			if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[1], &vec[2], &col[0], &col[1], &col[2]) == 6) {
				col[0] = std::min(1., std::max(0., col[0]));
				col[1] = std::min(1., std::max(0., col[1]));
				col[2] = std::min(1., std::max(0., col[2]));
				vertices.push_back(vec);
				vertexcolors.push_back(col);
			} else {
				sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
				vertices.push_back(vec);
			}
		}
		if (line[0] == 'v' && line[1] == 'n') {
			Vect vec;
			sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[1], &vec[2]);
			normals.push_back(vec);
		}
		if (line[0] == 'v' && line[1] == 't') {
			Vect vec;
			sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
			uvs.push_back(vec);
		}
		if (line[0] == 'f') {
			TriangleIndices t;
			int i0, i1, i2, i3;
			int j0, j1, j2, j3;
			int k0, k1, k2, k3;
			int nn;
			t.group = curGroup;
			char* consumedline = line + 1;
			int offset;
			nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
			if (nn == 9) {
				if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
				if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
				if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
				if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
				if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
				if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
				if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
				if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
				if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
				indices.push_back(t);
			} else {
				nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
				if (nn == 6) {
					if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
					if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
					if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
					if (j0 < 0) t.uvi = uvs.size() + j0; else	t.uvi = j0 - 1;
					if (j1 < 0) t.uvj = uvs.size() + j1; else	t.uvj = j1 - 1;
					if (j2 < 0) t.uvk = uvs.size() + j2; else	t.uvk = j2 - 1;
					indices.push_back(t);
				} else {
					nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
					if (nn == 3) {
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						indices.push_back(t);
					} else {
						nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
						if (i0 < 0) t.vtxi = vertices.size() + i0; else	t.vtxi = i0 - 1;
						if (i1 < 0) t.vtxj = vertices.size() + i1; else	t.vtxj = i1 - 1;
						if (i2 < 0) t.vtxk = vertices.size() + i2; else	t.vtxk = i2 - 1;
						if (k0 < 0) t.ni = normals.size() + k0; else	t.ni = k0 - 1;
						if (k1 < 0) t.nj = normals.size() + k1; else	t.nj = k1 - 1;
						if (k2 < 0) t.nk = normals.size() + k2; else	t.nk = k2 - 1;
						indices.push_back(t);
					}
				}
			}
			consumedline = consumedline + offset;
			while (true) {
				if (consumedline[0] == '\n') break;
				if (consumedline[0] == '\0') break;
				nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
				TriangleIndices t2;
				t2.group = curGroup;
				if (nn == 3) {
					if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
					if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
					if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
					if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
					if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
					if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
					if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
					if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
					if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;
					indices.push_back(t2);
					consumedline = consumedline + offset;
					i2 = i3;
					j2 = j3;
					k2 = k3;
				} else {
					nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
					if (nn == 2) {
						if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
						if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
						if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
						if (j0 < 0) t2.uvi = uvs.size() + j0; else	t2.uvi = j0 - 1;
						if (j2 < 0) t2.uvj = uvs.size() + j2; else	t2.uvj = j2 - 1;
						if (j3 < 0) t2.uvk = uvs.size() + j3; else	t2.uvk = j3 - 1;
						consumedline = consumedline + offset;
						i2 = i3;
						j2 = j3;
						indices.push_back(t2);
					} else {
						nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
						if (nn == 2) {
							if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
							if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
							if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
							if (k0 < 0) t2.ni = normals.size() + k0; else	t2.ni = k0 - 1;
							if (k2 < 0) t2.nj = normals.size() + k2; else	t2.nj = k2 - 1;
							if (k3 < 0) t2.nk = normals.size() + k3; else	t2.nk = k3 - 1;								
							consumedline = consumedline + offset;
							i2 = i3;
							k2 = k3;
							indices.push_back(t2);
						} else {
							nn = sscanf(consumedline, "%u%n", &i3, &offset);
							if (nn == 1) {
								if (i0 < 0) t2.vtxi = vertices.size() + i0; else	t2.vtxi = i0 - 1;
								if (i2 < 0) t2.vtxj = vertices.size() + i2; else	t2.vtxj = i2 - 1;
								if (i3 < 0) t2.vtxk = vertices.size() + i3; else	t2.vtxk = i3 - 1;
								consumedline = consumedline + offset;
								i2 = i3;
								indices.push_back(t2);
							} else {
								consumedline = consumedline + 1;
							}
						}
					}
				}
			}
		}
	}
	fclose(f);
}