// Compile: g++ -O3 -std=c++11 2d_triangles.cpp
#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <vector>

using namespace std;
using Real = float;

uniform_real_distribution<Real> uni_dist(0, 1);

// Some basic vector operations
template <typename T>
struct Vec2 {
    T x, y;
    Vec2(T x = 0, T y = 0) : x(x), y(y) {}
};
template <typename T>
struct Vec3 {
    T x, y, z;
    Vec3(T x = 0, T y = 0, T z = 0) : x(x), y(y), z(z) {}
};
using Vec2f = Vec2<Real>;
using Vec3i = Vec3<int>;
using Vec3f = Vec3<Real>;
Vec2f operator+(const Vec2f &v0, const Vec2f &v1) {return Vec2f{v0.x+v1.x, v0.y+v1.y};}
Vec2f& operator+=(Vec2f &v0, const Vec2f &v1) {v0.x += v1.x; v0.y += v1.y; return v0;}
Vec2f operator-(const Vec2f &v0, const Vec2f &v1) {return Vec2f{v0.x-v1.x, v0.y-v1.y};}
Vec2f operator*(Real s, const Vec2f &v) {return Vec2f{v.x * s, v.y * s};}
Vec2f operator*(const Vec2f &v, Real s) {return Vec2f{v.x * s, v.y * s};}
Vec2f operator/(const Vec2f &v, Real s) {return Vec2f{v.x / s, v.y / s};}
Real dot(const Vec2f &v0, const Vec2f &v1) {return v0.x * v1.x + v0.y * v1.y;}
Real length(const Vec2f &v) {return sqrt(dot(v, v));}
Vec2f normal(const Vec2f &v) {return Vec2f{-v.y, v.x};}
Vec3f& operator+=(Vec3f &v0, const Vec3f &v1) {v0.x += v1.x; v0.y += v1.y; v0.z += v1.z; return v0;}
Vec3f operator-(const Vec3f &v0, const Vec3f &v1) {return Vec3f{v0.x-v1.x, v0.y-v1.y, v0.z-v1.z};}
Vec3f operator-(const Vec3f &v) {return Vec3f{-v.x, -v.y, -v.z};}
Vec3f operator*(const Vec3f &v, Real s) {return Vec3f{v.x*s, v.y*s, v.z*s};}
Vec3f operator*(Real s, const Vec3f &v) {return Vec3f{v.x*s, v.y*s, v.z*s};}
Vec3f operator/(const Vec3f &v, Real s) {return Vec3f{v.x/s, v.y/s, v.z/s};}
Real dot(const Vec3f &v0, const Vec3f &v1) {return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;}

template <typename T>
T clamp(T v, T l, T u) {
    if (v < l) return l;
    else if (v > u) return u;
    return v;
}

struct TriangleMesh {
    vector<Vec2f> vertices;
    vector<Vec3i> indices;
    vector<Vec3f> colors;
};

struct DTriangleMesh {
    DTriangleMesh(int num_vertices, int num_colors) {
        vertices.resize(num_vertices, Vec2f{0, 0});
        colors.resize(num_colors, Vec3f{0, 0, 0});
    }

    vector<Vec2f> vertices;
    vector<Vec3f> colors;
};

struct Edge {
    int v0, v1;

    Edge(int v0, int v1) : v0(min(v0, v1)), v1(max(v0, v1)) {}

    bool operator<(const Edge &e) const {
        return this->v0 != e.v0 ? this->v0 < e.v0 : this->v1 < e.v1;
    }
};

struct Sampler {
    vector<Real> pmf, cdf;
};

Sampler build_edge_sampler(const TriangleMesh &mesh,
                           const vector<Edge> &edges) {
    vector<Real> pmf;
    vector<Real> cdf;
    pmf.reserve(edges.size());
    cdf.reserve(edges.size() + 1);
    cdf.push_back(0);
    for (auto edge : edges) {
        auto v0 = mesh.vertices[edge.v0];
        auto v1 = mesh.vertices[edge.v1];
        pmf.push_back(length(v1 - v0));
        cdf.push_back(pmf.back() + cdf.back());
    }
    auto length_sum = cdf.back();
    for_each(pmf.begin(), pmf.end(), [&](Real &p) {p /= length_sum;});
    for_each(cdf.begin(), cdf.end(), [&](Real &p) {p /= length_sum;});
    return Sampler{pmf, cdf};
}

// Binary search to invert the CDF in the sampler
int sample(const Sampler &sampler, const Real u) {
    return clamp<int>(upper_bound(
        sampler.cdf.begin(), sampler.cdf.end(), u) - sampler.cdf.begin() - 1,
        0, sampler.cdf.size() - 2);
}

// Given a triangle mesh, collect all edges.
vector<Edge> collect_edges(const TriangleMesh &mesh) {
    set<Edge> edges;
    for (auto index : mesh.indices) {
        edges.insert(Edge(index.x, index.y));
        edges.insert(Edge(index.y, index.z));
        edges.insert(Edge(index.z, index.x));
    }
    return vector<Edge>(edges.begin(), edges.end());
}

struct Img {
    Img(int width, int height, const Vec3f &val = Vec3f{0, 0, 0}) :
            width(width), height(height) {
        color.resize(width * height, val);
    }

    vector<Vec3f> color;
    int width;
    int height;
};

void save_img(const Img &img, const string &filename, bool flip = false) {
    fstream fs(filename.c_str(), fstream::out);
    fs << "P3" << endl << img.width << " " << img.height << " 255" << endl;
    for (int i = 0; i < img.width * img.height; i++) {
        auto tonemap = [](Real x) {
            return int(pow(clamp(x, Real(0), Real(1)), Real(1/2.2))*255 + Real(.5));};
        auto c = flip ? -img.color[i] : img.color[i];
        fs << tonemap(c.x) << " " << tonemap(c.y) << " " << tonemap(c.z) << " ";
    }
}

// Trace a single ray at screen_pos, intersect with the triangle mesh.
Vec3f raytrace(const TriangleMesh &mesh,
               const Vec2f &screen_pos,
               int *hit_index = nullptr) {
    // Loop over all triangles in a mesh, return the first one that hits
    for (int i = 0; i < (int)mesh.indices.size(); i++) {
        // Retrieve the three vertices of a triangle
        auto index = mesh.indices[i];
        auto v0 = mesh.vertices[index.x];
        auto v1 = mesh.vertices[index.y];
        auto v2 = mesh.vertices[index.z];
        // Form three half-planes: v1-v0, v2-v1, v0-v2
        // If a point is on the same side of all three half-planes, it's inside the triangle.
        auto n01 = normal(v1 - v0);
        auto n12 = normal(v2 - v1);
        auto n20 = normal(v0 - v2);
        auto side01 = dot(screen_pos - v0, n01) > 0;
        auto side12 = dot(screen_pos - v1, n12) > 0;
        auto side20 = dot(screen_pos - v2, n20) > 0;
        if ((side01 && side12 && side20) || (!side01 && !side12 && !side20)) {
            if (hit_index != nullptr) {
                *hit_index = i;
            }
            return mesh.colors[i];
        }
    }
    if (hit_index != nullptr) {
        *hit_index = -1;
    }
    return Vec3f{0, 0, 0};
}

void render(const TriangleMesh &mesh,
            int samples_per_pixel,
            mt19937 &rng,
            Img &img) {
    auto sqrt_num_samples = (int)sqrt((Real)samples_per_pixel);
    samples_per_pixel = sqrt_num_samples * sqrt_num_samples;
    for (int y = 0; y < img.height; y++) {
        for (int x = 0; x < img.width; x++) {
            for (int dy = 0; dy < sqrt_num_samples; dy++) {
                for (int dx = 0; dx < sqrt_num_samples; dx++) {
                    auto xoff = (dx + uni_dist(rng)) / sqrt_num_samples;
                    auto yoff = (dy + uni_dist(rng)) / sqrt_num_samples;
                    auto screen_pos = Vec2f{x + xoff, y + yoff};
                    auto color = raytrace(mesh, screen_pos);
                    img.color[y * img.width + x] += color / samples_per_pixel;
                }
            }
        }
    }
}

void compute_interior_derivatives(const TriangleMesh &mesh,
                                  int samples_per_pixel,
                                  const Img &adjoint,
                                  mt19937 &rng,
                                  vector<Vec3f> &d_colors) {
    auto sqrt_num_samples = (int)sqrt((Real)samples_per_pixel);
    samples_per_pixel = sqrt_num_samples * sqrt_num_samples;
    for (int y = 0; y < adjoint.height; y++) {
        for (int x = 0; x < adjoint.width; x++) {
            for (int dy = 0; dy < sqrt_num_samples; dy++) {
                for (int dx = 0; dx < sqrt_num_samples; dx++) {
                    auto xoff = (dx + uni_dist(rng)) / sqrt_num_samples;
                    auto yoff = (dy + uni_dist(rng)) / sqrt_num_samples;
                    auto screen_pos = Vec2f{x + xoff, y + yoff};
                    int hit_index = -1;
                    raytrace(mesh, screen_pos, &hit_index);
                    if (hit_index != -1) {
                        d_colors[hit_index] += 
                            adjoint.color[y * adjoint.width + x] / samples_per_pixel;
                    }
                }
            }
        }
    }
}

void compute_edge_derivatives(
        const TriangleMesh &mesh,
        const vector<Edge> &edges,
        const Sampler &edge_sampler,
        const Img &adjoint,
        const int num_edge_samples,
        mt19937 &rng,
        Img &screen_dx,
        Img &screen_dy,
        vector<Vec2f> &d_vertices) {
    for (int i = 0; i < num_edge_samples; i++) {
        // Pick an edge
        auto edge_id = sample(edge_sampler, uni_dist(rng));
        auto edge = edges[edge_id];
        auto pmf = edge_sampler.pmf[edge_id];
        // Pick a point p on the edge
        auto v0 = mesh.vertices[edge.v0];
        auto v1 = mesh.vertices[edge.v1];
        auto t = uni_dist(rng);
        auto p = v0 + t * (v1 - v0);
        auto xi = (int)p.x; auto yi = (int)p.y; // integer coordinates
        if (xi < 0 || yi < 0 || xi >= adjoint.width || yi >= adjoint.height) {
            continue;
        }
        // Sample the two sides of the edge
        auto n = normal((v1 - v0) / length(v1 - v0));
        auto color_in = raytrace(mesh, p - 1e-3f * n);
        auto color_out = raytrace(mesh, p + 1e-3f * n);
        // Get corresponding adjoint from the adjoint image,
        // multiply with the color difference and divide by the pdf & number of samples.
        auto weight = Real(length(v1 - v0) / (pmf * Real(num_edge_samples)));
        auto adj = dot(color_in - color_out, adjoint.color[yi * adjoint.width + xi]) * weight;
        // Reynolds transport theorem:
        // The boundary point is p = v0 + t * (v1 - v0)
        // The derivatives w.r.t. q is color_diff * dot(n, dp/dq)
        // For v0, dp/dv0.x = (1 - t, 0), dp/dv0.y = (0, 1 - t)
        // For v1, dp/dv1.x = (t, 0), dp/dv1.y = (0, t)
        auto d_v0 = Vec2f{(1 - t) * n.x, (1 - t) * n.y} * adj;
        auto d_v1 = Vec2f{     t  * n.x,      t  * n.y} * adj;
        // For the derivatives w.r.t. p, dp/dp.x = (1, 0) and dp/dp.y = (0, 1)
        // The screen space derivatives are the negation of this
        auto dx = -n.x * (color_in - color_out) * weight;
        auto dy = -n.y * (color_in - color_out) * weight;
        // In the parallel case, use atomic add here.
        screen_dx.color[yi * adjoint.width + xi] += dx;
        screen_dy.color[yi * adjoint.width + xi] += dy;
        d_vertices[edge.v0] += d_v0;
        d_vertices[edge.v1] += d_v1;
    }    
}

void d_render(const TriangleMesh &mesh,
              const Img &adjoint,
              const int interior_samples_per_pixel,
              const int edge_samples_in_total,
              mt19937 &rng,
              Img &screen_dx,
              Img &screen_dy,
              DTriangleMesh &d_mesh) {
    compute_interior_derivatives(mesh, interior_samples_per_pixel, adjoint, rng, d_mesh.colors);
    auto edges = collect_edges(mesh);
    auto edge_sampler = build_edge_sampler(mesh, edges);
    compute_edge_derivatives(mesh, edges, edge_sampler, adjoint, edge_samples_in_total,
                             rng, screen_dx, screen_dy, d_mesh.vertices);
}

int main(int argc, char *argv[]) {
    // Let's define two 2D triangles as a single triangle mesh
    TriangleMesh mesh{
        {{50.0, 25.0}, {200.0, 200.0}, {15.0, 150.0},
         {200.0, 15.0}, {150.0, 250.0}, {50.0, 100.0}},
        {{0, 1, 2}, {3, 4, 5}},
        {{0.3, 0.5, 0.3}, {0.3, 0.3, 0.5}}
    };
    Img img(256, 256);
    mt19937 rng(1234);
    render(mesh, 4 /* samples_per_pixel */, rng, img);
    save_img(img, "render.ppm");

    Img adjoint(img.width, img.height, Vec3f{1, 1, 1});
    Img dx(img.width, img.height), dy(img.width, img.height);
    DTriangleMesh d_mesh(mesh.vertices.size(), mesh.colors.size());
    d_render(mesh, adjoint, 4 /* interior_samples_per_pixel */,
             img.width * img.height /* edge_samples_in_total */, rng, dx, dy, d_mesh);
    save_img(dx, "dx_pos.ppm", false /*flip*/); save_img(dx, "dx_neg.ppm", true /*flip*/);
    save_img(dy, "dy_pos.ppm", false /*flip*/); save_img(dy, "dy_neg.ppm", true /*flip*/);
    return 0;
}
