#ifndef __utils_h__
#define __utils_h__

#include "parser.h"
#include <cmath>

using namespace parser;
using namespace std;
using Point3f = Vec3f;


Ray generate_ray(int i, 
                 int j, 
                 const Point3f &q, 
                 const Vec3f &u, 
                 const Vec3f &v, 
                 const float pixel_width, 
                 const float pixel_height,
                 const Camera &cam) {
    
    // Calculate the center point
    Point3f s;
    float s_u, s_v;
    s_u = (i + 0.5) * pixel_width;
    s_v = (j + 0.5) * pixel_height;
    s = q + s_u * unit_vector(u) - s_v * unit_vector(v);

    // Generate the ray
    Ray ray = Ray(cam.position, s - cam.position); 
    return ray;
}

float hit_sphere(const Point3f &center, float radius, const Ray &r) {
    Vec3f unit_direction = unit_vector(r.direction);
    Vec3f o_minus_c = r.origin - center;

    // Calculate A, B and C 
    float a = dot(unit_direction, unit_direction);
    float b = 2 * dot(unit_direction, o_minus_c);
    float c = dot(o_minus_c, o_minus_c) - r*r;
    float discriminant = b*b - 4*a*c;

    // Find t
    if (discriminant < 0)
        return -1;

    else if (discriminant == 0)
        return -b / (2 * a);

    else 
    {
        float delta = sqrt(discriminant);
        float t_1 = (-B + delta) / (2 * A);
        float t_2 = (-B - delta) / (2 * A);
        float t = t_1 < t_2 : t_1 ? t_2;
        return t; 
    } 
}

#endif 