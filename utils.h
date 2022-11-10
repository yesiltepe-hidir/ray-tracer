#ifndef __utils_h__
#define __utils_h__

#include "parser.h"
#include <cmath>
#include <vector>

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
    float c = dot(o_minus_c, o_minus_c) - radius*radius;
    float discriminant = b*b - 4*a*c;

    // Find t
    if (discriminant < 0)
        return -1;

    else if (discriminant == 0)
        return -b / (2 * a);

    else 
    {
        float delta = sqrt(discriminant);
        float t_1 = (-b + delta) / (2 * a);
        float t_2 = (-b - delta) / (2 * a);
        float t = t_1 < t_2 ? t_1 : t_2;
        return t; 
    } 
}

float determinant(const Vec3f &A, const Vec3f &B, const Vec3f &C) {
    // 3 x 3 Determinant
    return A.x * (B.y * C.z - B.z * C.y) 
        +  A.y * (B.z * C.x - B.x * C.z)
        +  A.z * (B.x * C.y - B.y * C.x);
}

float hit_triangle(const Vec3f &A, 
                   const Vec3f &B, 
                   const Vec3f &C,  
                   const Vec3f &normal_vector, 
                   const Ray &r) {
    /* 
    Barycentric Coordinates (alpha, beta, gamma) solution 
    
    Constraint: 
        alpha + beta + gamma = 1.0
        0 < alpha < 1 
        0 < beta  < 1 
        0 < gamma < 1 
    
    Trick: Calculate beta and gamma only. alpha = 1 - beta - gamma.
    */

    // 1. If triangle is parallel to the ray, then there is no hitting.
    if (dot(r.direction, normal_vector) == 0)
        return -1;

    Vec3f unit_direction    =   unit_vector(r.direction);

    // Calculate determinants
    float determinant_A     =   determinant(A - B, A - C, unit_direction);
    float determinant_beta  =   determinant(A - r.origin, A - C, unit_direction);
    float determinant_gamma =   determinant(A - B, A - r.origin, unit_direction);
    float determinant_t     =   determinant(A - B, A - C, A - r.origin);

    // Get alpha betta gamma and t.
    float alpha, beta, gamma, t;
    beta = determinant_beta / determinant_A;
    gamma = determinant_gamma / determinant_A;
    alpha = 1 - beta - gamma;
    t = determinant_t / determinant_A;
    
    // Hitting conditions
    bool c1, c2, c3;
    c1 = (beta + gamma <= 1); c2 = (beta  >= 0); c3 = (gamma >= 0);

    if ((c1 && c2) && c3) return t;
    else  return -1;
}

#endif 
