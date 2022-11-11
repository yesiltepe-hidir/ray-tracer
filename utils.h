#ifndef __utils_h__
#define __utils_h__

#include "parser.h"
#include <cmath>
#include <vector>

#define MAX 9999999999
#define EPSILON 0.0000000000001

using namespace parser;
using namespace std;
using Point3f = Vec3f;

float max(float a, float b) {
    return a > b ? a : b;
}

float min(float a, float b) {
    return a < b ? a : b;
}

Vec3f clip(Vec3f color_intensity) {
    return Vec3f(max(0, min(color_intensity.x, 255)),
                 max(0, min(color_intensity.y, 255)),
                 max(0, min(color_intensity.z, 255)));
}

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
    s = q + s_u*u - s_v*v;

    // Generate the ray
    Ray ray = Ray(cam.position, unit_vector(s - cam.position)); 
    return ray;
}

float hit_sphere(const Point3f &center, float radius, const Ray &r) {
    Vec3f o_minus_c = r.origin - center;

    // Calculate A, B and C 
    float a = dot(r.direction, r.direction);
    float b = 2 * dot(r.direction, o_minus_c);
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
        beta + gamma <_ 1
        0 <_ beta
        0 <_ gamma
        t_min <= t <= t_max
    
    Trick: Calculate beta and gamma only. alpha = 1 - beta - gamma.
    */

    // 1. If triangle is parallel to the ray, then there is no hitting.
    if (dot(r.direction, normal_vector) == 0)
        return -1;

    // Calculate determinants
    float determinant_A     =   determinant(A - B, A - C, r.direction);
    float determinant_beta  =   determinant(A - r.origin, A - C, r.direction);
    float determinant_gamma =   determinant(A - B, A - r.origin, r.direction);
    float determinant_t     =   determinant(A - B, A - C, A - r.origin);

    // Get alpha betta gamma and t.
    float alpha, beta, gamma, t;
    beta = determinant_beta / determinant_A;
    gamma = determinant_gamma / determinant_A;
    alpha = 1 - beta - gamma;
    t = determinant_t / determinant_A;
    
    // Hitting conditions
    bool c1, c2, c3, c4;
    c1 = (beta + gamma <= 1); c2 = (beta  >= 0); c3 = (gamma >= 0); c4 = (t >= 0); // TODO: Consider again c4.  

    if ((c1 && c2) && (c3 && c4)) return t;
    else  return -1;
}

Vec3f sphere_normal(const Vec3f &intersection_point, const Vec3f& center, float radius){
    return unit_vector((intersection_point - center) / radius);
}

Vec3f triangle_normal(const Vec3f &A, const Vec3f &B, const Vec3f &C){
    return unit_vector(cross(C - B, A - B));
}

Vec3f make_life_more_colorful(const Scene &scene, 
                              const Camera cam, 
                              const Ray r, 
                              bool &flag, 
                              int &recursion_depth) {
    
    // Initializations
    float min_intersection = MAX, t;
    int min_sphere_idx=-1, min_triangle_idx=-1, min_mesh_idx=-1, min_face_idx=-1;
    Vec3f color(0, 0, 0), rec_color(0, 0, 0);
    
    // Base condition
    if (flag == false && recursion_depth == 0)
    {
        // cout << "Base: " << color << endl;
        return color;
    }
        

    // Sphere
    Vec3f center;
    float radius;
    for (int i = 0; i < scene.spheres.size(); i++)
    {
        center = scene.vertex_data[scene.spheres[i].center_vertex_id-1];
        radius = scene.spheres[i].radius;
        t = hit_sphere(center, radius, r);
        if (t < min_intersection && t >= EPSILON) { min_sphere_idx = i; min_intersection = t; }
    }

    // Triangle
    Vec3f A, B, C, normal_vector;
    for (int i = 0; scene.triangles.size(); i++) 
    {
        A = scene.vertex_data[scene.triangles[i].indices.v0_id-1];
        B = scene.vertex_data[scene.triangles[i].indices.v1_id-1];
        C = scene.vertex_data[scene.triangles[i].indices.v2_id-1];
        normal_vector = triangle_normal(A, B, C);
        t = hit_triangle(A, B, C, normal_vector, r);
        if (t < min_intersection && t >= EPSILON) { min_triangle_idx = i; min_intersection = t; }
    }

    // Mesh
    for (int i = 0; i < scene.meshes.size(); i++)
    {
        for (int j = 0; j < scene.meshes[i].faces.size(); j++)
        {
            A = scene.vertex_data[scene.meshes[i].faces[j].v0_id-1];
            B = scene.vertex_data[scene.meshes[i].faces[j].v1_id-1];
            C = scene.vertex_data[scene.meshes[i].faces[j].v2_id-1];
            normal_vector = triangle_normal(A, B, C);
            t = hit_triangle(A, B, C, normal_vector, r);
            if (t < min_intersection && t >= EPSILON) { min_mesh_idx = i; min_face_idx = j; min_intersection = t; }
        }   
    }
    
    // No intersection
    if(min_intersection == MAX)
    {
        color = scene.background_color;
        // cout << "No intersection: " << color << endl;
        return color;
    }
        

    Vec3f intersection = r.origin + min_intersection*r.direction;
    Material material;

    // Intersection: Mesh
    if (min_mesh_idx > -1)
    {
        // cout << "Intersection: Mesh" << endl;
        A = scene.vertex_data[scene.meshes[min_mesh_idx].faces[min_face_idx].v0_id-1];
        B = scene.vertex_data[scene.meshes[min_mesh_idx].faces[min_face_idx].v1_id-1];
        C = scene.vertex_data[scene.meshes[min_mesh_idx].faces[min_face_idx].v2_id-1];
        normal_vector = triangle_normal(A, B, C);
        material = scene.materials[scene.meshes[min_mesh_idx].material_id-1];
    }

    // Intersection: Triangle
    else if (min_triangle_idx > -1)
    {
        // // cout << "Intersection: Triangle" << endl;
        A = scene.vertex_data[scene.triangles[min_triangle_idx].indices.v0_id-1];       
        B = scene.vertex_data[scene.triangles[min_triangle_idx].indices.v1_id-1];       
        C = scene.vertex_data[scene.triangles[min_triangle_idx].indices.v2_id-1];       
        normal_vector = triangle_normal(A, B, C);
        material = scene.materials[scene.triangles[min_triangle_idx].material_id-1];
    }

    // Intersection: Sphere
    else if (min_sphere_idx > -1)
    {
        // // cout << "Intersection: Sphere" << endl;
        center = scene.vertex_data[scene.spheres[min_sphere_idx].center_vertex_id-1];
        radius = scene.spheres[min_sphere_idx].radius;
        normal_vector = sphere_normal(intersection, center, radius);
        material = scene.materials[scene.spheres[min_sphere_idx].material_id-1];
    }

    // Shadow: Ambient
    Vec3f ambient = material.ambient * scene.ambient_light;
    color = color + ambient;

    // Shadow: Diffuse & Specularity
    bool under_shadow;
    Ray shadow_ray;
    Vec3f w, normal_shadow, w0, w0_w1, h;
    Vec3f intensity, diffuse(0, 0, 0), specularity(0, 0, 0);
    float squared_length_w, length_w, cos, cos_phong;

    // For each point light
    for (int i = 0; i < scene.point_lights.size(); i++)
    {
        // Initialize parameters
        min_sphere_idx = -1; min_triangle_idx = -1; min_mesh_idx = -1; min_face_idx = -1;
        under_shadow = false;

        // Vector between light and the intersection point. 
        w = scene.point_lights[i].position - intersection;
        squared_length_w = dot(w, w);
        length_w = sqrt(squared_length_w);
        w = unit_vector(w);

        // Shadow Ray
        shadow_ray.origin = intersection + normal_vector*scene.shadow_ray_epsilon;
        shadow_ray.direction = w;

        // Shadow: Sphere
        for(int i = 0; i < scene.spheres.size() && !under_shadow; i++){
            center = scene.vertex_data[scene.spheres[i].center_vertex_id-1];
            radius = scene.spheres[i].radius;
            t = hit_sphere(center, radius, r);
            if(t < length_w && t >= EPSILON )  
                under_shadow = true;
        }

        // Shadow: Triangle
        for(int i = 0; i < scene.triangles.size() && !under_shadow; i++){
            
            A = scene.vertex_data[scene.triangles[i].indices.v0_id-1];
            B = scene.vertex_data[scene.triangles[i].indices.v1_id-1];
            C = scene.vertex_data[scene.triangles[i].indices.v2_id-1];
            normal_shadow = triangle_normal(A, B, C);
            t = hit_triangle(A, B, C, normal_shadow, r);
            if(t < length_w && t >= EPSILON) 
                under_shadow = true;
        }

        // Shadow: Mesh
        for(int i = 0; i < scene.meshes.size() && !under_shadow; i++){
            for(int j = 0; j < scene.meshes[i].faces.size() && !under_shadow; j++){
                A = scene.vertex_data[scene.meshes[i].faces[j].v0_id-1];
                B = scene.vertex_data[scene.meshes[i].faces[j].v1_id-1];
                C = scene.vertex_data[scene.meshes[i].faces[j].v2_id-1];
                normal_shadow = triangle_normal(A, B, C);
                t = hit_triangle(A, B, C, normal_shadow, r);
                if(t < length_w && t >= EPSILON) 
                under_shadow = true;
            }
        }

        if(under_shadow) continue;

        cos =  max(0, dot(w, normal_vector));
        if(cos == 0) continue;
        
        // Intensity Arrangement
        intensity = scene.point_lights[i].intensity;
        intensity = intensity * (1.0 / squared_length_w);
        diffuse = diffuse + material.diffuse * (intensity * cos);

        if(flag == false)
            w0 = unit_vector(-(r.origin + r.direction));
        else
            w0 = unit_vector(cam.position - intersection);
        
        w0_w1 = unit_vector(w + w0);
        h = unit_vector(w0_w1 * sqrt(dot(w0_w1, w0_w1)));
        cos_phong = (float) pow(max(0, dot(normal_vector, h)),material.phong_exponent);
        specularity = specularity + material.specular * (intensity * cos_phong);
    }
 
    if(material.is_mirror == true){
        Ray* reflected_ray = new Ray;
        reflected_ray->direction =  unit_vector(-w0 + ((2 * normal_vector) * dot(normal_vector, w0)));
        reflected_ray->origin = intersection;
        flag = false;
        recursion_depth-=1;
        rec_color = material.mirror * make_life_more_colorful(scene,
                                                                cam,
                                                                *reflected_ray,
                                                                flag,
                                                                recursion_depth);
    }
    // Clip the color
    Vec3f temp_color = diffuse + specularity + rec_color + color;
    // cout << "Before clipping: " <<  temp_color << endl;
    color = clip(diffuse + specularity + rec_color + color);
    // cout << "After clipping: " << color << endl;
    
    return color;                     
}  

#endif 