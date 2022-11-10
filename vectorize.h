#ifndef __vectorize_h__
#define __vectorize_h__

#include "parser.h"
#include <cmath>

using namespace parser;
using namespace std;

inline ostream& operator<<(ostream &out, const Vec3f &v) {
    return out << v.x << ' ' << v.y << ' ' << v.z;
}

inline Vec3f operator+(const Vec3f &u, const Vec3f &v){
    return Vec3f(u.x + v.x,  u.y + v.y, u.z + v.z);
}

inline Vec3f operator-(const Vec3f &u, const Vec3f &v){
    return Vec3f(u.x - v.x,  u.y - v.y, u.z - v.z);
}

inline Vec3f operator-(const Vec3f &u){
    return Vec3f(-u.x, -u.y, -u.z);
}

inline Vec3f operator*(const Vec3f &u, const Vec3f &v){
    return Vec3f(u.x * v.x,  u.y * v.y, u.z * v.z);
}

inline Vec3f operator*(float t, const Vec3f &v){
    return Vec3f(t * v.x,  t * v.y, t * v.z);
}

inline Vec3f operator*(const Vec3f &v, float t){
    return t * v;
}

inline Vec3f operator/(const Vec3f &v, float t){
    return (1/t) * v;
}

inline float dot(const Vec3f &u, const Vec3f &v){
    return u.x * v.x
        +  u.y * v.y
        +  u.z * v.z; 
}

inline Vec3f cross(const Vec3f &u, const Vec3f &v) {
    return Vec3f(u.y * v.z - u.z * v.y,
                u.z * v.x - u.x * v.z,
                u.x * v.y - u.y * v.x);
}

inline float length(Vec3f v) {
    return sqrt(dot(v, v));
}

inline Vec3f unit_vector(Vec3f v) {
    return v / length(v);
}

inline Vec3f at(Ray r, float t) {
    return r.origin + t * unit_vector(r.direction);
}

#endif 
