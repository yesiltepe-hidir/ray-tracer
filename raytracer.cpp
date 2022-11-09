// PPM Example
#include <iostream>
#include <cmath>
using namespace std;

// Define Vector class

class vec3{
    public:
        // Default Constructor
        vec3() : e{0, 0, 0} {}
        
        // Copy Constructor
        vec3(double e0, double e1, double e2) : e{e0, e1, e2}{}

        // Getters for coordinates
        double x() const { return e[0]; }
        double y() const { return e[1]; }
        double z() const { return e[2]; }

        // Operators
        vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
        double operator[](int i) const { return e[i]; }
        double& operator[](int i) { return e[i]; }

        vec3& operator+=(const vec3& v){
            e[0] += v.e[0];
            e[1] += v.e[1];
            e[2] += v.e[2];
            return *this;
        }

        vec3& operator*=(const double t){
            e[0] *= t;
            e[1] *= t;
            e[2] *= t;
            return *this;
        }

        vec3& operator/=(const double t){
            return *this *= 1/t;
        }

        // Length of the vector
        double length() const {
            return sqrt(length_squared());
        }

        double length_squared() const {
            return e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
        }

    public:
        double e[3];
};       

// Aliases for vec3
using point3 = vec3; // 3D point
using color = vec3; //  RGB color


// Utilities for vec3
inline ostream& operator<<(ostream &out, const vec3 &v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}

inline vec3 operator+(const vec3 &u, const vec3 &v){
    return vec3(u.e[0] + v.e[0],  u.e[1] + v.e[1], u.e[2] + v.e[2]);
}

inline vec3 operator-(const vec3 &u, const vec3 &v){
    return vec3(u.e[0] - v.e[0],  u.e[1] - v.e[1], u.e[2] - v.e[2]);
}

inline vec3 operator*(const vec3 &u, const vec3 &v){
    return vec3(u.e[0] * v.e[0],  u.e[1] * v.e[1], u.e[2] * v.e[2]);
}

inline vec3 operator*(double t, const vec3 &v){
    return vec3(t * v.e[0],  t * v.e[1], t * v.e[2]);
}

inline vec3 operator*(const vec3 &v, double t){
    return t * v;
}

inline vec3 operator/(const vec3 &v, double t){
    return (1/t) * v;
}

inline double dot(const vec3 &u, const vec3 &v){
    return u.e[0] * v.e[0]
        +  u.e[1] * v.e[1]
        +  u.e[2] * v.e[2]; 
}

inline vec3 cross(const vec3 &u, const vec3 &v) {
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}

inline vec3 unit_vector(vec3 v) {
    return v / v.length();
}

// Color Utilities
void write_color(ostream &out, color pixel_color){
    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(255.999 * pixel_color.x()) << ' '
        << static_cast<int>(255.999 * pixel_color.y()) << ' '
        << static_cast<int>(255.999 * pixel_color.z()) << '\n';
}

// Ray 
class ray {
    public:
        // Default constructor
        ray() {}
        
        // Constructor
        ray(const point3 &origin, const vec3 &direction) : orig(origin), dir(direction) {}
        
        // Getters
        point3 origin() const { return orig; }
        vec3 direction() const { return dir; }

        // At time t
        point3 at(double t) const {
            return orig + t * dir;
        }
  
    public:
        point3 orig;
        vec3 dir;
};

// Hitting Sphere 
bool hit_sphere(const point3& center, double radius, const ray& r) {
    vec3 oc = r.direction() - center;
    auto a = dot(r.direction(), r.direction());
    auto b = 2.0 * dot(oc, r.direction());
    auto c = dot(oc, oc) - radius*radius;
    auto discriminant = b*b - 4*a*c;
    return (discriminant > 0);
}

// Ray color
color ray_color(const ray& r) {
    if (hit_sphere(point3(0,0,-1.0), 0.5, r))
        return color(1, 0, 0);
    
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}


int main(){

    // Image
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 400;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    // Camera
    auto viewport_height = 2.0;
    auto viewport_width = aspect_ratio * viewport_height;
    auto focal_length = 1.0;


    auto origin = point3(0, 0, 0);
    auto horizontal = vec3(viewport_width, 0, 0);
    auto vertical = vec3(0, viewport_height, 0);
    auto lower_left_corner = origin - horizontal/2 - vertical/2 - vec3(0, 0, focal_length);
    

    // Render
    cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    // Row
    for (int j = image_height-1; j >= 0; --j)
    {
        // Progress indicator
        cerr << "\rScanlines remaining: " << j << ' ' << '\n';
        // Column 
        for (int i = 0; i < image_width; ++i)
        {
            auto u = double(i) / (image_width-1);
            auto v = double(j) / (image_height-1);
            ray r(origin, lower_left_corner + u*horizontal + v*vertical - origin);
            color pixel_color = ray_color(r);
            write_color(cout, pixel_color);
        }
    }
    cerr << "\nDone.\n";
}
