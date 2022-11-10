#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "vectorize.h"
#include "utils.h"

using namespace parser;
using namespace std;
using Point3f = Vec3f;

#define COLOR_CHANNELS 3


int main(int argc, char* argv[])
{    
    parser::Scene scene;
    scene.loadFromXml(argv[1]);
    
    int image_width, image_height;
  
    unsigned char* image;
    
    // For each camera
    for(int cam_idx = 0; cam_idx < scene.cameras.size(); cam_idx++){
        
        // Take current camera
        Camera cam = scene.cameras[cam_idx];
        
        // Image width & height
        image_width  = cam.image_width;
        image_height = cam.image_height;

        // Create the image
        image = new unsigned char[image_width * image_height * COLOR_CHANNELS];

        // Precompute points m, q and unit vector u.
        Point3f m, q;
        Vec3f u, v, w;
        float l, r, b, t;

        l = cam.near_plane.x; t = cam.near_plane.w;
        w = -cam.gaze; v = unit_vector(cam.up); u = unit_vector(cross(v, w));
        m = cam.position + unit_vector(cam.gaze) * cam.near_distance;
        q = m + l*u + t*v;
        
        // Precompute pixel width and pixel height
        float pixel_width, pixel_height;
        pixel_width = (r - l) / image_width;
        pixel_height = (t - b) / image_height;

        // Render 
        int k = 0; 
        for (int i = 0; i < image_height; i++){
            for (int j = 0; j < image_width; j++){
                // Generate the ray
                Ray ray = generate_ray(i, j, q, u, v, pixel_width, pixel_height, cam);  

            
            } // j
        }  // i
    } // Camera
    
}
