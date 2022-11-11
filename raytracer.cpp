#include <iostream>
#include "parser.h"
#include "ppm.h"
#include "vectorize.h"
#include "utils.h"
#include <cstring>

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

        l = cam.near_plane.x; r = cam.near_plane.y; 
        b = cam.near_plane.z; t = cam.near_plane.w;
        w = -cam.gaze; v = unit_vector(cam.up); u = unit_vector(cross(v, w));
        m = cam.position + cam.gaze * cam.near_distance;
        q = m + l*u + t*v;
        
        // Precompute pixel width and pixel height
        float pixel_width, pixel_height;
        pixel_width  = (r - l) / image_width;
        pixel_height = (t - b) / image_height;

        // Render 
        int k = 0; 
        for (int i = 0; i < image_height; i++){
            cout << i << '/' << image_height << endl;
            for (int j = 0; j < image_width; j++){
                // Generate the ray
                Ray ray = generate_ray(j, i, q, u, v, pixel_width, pixel_height, cam);  
                int recursion_depth = scene.max_recursion_depth;
                bool flag = true;
                Vec3f color = make_life_more_colorful(scene, cam, ray, flag, recursion_depth);
                // cout << "Final Color: " <<  (int) (color.x+0.5) << ' ' << (int) (color.y+0.5) << ' ' << (int) (color.z+0.5) << endl;
                image[k++] = (int) (color.x+0.5);
                image[k++] = (int) (color.y+0.5);
                image[k++] = (int) (color.z+0.5);
            
            } // j
        }  // i
        cout << endl;
        cout << endl;
        cout << endl;
        char image_name[cam.image_name.size() + 1];
        strcpy(image_name,(cam.image_name).c_str());
        write_ppm(image_name, image, image_width, image_height); 

    } // Camera
    
}
