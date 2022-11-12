#include <iostream>
#include <thread>
#include <cstring>

#include "parser.h"
#include "ppm.h"
#include "vectorize.h"
#include "utils.h"
#include "render.h"

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
        Vec3f* result = new Vec3f[image_width * image_height];
        const int number_of_cores = thread::hardware_concurrency();
        
        if (number_of_cores == 0 || image_height < number_of_cores)
            render_scene(q, u, v, pixel_width, pixel_height, image_width, 0, image_height, scene, cam, result);

        else {
            thread* threads = new thread[number_of_cores];
            const int height_increase = image_height / number_of_cores;
            
            for (int i = 0; i < number_of_cores; i++) {
                const int min_height = i * height_increase;
                const int max_height = (i == number_of_cores - 1) ? image_height : (i + 1) * height_increase;
                threads[i] = thread(&render_scene, q, u, v, pixel_width, pixel_height, image_width, min_height, max_height, scene, cam, result);
            }
        
            for (int i = 0; i < number_of_cores; i++) threads[i].join();
            delete[] threads;
        }
        
        int k = 0;
        for (int i = 0; i < image_width; i++) {
            for (int j = 0; j < image_height; j++) {
                const Vec3f pixel = result[i * image_height + j];
                image[k++] = pixel.x;
                image[k++] = pixel.y;
                image[k++] = pixel.z;
            }
        }
        delete[] result;
        char image_name[cam.image_name.size() + 1];
        strcpy(image_name,(cam.image_name).c_str());
        write_ppm(image_name, image, image_width, image_height); 
        delete[] image;
    } // Camera
    
}
