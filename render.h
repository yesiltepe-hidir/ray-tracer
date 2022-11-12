#ifndef __render_h__
#define __render_h__

#include "utils.h"
#include "parser.h"
using namespace parser;

Vec3f render_pixel(int i, 
                  int j, 
                  const Point3f &q, 
                  const Vec3f &u, 
                  const Vec3f &v, 
                  const float pixel_width, 
                  const float pixel_height,
                  const Camera &cam,
                  const Scene &scene) {
    
    Ray ray = generate_ray(j, i, q, u, v, pixel_width, pixel_height, cam);  
    int recursion_depth = scene.max_recursion_depth;
    bool flag = true;
    Vec3f color = make_life_more_colorful(scene, cam, ray, flag, recursion_depth);
    return color;

}

void render_scene(const Point3f &q, 
                  const Vec3f &u, 
                  const Vec3f &v, 
                  const float pixel_width, 
                  const float pixel_height,
                  const int width,
                  const int row_start_idx,
                  const int row_end_idx,
                  const Scene &scene,
                  const Camera &cam,
                  Vec3f* result) {

    for (int i = row_start_idx; i < row_end_idx; i++) {
        for (int j = 0; j < width; j++) {
            result[i * width + j] = render_pixel(i, j, 
                                                 q, u, v, 
                                                 pixel_width, pixel_height, 
                                                 cam, 
                                                 scene);
        }
    }
}


#endif
