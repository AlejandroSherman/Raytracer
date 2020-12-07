#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"


Render_World::Render_World()
    :background_shader(0),ambient_intensity(0),enable_shadows(true),
    recursion_depth_limit(3)
{}

Render_World::~Render_World()
{
    delete background_shader;
    for(size_t i=0;i<objects.size();i++) delete objects[i];
    for(size_t i=0;i<lights.size();i++) delete lights[i];
}

// Find the closest object of intersection and return the object that was
// intersected.  Record the Hit structure in hit.  If no intersection occurred,
// return NULL.  Note that in the case of a Boolean, the object returned will be
// the Boolean, but the object stored in hit will be the underlying primitive.
// Any intersection with t<=small_t should be ignored.
Object* Render_World::Closest_Intersection(const Ray& ray, Hit& hit)
{
    double min_t = std::numeric_limits<double>::max();

    Hit closest_int;
    closest_int.object = NULL;
    closest_int.t = 0;
    std::vector<Hit> hits;
    Object* obj_to_return = NULL;
    unsigned save_i;

    for(unsigned i = 0; i < objects.size(); ++i){
      objects[i]->Intersection(ray, hits);
      for(unsigned j = 0; j < hits.size(); ++j){
        Hit curr_hit = hits[j];
        if (curr_hit.object != NULL){
          if (curr_hit.t > small_t && curr_hit.t < min_t){
            min_t = curr_hit.t;
            closest_int = curr_hit;
            hit = closest_int;
            obj_to_return = objects[i];
            save_i = i;
          }
        }
      }
    }
    //Object* obj_to_return = closest_int.object;

    if (debug_pixel){
      std::cout << "closest intersection: return obj[" << save_i << "]; hit =";
      const Object* obj = hit.object;
      std::cout << " { "<< obj << ", " << hit.t << ", " << hit.ray_exiting << " },";
      std::cout << " } " << std::endl;
    }

    return obj_to_return;
}

// set up the initial view ray and call
void Render_World::Render_Pixel(const ivec2& pixel_index)
{
    // set up the initial view ray here
    Ray ray;

    ray.endpoint = camera.position;
    vec3 direction = (camera.World_Position(pixel_index) - camera.position).normalized();
    ray.direction = direction;

    vec3 color=Cast_Ray(ray,1);
    camera.Set_Pixel(pixel_index,Pixel_Color(color));
}

void Render_World::Render()
{
    for(int j=0;j<camera.number_pixels[1];j++)
        for(int i=0;i<camera.number_pixels[0];i++)
            Render_Pixel(ivec2(i,j));
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray,int recursion_depth)
{
    vec3 color;
    Hit hit;
    Object* obj = Closest_Intersection(ray, hit);

    if (obj != NULL){
      vec3 intersection_point = ray.Point(hit.t);
      color = obj->material_shader->Shade_Surface(ray, intersection_point, obj->Normal(intersection_point), recursion_depth);
    }

    else{
      vec3 temp1;
      vec3 temp2;
      color = background_shader->Shade_Surface(ray, temp1, temp2, recursion_depth);
    }

    if (debug_pixel){
      std::cout << "cast ray: " << "end = " << ray.endpoint << std::endl;
      std::cout << "dir = " << ray.direction << std::endl;
      std::cout << "color = " << color << std::endl;
      std::cout << "call Shade_Surface with: " << std::endl;
      vec3 intersection_point = ray.Point(hit.t);
      std::cout << "location = " << intersection_point[0] << " " << intersection_point[1] << " " << intersection_point[2] << "; ";
      vec3 normal_test = obj->Normal(intersection_point);
      std::cout << "normal = " << normal_test[0] << " " << normal_test[1] << " " << normal_test[2] << std::endl;
    }

    return color;
}
