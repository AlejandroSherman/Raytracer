#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color = shader->Shade_Surface(ray, intersection_point, normal, recursion_depth);

    vec3 view_ray = -ray.direction;
  	vec3 reflection_direction = 2.0 * dot(view_ray.normalized(), normal) * normal - view_ray.normalized();
    Ray relection_ray(intersection_point, reflection_direction);

    if(recursion_depth < world.recursion_depth_limit){
       color = (1 - reflectivity) * color + reflectivity * world.Cast_Ray(relection_ray, recursion_depth + 1);
    }
    else
        color = (1 - reflectivity) * color;

    return color;
}
