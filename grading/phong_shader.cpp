#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"
#include <cmath>




vec3 Phong_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color;
    //determine the color
    vec3 diffuse, specular, reflection, c, ambient;

    ambient = world.ambient_color * world.ambient_intensity * this->color_ambient;

    color = ambient;

    for (unsigned i = 0; i < world.lights.size(); i++){
    vec3 light_direction = (world.lights.at(i)->position - intersection_point);
    reflection = 2.0 * dot(light_direction.normalized(), normal) * normal + -light_direction.normalized();

    if(world.enable_shadows) { //Shadows present, see if object in the way
    Ray shadow_ray (intersection_point, light_direction);
    Hit shadow_hit;
    Object* shadow_obj = world.Closest_Intersection(shadow_ray, shadow_hit);

    if((!shadow_obj) || shadow_hit.t > light_direction.magnitude()){

    diffuse = ((world.lights.at(i)->Emitted_Light(ray.direction))/light_direction.magnitude_squared()) * this->color_diffuse * std::max((dot(normal, light_direction.normalized())), 0.0);
    color += diffuse;

    c = -ray.direction;
    specular = ((world.lights.at(i)->Emitted_Light(ray.direction))/light_direction.magnitude_squared()) * this->color_specular * pow(std::max(dot(reflection.normalized(), c), 0.0), specular_power);
    color += specular;
    }

    }

    else{//No shadows, calculate diffuse and specular no matter what
    diffuse = ((world.lights.at(i)->Emitted_Light(ray.direction))/light_direction.magnitude_squared()) * this->color_diffuse * std::max((dot(normal, light_direction.normalized())), 0.0);
    color += diffuse;

    reflection = 2.0 * dot(light_direction.normalized(), normal) * normal + -light_direction.normalized();
    c = -ray.direction;
    specular = ((world.lights.at(i)->Emitted_Light(ray.direction))/light_direction.magnitude_squared()) * this->color_specular * pow(std::max(dot(reflection.normalized(), c), 0.0), specular_power);
    color += specular;
    }

    }

    return color;

  }
