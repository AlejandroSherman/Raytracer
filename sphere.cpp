#include "sphere.h"
#include "ray.h"

// Determine if the ray intersects with the sphere
void Sphere::Intersection(const Ray& ray, std::vector<Hit>& hits) const
{
  double discrim = pow(dot(ray.direction, ray.endpoint - center), 2) - (dot(ray.direction, ray.direction) * (dot(ray.endpoint - center, ray.endpoint - center) - pow(radius, 2)));
  double t1 = 0;
  double t2 = 0;
  Hit hit;

  if (discrim > 0) {
      t1 = -dot(ray.direction, ray.endpoint - center) + sqrt(discrim) / dot(ray.direction, ray.direction);
      t2 = -dot(ray.direction, ray.endpoint - center) - sqrt(discrim) / dot(ray.direction, ray.direction);

      if (t1 < t2 && t1 >= small_t) {
        hit.t = t1;
        hit.object = this;
        hit.ray_exiting = true;
        hits.push_back(hit);
        hit.t = t2;
        hit.ray_exiting = false;
        hits.push_back(hit);
      }

      else if (t1 >= t2 && t2 >= small_t) {
        hit.t = t2;
        hit.object = this;
        hit.ray_exiting = false;
        hits.push_back(hit);
        hit.t = t1;
        hit.ray_exiting = true;
        hits.push_back(hit);
      }

      else if (t1 < t2 && t1 < small_t) {
        hit.t = 0;
        hit.object = this;
        hit.ray_exiting = false;
        hits.push_back(hit);
        hit.t = t1;
        hit.ray_exiting = true;
        hits.push_back(hit);
      }

      else if (t1 >= t2 && t2 < small_t) {
        hit.t = 0;
        hit.object = this;
        hit.ray_exiting = false;
        hits.push_back(hit);
        hit.t = t2;
        hit.ray_exiting = true;
        hits.push_back(hit);
      }

      if (debug_pixel){
        std::cout << "intersect test with obj[index below]: hits = {";
        for (size_t j = 0; j < hits.size(); j++){
          const Object* obj = hits[j].object;
          std::cout << " { "<< obj << ", " << hits[j].t << ", " << hits[j].ray_exiting << " },";
        }
        std::cout << " } " << std::endl;
      }

  }
}

vec3 Sphere::Normal(const vec3& point) const
{
    vec3 normal;
    // compute the normal direction
    normal = (point - center).normalized();
    return normal;
}

// set this->bounding_box so that it contains the sphere as tightly as possible.
void Sphere::Update_Bounding_Box()
{
    //TODO;
    infinite_box=false;
}
