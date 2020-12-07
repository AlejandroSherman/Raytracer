#include "plane.h"
#include "ray.h"
#include <cfloat>

// Intersect with the half space defined by the plane.  The plane's normal
// points outside.  If the ray starts on the "inside" side of the plane, be sure
// to record a hit with t=0 as the first entry in hits.
void Plane::
Intersection(const Ray& ray, std::vector<Hit>& hits) const
{
  double t;
  double denom;
  Hit hit;

  denom = dot(ray.direction, normal);

  if (denom != 0) {
  t = dot(x1 - ray.endpoint, normal) / denom;

  if (t < small_t) {
    hit.t = 0;
    hit.object = this;
    hits.push_back(hit);
    hit.t = t;
    hits.push_back(hit);
  }

  else{
    hit.t = t;
    hit.object = this;
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

vec3 Plane::
Normal(const vec3& point) const
{
    return normal;
}

void Plane::Update_Bounding_Box()
{
    bounding_box.Make_Full();
    infinite_box=true;
}
