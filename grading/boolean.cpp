#include "boolean.h"

// Determine if the ray intersects with the boolean of A and B.
void Boolean::Intersection(const Ray& ray, std::vector<Hit>& hits) const
{
    //TODO;
    // std::vector<Hit> a_hits;
    // A->Intersection(ray, a_hits);
    //
    // std::vector<Hit> b_hits;
    // B->Intersection(ray, b_hits);
    //
    // //for (unsigned ai = 0; ai < a_hits.length(); ++ai){
    //   //for (unsigned bi = 0; bi < b_hits.length(); ++bi){
    //
    //   //}
    // //}
    //
    // hits.push_back(a_hits[0]);
    // hits.push_back(b_hits[0]);


}

// This should never be called.
vec3 Boolean::Normal(const vec3& point) const
{
    assert(false);
    return vec3();
}

void Boolean::Update_Bounding_Box()
{
    A->Update_Bounding_Box();
    B->Update_Bounding_Box();
    // Compute bounding_box from A->bounding_box and B->bounding_box.
    // Note that this should depend on the type of boolean being performed.
    //TODO;
}
