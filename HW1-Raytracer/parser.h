#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <string>
#include <vector>
#include <cmath>

namespace parser
{
    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    struct Vec3f
    {
        float x, y, z;

        Vec3f() = default;

        Vec3f(float a, float b, float c) : x{a} , y{b}, z{c} {};

        Vec3f(const Vec3f& rhs) : x{rhs.x}, y{rhs.y}, z{rhs.z} {};

        Vec3f operator*(float mult){
            return {x*mult, y*mult, z*mult};
        }

        Vec3f& operator=(const Vec3f& rhs){
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            return *this;
        }

        Vec3f& operator=(Vec3f&& rhs){
            x = rhs.x;
            y = rhs.y;
            z = rhs.z;
            return *this;
        }

        Vec3f(Vec3f&& a) : x{a.x}, y{a.y}, z{a.z} {};

        Vec3f elemwise_prod(const Vec3f& mult){
            return {x*mult.x, y*mult.y, z*mult.z};
        }

        Vec3f elemwise_prod(Vec3f&& mult){
            return {x*mult.x, y*mult.y, z*mult.z};
        }

        void operator*=(float mult){
            x *= mult;
            y *= mult;
            z *= mult;
        }


        Vec3f operator+(const Vec3f& rhs){
            return {x+rhs.x, y+rhs.y, z+rhs.z};
        }

        Vec3f operator+(Vec3f&& rhs){
            return {x+rhs.x, y+rhs.y, z+rhs.z};
        }


        void operator+=(const Vec3f& rhs){
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
        
        }

        void operator+=(Vec3f&& rhs){
            x += rhs.x;
            y += rhs.y;
            z += rhs.z;
        }

        Vec3f operator-(const Vec3f& rhs){
            return {x-rhs.x, y-rhs.y, z-rhs.z};
        }

        Vec3f operator-(Vec3f&& rhs){
            return {x-rhs.x, y-rhs.y, z-rhs.z};
        }

        Vec3f operator/(float mult){
            return {x/mult, y/mult, z/mult};
        }

        void operator/=(float mult){
            x /= mult;
            y /= mult;
            z /= mult;
        }

        float dot_prod(const Vec3f& vec){

            return x*vec.x + y*vec.y + z*vec.z;
        }

        float dot_prod(Vec3f&& vec){

            return x*vec.x + y*vec.y + z*vec.z;
        }

        Vec3f cross_prod(const Vec3f& vec){
            float x_comp = y * vec.z - z * vec.y;
            float y_comp = z * vec.x - x * vec.z;
            float z_comp = x * vec.y - y * vec.x;
            return {x_comp, y_comp, z_comp};
        } 

        Vec3f cross_prod(Vec3f&& vec){
            float x_comp = y * vec.z - z * vec.y;
            float y_comp = z * vec.x - x * vec.z;
            float z_comp = x * vec.y - y * vec.x;
            return {x_comp, y_comp, z_comp};
        }


        void normalize(float len){
            x /= len;
            y /= len;
            z /= len;
        }

        void normalize(){
            float len = sqrt(this->dot_prod(*this));
            (*this) /= len;
        }

    };

    struct Vec3i
    {
        int x, y, z;
    };

    struct Vec4f
    {
        float x, y, z, w;
    };

    struct Camera
    {
        Vec3f position;
        Vec3f gaze;
        Vec3f up;
        Vec4f near_plane;
        float near_distance;
        int image_width, image_height;
        std::string image_name;
    };

    struct PointLight
    {
        Vec3f position;
        Vec3f intensity;
    };

    struct Material
    {
        bool is_mirror;
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        float phong_exponent;
    };

    struct Face
    {
        int v0_id;
        int v1_id;
        int v2_id;
    };

    struct Mesh
    {
        int material_id;
        std::vector<Face> faces;
        std::vector<Vec3f> normals;
    };

    struct Triangle
    {
        int material_id;
        Face indices;
        Vec3f normal;
    };

    struct Sphere
    {
        int material_id;
        int center_vertex_id;
        float radius;
        float radius_sqd;
    };


    struct Ray{
        Vec3f orig;
        Vec3f dir;
    };

    struct Scene
    {
        //Data
        Vec3i background_color;
        float shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;
        Vec3f bg_color;
    
    
        //Functions
        void loadFromXml(const std::string &filepath);

        void main_render();

        void initialize_mesh_normals();

        void initialize_tri_normals();

        void initialize_radius_sqd();

        void findIntersection(Ray&, int&, int&, int&, float, float&);

        bool intersectSphere(Ray& ray, Vec3f& center, float radius_sqd, float, float& tmax);
        bool intersectTriBFC(Ray& ray, Vec3f& a, Vec3f& b, Vec3f& c, Vec3f& normal, float, float& tmax);

        bool shadowRayIntersection(Ray& ray, float tmax);
        bool intersectTriShadow(Ray&, Vec3f&, Vec3f&, Vec3f&, float, float);

        Vec3f computePixelColor(Ray& ray, float);
        Vec3f recursiveRayTrace(Ray& , int);

        inline static unsigned char my_clamp(float);
    };

    ////////////////


}

#endif
