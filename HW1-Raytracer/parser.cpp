#include "parser.h"
#include "tinyxml2.h"
#include <sstream>
#include <stdexcept>
#include <limits>
#include <cmath>
#include "ppm.h"

void parser::Scene::loadFromXml(const std::string &filepath)
{
    tinyxml2::XMLDocument file;
    std::stringstream stream;

    auto res = file.LoadFile(filepath.c_str());
    if (res)
    {
        throw std::runtime_error("Error: The xml file cannot be loaded.");
    }

    auto root = file.FirstChild();
    if (!root)
    {
        throw std::runtime_error("Error: Root is not found.");
    }

    //Get BackgroundColor
    auto element = root->FirstChildElement("BackgroundColor");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0 0 0" << std::endl;
    }
    stream >> background_color.x >> background_color.y >> background_color.z;

    //Get ShadowRayEpsilon
    element = root->FirstChildElement("ShadowRayEpsilon");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0.001" << std::endl;
    }
    stream >> shadow_ray_epsilon;

    //Get MaxRecursionDepth
    element = root->FirstChildElement("MaxRecursionDepth");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0" << std::endl;
    }
    stream >> max_recursion_depth;

    //Get Cameras
    element = root->FirstChildElement("Cameras");
    element = element->FirstChildElement("Camera");
    Camera camera;
    while (element)
    {
        auto child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Gaze");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Up");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("NearPlane");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("NearDistance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageResolution");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageName");
        stream << child->GetText() << std::endl;

        stream >> camera.position.x >> camera.position.y >> camera.position.z;
        stream >> camera.gaze.x >> camera.gaze.y >> camera.gaze.z;
        stream >> camera.up.x >> camera.up.y >> camera.up.z;
        stream >> camera.near_plane.x >> camera.near_plane.y >> camera.near_plane.z >> camera.near_plane.w;
        stream >> camera.near_distance;
        stream >> camera.image_width >> camera.image_height;
        stream >> camera.image_name;

        cameras.push_back(camera);
        element = element->NextSiblingElement("Camera");
    }

    //Get Lights
    element = root->FirstChildElement("Lights");
    auto child = element->FirstChildElement("AmbientLight");
    stream << child->GetText() << std::endl;
    stream >> ambient_light.x >> ambient_light.y >> ambient_light.z;
    element = element->FirstChildElement("PointLight");
    PointLight point_light;
    while (element)
    {
        child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Intensity");
        stream << child->GetText() << std::endl;

        stream >> point_light.position.x >> point_light.position.y >> point_light.position.z;
        stream >> point_light.intensity.x >> point_light.intensity.y >> point_light.intensity.z;

        point_lights.push_back(point_light);
        element = element->NextSiblingElement("PointLight");
    }

    //Get Materials
    element = root->FirstChildElement("Materials");
    element = element->FirstChildElement("Material");
    Material material;
    while (element)
    {
        material.is_mirror = (element->Attribute("type", "mirror") != NULL);

        child = element->FirstChildElement("AmbientReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("DiffuseReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("SpecularReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("MirrorReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("PhongExponent");
        stream << child->GetText() << std::endl;

        stream >> material.ambient.x >> material.ambient.y >> material.ambient.z;
        stream >> material.diffuse.x >> material.diffuse.y >> material.diffuse.z;
        stream >> material.specular.x >> material.specular.y >> material.specular.z;
        stream >> material.mirror.x >> material.mirror.y >> material.mirror.z;
        stream >> material.phong_exponent;

        materials.push_back(material);
        element = element->NextSiblingElement("Material");
    }

    //Get VertexData
    element = root->FirstChildElement("VertexData");
    stream << element->GetText() << std::endl;
    Vec3f vertex;
    while (!(stream >> vertex.x).eof())
    {
        stream >> vertex.y >> vertex.z;
        vertex_data.push_back(vertex);
    }
    stream.clear();

    //Get Meshes
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Mesh");
    Mesh mesh;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> mesh.material_id;

        child = element->FirstChildElement("Faces");
        stream << child->GetText() << std::endl;
        Face face;
        while (!(stream >> face.v0_id).eof())
        {
            stream >> face.v1_id >> face.v2_id;
            mesh.faces.push_back(face);
        }
        stream.clear();

        meshes.push_back(mesh);
        mesh.faces.clear();
        element = element->NextSiblingElement("Mesh");
    }
    stream.clear();

    //Get Triangles
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Triangle");
    Triangle triangle;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> triangle.material_id;

        child = element->FirstChildElement("Indices");
        stream << child->GetText() << std::endl;
        stream >> triangle.indices.v0_id >> triangle.indices.v1_id >> triangle.indices.v2_id;

        triangles.push_back(triangle);
        element = element->NextSiblingElement("Triangle");
    }

    //Get Spheres
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Sphere");
    Sphere sphere;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> sphere.material_id;

        child = element->FirstChildElement("Center");
        stream << child->GetText() << std::endl;
        stream >> sphere.center_vertex_id;

        child = element->FirstChildElement("Radius");
        stream << child->GetText() << std::endl;
        stream >> sphere.radius;

        spheres.push_back(sphere);
        element = element->NextSiblingElement("Sphere");
    }
}

void parser::Scene::initialize_mesh_normals(){
    for(int i = 0; i < meshes.size(); ++i){
        Mesh& curr_mesh = meshes[i];
        curr_mesh.normals.reserve(curr_mesh.faces.size());
        for(int j = 0; j < curr_mesh.faces.size(); ++j){
            Vec3f& a = vertex_data[curr_mesh.faces[j].v0_id - 1];
            Vec3f& b = vertex_data[curr_mesh.faces[j].v1_id - 1];
            Vec3f& c = vertex_data[curr_mesh.faces[j].v2_id - 1];
            Vec3f b_a = b - a;
            Vec3f c_a = c - a;
            Vec3f tri_normal = b_a.cross_prod(c_a);
            tri_normal.normalize();
            curr_mesh.normals.push_back(tri_normal);
        }
    }
}

void parser::Scene::initialize_tri_normals(){
    for(Triangle& tri : triangles){
        Vec3f& a = vertex_data[tri.indices.v0_id-1];
        Vec3f& b = vertex_data[tri.indices.v1_id-1];
        Vec3f& c = vertex_data[tri.indices.v2_id-1];
        Vec3f b_a = b-a;
        Vec3f c_a = c-a;
        tri.normal = b_a.cross_prod(c_a);
        tri.normal.normalize();
    }
}

void parser::Scene::initialize_radius_sqd(){
    for(Sphere& sphere : spheres){
        sphere.radius_sqd = sphere.radius*sphere.radius;
    }
}

bool parser::Scene::intersectSphere(Ray& ray, Vec3f& center, float radius_sqd, float tmin, float& tmax){
    Vec3f center_to_orig = ray.orig - center;
    float b = ray.dir.dot_prod(center_to_orig);
    float ac = center_to_orig.dot_prod(center_to_orig) - radius_sqd;
    float discriminant = b*b-ac;
    if(discriminant < 0){
        return false;
    }
    float d = sqrt(discriminant);
    float t = -b - d;
    if(t > tmin && t < tmax){
        tmax = t;
        return true;
    }
    return false; 
}


bool parser::Scene::intersectTriBFC(Ray& ray, Vec3f& a, Vec3f& b, Vec3f& c, Vec3f& normal, float tmin, float& tmax){

    if(ray.dir.dot_prod(normal) > 0){
        return false;
    }

    Vec3f b_a = b-a;
    Vec3f c_a = c-a;

    Vec3f dirxc_a = ray.dir.cross_prod(c_a);
    float a_det = b_a.dot_prod(dirxc_a);

    if(std::fabs(a_det) < 0.0000001f){
        return false;
    }


    float a_inverse = 1.f/a_det;

    Vec3f a_to_orig = ray.orig - a;
    float beta = a_to_orig.dot_prod(dirxc_a) * a_inverse;

    if(beta < 0.f || beta > 1.f){
        return false;
    }

    Vec3f atoxb_a = a_to_orig.cross_prod(b_a);
    float gamma = ray.dir.dot_prod(atoxb_a) * a_inverse;

    if( gamma < 0 || beta+gamma > 1.f){
        return false;
    }

    float t = c_a.dot_prod(atoxb_a) * a_inverse;
    if(t > tmin && t < tmax){
        tmax = t;
        return true;
    }

    return false;
    
}


bool parser::Scene::intersectTriShadow(Ray& ray, Vec3f& a, Vec3f& b, Vec3f& c, float tmin, float tmax){

    Vec3f b_a = b-a;
    Vec3f c_a = c-a;

    Vec3f dirxc_a = ray.dir.cross_prod(c_a);
    float a_det = b_a.dot_prod(dirxc_a);

    if(std::fabs(a_det) < 0.0000001f){
        return false;
    }

    float a_inverse = 1.f/a_det;

    Vec3f a_to_orig = ray.orig - a;
    float beta = a_to_orig.dot_prod(dirxc_a) * a_inverse;

    if(beta < 0.f || beta > 1.f){
        return false;
    }

    Vec3f atoxb_a = a_to_orig.cross_prod(b_a);
    float gamma = ray.dir.dot_prod(atoxb_a) * a_inverse;

    if( gamma < 0.f || beta+gamma > 1.f){
        return false;
    }

    float t = c_a.dot_prod(atoxb_a) * a_inverse;

    if( t > tmin && t < tmax ){
        return true;
    }

    return false;
    
}


void parser::Scene::findIntersection(Ray& ray, int& obj_type, int& obj_index, int& tri_index, float tmin, float& t_isect){
    for(int i = 0; i < spheres.size(); ++i){
        if(intersectSphere(ray, vertex_data[spheres[i].center_vertex_id - 1], spheres[i].radius_sqd, tmin, t_isect)){
            obj_type = 0;
            obj_index = i;
        }
    }
    for(int j = 0; j < triangles.size(); ++j){
        if(intersectTriBFC(ray, vertex_data[triangles[j].indices.v0_id - 1], vertex_data[triangles[j].indices.v1_id - 1], vertex_data[triangles[j].indices.v2_id - 1], triangles[j].normal, tmin, t_isect)){
            obj_type = 1;
            obj_index = j;
        }
    }
    for(int k = 0; k < meshes.size(); ++k){
        Mesh& curr_mesh = meshes[k];
        for(int i = 0; i < curr_mesh.faces.size(); ++i){
            if(intersectTriBFC(ray, vertex_data[curr_mesh.faces[i].v0_id - 1], vertex_data[curr_mesh.faces[i].v1_id - 1], vertex_data[curr_mesh.faces[i].v2_id - 1], curr_mesh.normals[i], tmin, t_isect)){
                obj_type = 2;
                obj_index = k;
                tri_index = i;
            }
        }
    }
}

bool parser::Scene::shadowRayIntersection(Ray& ray, float t_max){
    for(int i = 0; i < spheres.size(); ++i){
        if(intersectSphere(ray, vertex_data[spheres[i].center_vertex_id - 1], spheres[i].radius_sqd, shadow_ray_epsilon, t_max)){
            return true;
        }
    }

    for(int j = 0; j < triangles.size(); ++j){
        if(intersectTriShadow(ray, vertex_data[triangles[j].indices.v0_id - 1], vertex_data[triangles[j].indices.v1_id - 1], vertex_data[triangles[j].indices.v2_id - 1], shadow_ray_epsilon, t_max)){
            return true;
        }
    }

    for(int k = 0; k < meshes.size(); ++k){
        Mesh& curr_mesh = meshes[k];
        for(int i = 0; i < curr_mesh.faces.size(); ++i){
            if(intersectTriShadow(ray, vertex_data[curr_mesh.faces[i].v0_id - 1], vertex_data[curr_mesh.faces[i].v1_id - 1], vertex_data[curr_mesh.faces[i].v2_id - 1], shadow_ray_epsilon, t_max)){
                return true;
            }
        }
    }
    return false;
}


parser::Vec3f parser::Scene::recursiveRayTrace(Ray& ray, int rec_depth){
    if(rec_depth > max_recursion_depth){
        return {0,0,0};
    }
    int which_obj = -1;
    int obj_index = -1;
    int tri_index = -1;
    float t_isect = std::numeric_limits<float>::infinity();

    findIntersection(ray, which_obj, obj_index, tri_index, shadow_ray_epsilon, t_isect);


    
    if(which_obj == -1){
        return bg_color;
    }else{
        Material* mat_type;

        Vec3f surface_normal;
        Vec3f isect_point = ray.orig + ray.dir*t_isect;
        switch(which_obj){

        case 0:
            mat_type = &materials[spheres[obj_index].material_id -1];

            surface_normal = isect_point - vertex_data[spheres[obj_index].center_vertex_id - 1];
            surface_normal /= spheres[obj_index].radius;
            break;

        case 1:
            mat_type = &materials[triangles[obj_index].material_id -1];

            surface_normal = triangles[obj_index].normal;
            break;

        case 2:
            mat_type = &materials[meshes[obj_index].material_id -1];

            surface_normal = meshes[obj_index].normals[tri_index];
            break;

            default: break;
        }
        
        Vec3f ambi = ambient_light.elemwise_prod(mat_type->ambient);

        for(PointLight& light_source : point_lights){
            Vec3f to_light = light_source.position - isect_point;
            
            float distance_sqd = to_light.dot_prod(to_light);
            float distance = sqrt(distance_sqd);
            to_light /= distance;
            Ray shadow_ray{isect_point, to_light};
            if(shadowRayIntersection(shadow_ray, distance)){
                continue;
            }else{

                float cosine = to_light.dot_prod(surface_normal);
                if(cosine > 0){
                    Vec3f diffuse = mat_type->diffuse.elemwise_prod(light_source.intensity);
                    diffuse *= cosine;
                    diffuse /= distance_sqd;
                    ambi += diffuse;

                    Vec3f half_vector = to_light - ray.dir;
                    half_vector.normalize();
                    float half_v_cosine = half_vector.dot_prod(surface_normal);

                    if(half_v_cosine > 0){
                        Vec3f specular = mat_type->specular.elemwise_prod(light_source.intensity);
                        specular *= std::pow(half_v_cosine, mat_type->phong_exponent);
                        specular /= distance_sqd;
                        ambi += specular;
                    }
                }

            }
        }

        if(mat_type->is_mirror){
            float proj_coeff = surface_normal.dot_prod(ray.dir);
            Vec3f refl_ray = ray.dir - surface_normal * 2 * proj_coeff;
            Ray reflected{isect_point, refl_ray};
            ambi += mat_type->mirror.elemwise_prod(recursiveRayTrace(reflected, rec_depth+1));
        }

        return ambi;

    }        
}

parser::Vec3f parser::Scene::computePixelColor(Ray& cam_ray, float dist){
    int which_obj = -1;
    int obj_index = -1;
    int tri_index = -1;
    float t_isect = std::numeric_limits<float>::infinity();
    findIntersection(cam_ray, which_obj, obj_index, tri_index, dist, t_isect);


    if(which_obj == -1){
        return bg_color;
    }else{
        Material* mat_type;
        Vec3f surface_normal;
        Vec3f isect_point = cam_ray.orig + cam_ray.dir*t_isect;
        switch(which_obj){

            case 0:
                mat_type = &materials[spheres[obj_index].material_id-1];

                surface_normal = isect_point - vertex_data[spheres[obj_index].center_vertex_id - 1];
                surface_normal /= spheres[obj_index].radius;
            break;

            case 1:
                mat_type = &materials[triangles[obj_index].material_id-1];
                surface_normal = triangles[obj_index].normal;
            break;

            case 2:
                mat_type = &materials[meshes[obj_index].material_id-1];
                surface_normal = meshes[obj_index].normals[tri_index];
            break;
    
            default: break;
        }

        Vec3f ambi = ambient_light.elemwise_prod(mat_type->ambient);

        for(PointLight& light_source : point_lights){
            Vec3f to_light = light_source.position - isect_point;
            float distance_sqd = to_light.dot_prod(to_light);
            float distance = sqrt(distance_sqd);
            to_light /= distance;            
            Ray isect_light{isect_point, to_light};
            if(shadowRayIntersection(isect_light, distance)){
                continue;
            }else{

                float cosine = to_light.dot_prod(surface_normal);
                if(cosine > 0){
                    Vec3f diffuse = mat_type->diffuse.elemwise_prod(light_source.intensity);
                    diffuse *= cosine;
                    diffuse /= distance_sqd;
                    ambi += diffuse;

                    Vec3f half_vector = to_light - cam_ray.dir;
                    half_vector.normalize();
                    float half_v_cosine = half_vector.dot_prod(surface_normal);

                    if(half_v_cosine > 0){
                        Vec3f specular = mat_type->specular.elemwise_prod(light_source.intensity);
                        specular *= std::pow(half_v_cosine, mat_type->phong_exponent);
                        specular /= distance_sqd;
                        ambi += specular;
                    }
                }

            }
        }

        if(mat_type->is_mirror){
            float proj_coeff = surface_normal.dot_prod(cam_ray.dir);
            Vec3f refl_ray = cam_ray.dir - surface_normal * 2 * proj_coeff;
            Ray reflected{isect_point, refl_ray};
            ambi += mat_type->mirror.elemwise_prod(recursiveRayTrace(reflected, 1));
        }

        return ambi;

    }
}

unsigned char parser::Scene::my_clamp(float a){            
    return a > 255.f ? 255 : (unsigned char) (a+0.5);
}

void parser::Scene::main_render(){
    using namespace parser;
    initialize_tri_normals();
    initialize_mesh_normals();
    initialize_radius_sqd();
    bg_color.x = background_color.x;
    bg_color.y = background_color.y;
    bg_color.z = background_color.z;
    for(Camera& curr_cam : cameras){
        unsigned char* pixels = new unsigned char[curr_cam.image_width*curr_cam.image_height*3];
        curr_cam.gaze.normalize();
        curr_cam.up.normalize();
        Vec3f u_vec = curr_cam.gaze.cross_prod(curr_cam.up);
        float pixel_width = (curr_cam.near_plane.y - curr_cam.near_plane.x) / curr_cam.image_width;
        float pixel_height = (curr_cam.near_plane.w - curr_cam.near_plane.z) / curr_cam.image_height;

        Vec3f top_left = curr_cam.position + (curr_cam.gaze * curr_cam.near_distance);
        top_left += u_vec * curr_cam.near_plane.x;
        top_left += curr_cam.up * curr_cam.near_plane.w;

        Ray pixel_ray{curr_cam.position, top_left};
        for(int j = 0; j < curr_cam.image_height; j++){
            float up_dist = (j+0.5) * pixel_height;
            size_t up_to_row_pix = j*curr_cam.image_width;
            for(int i = 0; i < curr_cam.image_width; i++){
                float u_vec_dist = (i+0.5) * pixel_width;
                Vec3f pix_center = top_left + (u_vec*u_vec_dist) - (curr_cam.up*up_dist);
                pixel_ray.dir = pix_center-curr_cam.position;
                float distance = sqrt(pixel_ray.dir.dot_prod(pixel_ray.dir));
                pixel_ray.dir /= distance;
                Vec3f pixel_color = computePixelColor(pixel_ray, distance);
                size_t pixel_loc = 3*(up_to_row_pix+i);

                pixels[pixel_loc]   = my_clamp(pixel_color.x);
                pixels[pixel_loc+1] = my_clamp(pixel_color.y);
                pixels[pixel_loc+2] = my_clamp(pixel_color.z);

            }
        }
        write_ppm(curr_cam.image_name.c_str(), pixels, curr_cam.image_width, curr_cam.image_height);
        delete [] pixels;
    };
}

