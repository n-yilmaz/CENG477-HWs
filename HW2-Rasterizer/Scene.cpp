#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"

using namespace tinyxml2;
using namespace std;


void Scene::initialize_tform_mat(Mesh* mesh){
    Matrix4 last_form = id_matrix();
    int curr_id = 0;
    while(curr_id < mesh->numberOfTransformations){
        Matrix4 curr_tform;
        switch(mesh->transformationTypes[curr_id]){
            case 't':
                curr_tform = makeTMat(translations[mesh->transformationIds[curr_id]-1]);
                break;
            case 'r':
                curr_tform = makeRMat(rotations[mesh->transformationIds[curr_id]-1]);
                break;
            case 's':
                curr_tform = makeSMat(scalings[mesh->transformationIds[curr_id]-1]);
                break;
            default:
                curr_tform = id_matrix();
                break;
        }
        last_form = matxmat_mult(curr_tform, last_form);
        curr_id++;
    }
    mesh->tform_mat = last_form;
}

bool Scene::visible(double num, double denom, double& t_e, double& t_l){

    double t;
    if(denom < 0){
        t = num/denom;
        if(t > t_l){
            return false;
        }
        t_e = t > t_e ? t : t_e;
    }else if(denom > 0){
        t = num/denom;
        if(t < t_e){
            return false;
        }
        t_l = t < t_l ? t : t_l;
        
    }else if(denom == 0){
        if(num < 0){
            return false;
        }
        return true;
    }
    return true;

}

bool Scene::clipLineOrtho(Vec4& v1, Vec4& v2, Color& c1, Color& c2){
    double t_enter = 0.0;
    double t_exit = 1.0;
    double dx = v2.x-v1.x;
    double dy = v2.y-v1.y;
    double dz = v2.z-v1.z;
    
    if(visible(v1.x+1.0, -dx, t_enter, t_exit)){
        if(visible(1.0-v1.x, dx, t_enter, t_exit)){
            if(visible(v1.y+1.0, -dy, t_enter, t_exit)){
                if(visible(1.0-v1.y, dy, t_enter, t_exit)){
                    if(visible(v1.z+1.0, -dz, t_enter, t_exit)){
                        if(visible(1.0-v1.z, dz, t_enter, t_exit)){
                            Vec4 v1_new = v1;
                            Vec4 v2_new = v2;
                            c1 = *colorsOfVertices[v1.colorId];
                            c2 = *colorsOfVertices[v2.colorId];
                            Vec4 diff{dx, dy, dz, 0.0, -1};
                            Color c_diff = c2-c1;
                            if(t_exit < 1.0){
                                v2_new = v1 + (diff * t_exit);
                                c2     = *colorsOfVertices[v1.colorId] + (c_diff * t_exit);
                            }
                            if(t_enter > 0.0){
                                v1_new = v1 + (diff * t_enter);
                                c1     = *colorsOfVertices[v1.colorId] + (c_diff * t_enter);
                            }
                            v1 = v1_new;
                            v2 = v2_new;
                            return true;
                        }
                    }
                }
            }
        }
    }

    return false;
    
        

}


bool Scene::clipLinePers(Vec4& v1, Vec4& v2, Color& c1, Color& c2){

    bool v1_out = v1.w < 0.0;
    bool v2_out = v2.w < 0.0;
    if(v1_out < 0 && v2_out < 0){
        return false;
    }


    bool v1_or_v2_out = v1_out ^ v2_out;
    c1 = *colorsOfVertices[v1.colorId];
    c2 = *colorsOfVertices[v2.colorId];

    if(v1_or_v2_out){
        Vec4 diff = v2-v1;
        Color c_diff = c2-c1;
        double t = (0.00001 - v1.w) / (v2.w-v1.w);
        if(v1_out){
            Vec4 v1_new = v1 + (diff * t);
            v1.x = v1_new.x;
            v1.y = v1_new.y;
            v1.z = v1_new.z;
            v1.w = v1_new.w;
            c1 = *colorsOfVertices[v1.colorId] + (c_diff * t);
        }else{
            Vec4 v2_new = v1 + (diff * t);
            v2.x = v2_new.x;
            v2.y = v2_new.y;
            v2.z = v2_new.z;
            v2.w = v2_new.w;
            c2 = *colorsOfVertices[v1.colorId] + (c_diff * t);
        }
    }
        

            
    double t_e = 0.0;
    double t_l = 1.0;

    Vec4 diff = v2-v1;
    if(visible(v1.w-v1.x, diff.x-diff.w, t_e, t_l)){
        if(visible(v1.w+v1.x, -diff.w-diff.x, t_e, t_l)){
            if(visible(v1.w-v1.y, diff.y-diff.w, t_e, t_l)){
                if(visible(v1.w+v1.y, -diff.w-diff.y, t_e, t_l)){
                    if(visible(v1.w-v1.z, diff.z-diff.w, t_e, t_l)){
                        if(visible(v1.w+v1.z, -diff.w-diff.z, t_e, t_l)){
                            Vec4 v1_new = v1;
                            Vec4 v2_new = v2;

                            
                            Color c_diff = c2-c1;
                            Color c1_new = c1;
                            Color c2_new = c2;
                            if(t_l < 1){
                                v2_new     = v1 + (diff * t_l);
                                c2_new     = c1 + (c_diff * t_l);
                            }
                            if(t_e > 0){
                                v1_new     = v1 + (diff * t_e);
                                c1_new        = c1 + (c_diff * t_e);
                            }
                            v1 = v1_new;
                            v2 = v2_new;
                            c1 = c1_new;
                            c2 = c2_new;

                            return true;
                        }
                    }
                }
            }
        }
    }

    return false;
}

void Scene::initialize_tri_normals(Mesh* mesh){
    mesh->tri_normals.reserve(mesh->numberOfTriangles);
    for(int i = 0; i < mesh->numberOfTriangles; ++i){
        Vec3 v1_v0 = *(vertices[mesh->triangles[i].vertexIds[1]-1]) - *(vertices[mesh->triangles[i].vertexIds[0]-1]);
        Vec3 v2_v0 = *(vertices[mesh->triangles[i].vertexIds[2]-1]) - *(vertices[mesh->triangles[i].vertexIds[0]-1]);
        Vec3 normal = cross_prod(v1_v0, v2_v0);
        normal.normalize();
        mesh->tri_normals.push_back(normal);
    }
} 

void Scene::main_loop(){

    if(cullingEnabled){
        for(int i = 0; i < meshes.size(); ++i){
            initialize_tri_normals(meshes[i]);
            initialize_tform_mat(meshes[i]);
        }
    }else{
        for(int i = 0; i < meshes.size();++i){
            initialize_tform_mat(meshes[i]);
        }
    }



    for(int i = 0; i < cameras.size(); ++i){

        cameras[i]->initialize_camera();
        initializeImage(cameras[i]);
        render_image(*cameras[i]);
        writeImageToPPMFile(cameras[i]);
    }

}




void Scene::draw_solid(Color* c_0, Color* c_1, Color* c_2, int x_0, int y_0, int x_1, int y_1, int x_2, int y_2, float z_0, float z_1, float z_2){
    int ymin = min(y_0, min(y_1, y_2));
    int ymax = max(y_0, max(y_1, y_2));
    int xmin = min(x_0, min(x_1, x_2));
    int xmax = max(x_0, max(x_1, x_2));

    double a = y_0-y_1;
    double b = y_1-y_2;
    double c = y_2-y_0;

    double v = x_1-x_0;
    double w = x_2-x_1;
    double z = x_0-x_2;

    double lit_0 = (x_0+0.5)*(y_1+0.5) - (y_0+0.5)*(x_1+0.5);
    double lit_1 = (x_1+0.5)*(y_2+0.5) - (y_1+0.5)*(x_2+0.5);
    double lit_2 = (x_2+0.5)*(y_0+0.5) - (y_2+0.5)*(x_0+0.5);

    double alpha_denom = (x_0+0.5)*b+(y_0+0.5)*w+lit_1;
    double beta_denom = (x_1+0.5)*c+(y_1+0.5)*z+lit_2;
    double gamma_denom = (x_2+0.5)*a+(y_2+0.5)*v+lit_0;
    for(int j = ymin; j <= ymax; ++j){
        for(int i = xmin; i <= xmax; ++i){
            double alpha = ((i+0.5)*b + (j+0.5)*w + lit_1)/alpha_denom;
            double beta = ((i+0.5)*c + (j+0.5)*z + lit_2)/beta_denom;
            double gamma = ((i+0.5)*a + (j+0.5)*v + lit_0)/gamma_denom;
            if(alpha >= 0 && beta >= 0 && gamma >= 0){

                float depth = interpolate_depth(z_0, z_1, z_2, alpha, beta, gamma);

                if(depth < depth_buffer[i][j]){ 
                    image[i][j] = interpolate_colors(c_0, c_1, c_2, alpha, beta, gamma);
                    depth_buffer[i][j] = depth;
                }
                
            }
        }
    }

}




void Scene::mpLowHelper(int x0, int y0, int x1, int y1, float z0, float z1, Color* c0, Color* c1){
    int dx = x1-x0;
    int dy = y1-y0;
    int y_inc = 1;
    if(dy < 0){
        y_inc = -1;
        dy = -dy;
    }

    Color base = *c0;
    Color inc = (*c1 - *c0);
    inc.r /= dx;
    inc.g /= dx;
    inc.b /= dx;

    float z_inc = (z1 - z0) / (float) dx;
    int d = 2*dy - dx;
    int y = y0;
    for(int i = x0; i <= x1; ++i){
        if(depth_buffer[i][y] > z0){
            image[i][y] = base;
            depth_buffer[i][y] = z0;
        }
        if(d > 0){
            y += y_inc;
            d += 2*(dy-dx);
        }else{
            d += 2*dy;
        }
        base.r += inc.r;
        base.g += inc.g;
        base.b += inc.b;
        z0 += z_inc;
    }

}

void Scene::mpHighHelper(int x0, int y0, int x1, int y1, float z0, float z1, Color* c0, Color* c1){
    int dx = x1-x0;
    int dy = y1-y0;
    int x_inc = 1;

    if(dx < 0){
        x_inc = -1;
        dx = -dx;
    }

    Color base = *c0;
    Color inc = (*c1 - *c0);
    inc.r /= dy;
    inc.g /= dy;
    inc.b /= dy;

    float z_inc = (z1 - z0) / (float) dy;

    int d = 2*dx-dy;
    int x = x0;
    for(int j = y0; j <= y1; ++j){
	if(depth_buffer[x][j] > z0){
            image[x][j] = base;
            depth_buffer[x][j] = z0;
        }
        if(d > 0){
            x += x_inc;
            d += 2*(dx-dy);
        }else{
            d += 2*dx;
        }
        base.r += inc.r;
        base.g += inc.g;
        base.b += inc.b;
        z0 += z_inc;
    }

}


void Scene::drawMidpoint(int x0, int y0, int x1, int y1, float z0, float z1, Color* c_0, Color* c_1){
    if(abs(y1-y0) < abs(x1-x0)){
        if(x0 > x1){
            mpLowHelper(x1, y1, x0, y0, z1, z0, c_1, c_0);
        }else{
            mpLowHelper(x0, y0, x1, y1, z0, z1, c_0, c_1);
        }
    }else{
        if(y0 > y1){
            mpHighHelper(x1, y1, x0, y0, z1, z0, c_1, c_0);
        }else{
            mpHighHelper(x0, y0, x1, y1, z0, z1, c_0, c_1);
        }
    }
}

void Scene::drawSolidPers(Camera& camera, Mesh* mesh){
    Matrix4 inv_tform;
    if(cullingEnabled){
        inv_tform = invTpose(mesh->tform_mat);
    }
    Matrix4 complete_mat = matxmat_mult(camera.projxcam_mat, mesh->tform_mat);
    double a = ((double)camera.horRes)/2.0;
    double b = ((double)camera.verRes)/2.0;
    int e = camera.horRes-1;
    int f = camera.verRes-1;
    double c = ((double)e)/2.0;
    double d = ((double)f)/2.0;
    for(int i = 0; i < mesh->numberOfTriangles; ++i){

        Vec4 v0{*vertices[mesh->triangles[i].vertexIds[0]-1], mesh->triangles[i].vertexIds[0]-1}; 

        v0 = matxvec_mult(mesh->tform_mat, v0);

        if(cullingEnabled){
            Vec3 normal_new = findNewNormal(inv_tform, mesh->tri_normals[i]);
            Vec3 ctov{v0.x-camera.pos.x,v0.y-camera.pos.y,v0.z-camera.pos.z,-1};
            if(!(dot_prod(ctov, normal_new) < 0)){ 
                continue;
            }
        }

        Vec4 v1{*vertices[mesh->triangles[i].vertexIds[1]-1], mesh->triangles[i].vertexIds[1]-1};
        Vec4 v2{*vertices[mesh->triangles[i].vertexIds[2]-1], mesh->triangles[i].vertexIds[2]-1};

        v0 = matxvec_mult(camera.projxcam_mat, v0);
        v1 = matxvec_mult(complete_mat, v1);
        v2 = matxvec_mult(complete_mat, v2);

        double x0 = v0.x/v0.w;
        double y0 = v0.y/v0.w;
        double x1 = v1.x/v1.w;
        double y1 = v1.y/v1.w;
        double x2 = v2.x/v2.w;
        double y2 = v2.y/v2.w;

        double z0 = v0.z/v0.w;
        double z1 = v1.z/v1.w;
        double z2 = v2.z/v2.w;
        z0 = 0.5*z0 + 0.5;
        z1 = 0.5*z1 + 0.5;
        z2 = 0.5*z2 + 0.5;

        x0 = a*x0 + c;
        y0 = b*y0 + d;
        x1 = a*x1 + c;
        y1 = b*y1 + d;
        x2 = a*x2 + c;
        y2 = b*y2 + d;
        int x0_pix = max(0, min((int)(x0+0.5), e));
        int y0_pix = max(0, min((int)(y0+0.5), f));
        int x1_pix = max(0, min((int)(x1+0.5), e));
        int y1_pix = max(0, min((int)(y1+0.5), f));
        int x2_pix = max(0, min((int)(x2+0.5), e));
        int y2_pix = max(0, min((int)(y2+0.5), f));
        draw_solid(colorsOfVertices[v0.colorId], colorsOfVertices[v1.colorId], colorsOfVertices[v2.colorId], x0_pix, y0_pix, x1_pix, y1_pix, x2_pix, y2_pix, z0, z1, z2);
    }
}



void Scene::drawWFPers(Camera& camera, Mesh* mesh){
    Matrix4 inv_tform;
    if(cullingEnabled){
        inv_tform = invTpose(mesh->tform_mat);
    }
    double a = ((double)camera.horRes)/2.0;
    double b = ((double)camera.verRes)/2.0;
    int e = camera.horRes-1;
    int f = camera.verRes-1;
    double c = ((double)e)/2.0;
    double d = ((double)f)/2.0;


    Matrix4 complete_mat = matxmat_mult(camera.projxcam_mat, mesh->tform_mat);
    for(int i = 0; i < mesh->numberOfTriangles; ++i){

        Vec4 v0{*vertices[mesh->triangles[i].vertexIds[0]-1], mesh->triangles[i].vertexIds[0]-1}; 

        v0 = matxvec_mult(mesh->tform_mat, v0);

        if(cullingEnabled){
            Vec3 normal_new = findNewNormal(inv_tform, mesh->tri_normals[i]);
            Vec3 ctov{v0.x-camera.pos.x,v0.y-camera.pos.y,v0.z-camera.pos.z,-1};
            if(!(dot_prod(ctov, normal_new) < 0)){ 
                continue;
            }
        }

        Vec4 v1{*vertices[mesh->triangles[i].vertexIds[1]-1], mesh->triangles[i].vertexIds[1]-1};
        Vec4 v2{*vertices[mesh->triangles[i].vertexIds[2]-1], mesh->triangles[i].vertexIds[2]-1};

        Color c1;
        Color c2;

        v0 = matxvec_mult(camera.projxcam_mat, v0);
        v1 = matxvec_mult(complete_mat, v1);
        v2 = matxvec_mult(complete_mat, v2);

        Vec4 tbc_v0 = v0;
        Vec4 tbc_v1 = v1;
        if(clipLinePers(tbc_v0, tbc_v1, c1, c2)){
        


            double x0 = tbc_v0.x/tbc_v0.w;
            double y0 = tbc_v0.y/tbc_v0.w;
            double x1 = tbc_v1.x/tbc_v1.w;
            double y1 = tbc_v1.y/tbc_v1.w;
            
            float z0 = tbc_v0.z/tbc_v0.w;
            float z1 = tbc_v1.z/tbc_v1.w;
            z0 = 0.5f*z0 + 0.5f;
            z1 = 0.5f*z1 + 0.5f;


            x0 = a*x0 + c;
            y0 = b*y0 + d;
            x1 = a*x1 + c;
            y1 = b*y1 + d;

            int x0_pix = max(0, min((int)(x0+0.5), e));
            int y0_pix = max(0, min((int)(y0+0.5), f));
            int x1_pix = max(0, min((int)(x1+0.5), e));
            int y1_pix = max(0, min((int)(y1+0.5), f));
            drawMidpoint(x0_pix, y0_pix, x1_pix, y1_pix, z0, z1, &c1, &c2); 
        }
        Vec4 v1_clone = v1;
        Vec4 v2_clone = v2;
        if(clipLinePers(v1_clone, v2_clone, c1, c2)){

            double x1 = v1_clone.x/v1_clone.w;
            double y1 = v1_clone.y/v1_clone.w;
            double x2 = v2_clone.x/v2_clone.w;
            double y2 = v2_clone.y/v2_clone.w;
            
            float z1 = v1_clone.z/v1_clone.w;
            float z2 = v2_clone.z/v2_clone.w;
            z1 = 0.5f*z1 + 0.5f;
            z2 = 0.5f*z2 + 0.5f;

            x1 = a*x1 + c;
            y1 = b*y1 + d;
            x2 = a*x2 + c;
            y2 = b*y2 + d;


            int x1_pix = max(0, min((int)(x1+0.5), e));
            int y1_pix = max(0, min((int)(y1+0.5), f));
            int x2_pix = max(0, min((int)(x2+0.5), e));
            int y2_pix = max(0, min((int)(y2+0.5), f));
            drawMidpoint(x1_pix, y1_pix, x2_pix, y2_pix, z1, z2, &c1, &c2);
        }
        Vec4 v0_cl = v0;
        Vec4 v2_cl = v2;

        if(clipLinePers(v2_cl, v0_cl, c2, c1)){

            double x0 = v0_cl.x/v0_cl.w;
            double y0 = v0_cl.y/v0_cl.w;
            double x2 = v2_cl.x/v2_cl.w;
            double y2 = v2_cl.y/v2_cl.w;
            
            float z0 = v0_cl.z/v0_cl.w;
            float z2 = v2_cl.z/v2_cl.w;
            z0 = 0.5f*z0 + 0.5f;
            z2 = 0.5f*z2 + 0.5f;

            x0 = a*x0 + c;
            y0 = b*y0 + d;
            x2 = a*x2 + c;
            y2 = b*y2 + d;


            int x0_pix = max(0, min((int)(x0+0.5), e));
            int y0_pix = max(0, min((int)(y0+0.5), f));
            int x2_pix = max(0, min((int)(x2+0.5), e));
            int y2_pix = max(0, min((int)(y2+0.5), f));
            drawMidpoint(x2_pix, y2_pix, x0_pix, y0_pix, z2, z0, &c2, &c1);
        }
    }

}

void Scene::drawSolidOrtho(Camera& camera, Mesh* mesh){
    Matrix4 inv_tform;
    if(cullingEnabled){
        inv_tform = invTpose(mesh->tform_mat);
    }
    Matrix4 complete_mat = matxmat_mult(camera.projxcam_mat, mesh->tform_mat);
    double a = ((double)camera.horRes)/2.0;
    double b = ((double)camera.verRes)/2.0;
    int e = camera.horRes-1;
    int f = camera.verRes-1;
    double c = ((double)e)/2.0;
    double d = ((double)f)/2.0;
    for(int i = 0; i < mesh->numberOfTriangles; ++i){
        Vec4 v0{*vertices[mesh->triangles[i].vertexIds[0]-1], mesh->triangles[i].vertexIds[0]-1}; 


        v0 = matxvec_mult(mesh->tform_mat, v0);

        if(cullingEnabled){
            Vec3 normal_new = findNewNormal(inv_tform, mesh->tri_normals[i]);
            normal_new.normalize();
            Vec3 ctov{v0.x-camera.pos.x,v0.y-camera.pos.y,v0.z-camera.pos.z,-1};
            ctov.normalize();
            if(!(dot_prod(ctov, normal_new) < 0)){ 
                continue;
            }
        }

        Vec4 v1{*vertices[mesh->triangles[i].vertexIds[1]-1], mesh->triangles[i].vertexIds[1]-1};
        Vec4 v2{*vertices[mesh->triangles[i].vertexIds[2]-1], mesh->triangles[i].vertexIds[2]-1};

        v0 = matxvec_mult(camera.projxcam_mat, v0);
        v1 = matxvec_mult(complete_mat, v1);
        v2 = matxvec_mult(complete_mat, v2);

        double x0 = v0.x;
        double y0 = v0.y;
        double x1 = v1.x;
        double y1 = v1.y;
        double x2 = v2.x;
        double y2 = v2.y;

        double z0 = v0.z;
        double z1 = v1.z;
        double z2 = v2.z;
        z0 = 0.5*z0 + 0.5;
        z1 = 0.5*z1 + 0.5;
        z2 = 0.5*z2 + 0.5;

        x0 = a*x0 + c;
        y0 = b*y0 + d;
        x1 = a*x1 + c;
        y1 = b*y1 + d;
        x2 = a*x2 + c;
        y2 = b*y2 + d;
        int x0_pix = max(0, min((int)(x0+0.5), e));
        int y0_pix = max(0, min((int)(y0+0.5), f));
        int x1_pix = max(0, min((int)(x1+0.5), e));
        int y1_pix = max(0, min((int)(y1+0.5), f));
        int x2_pix = max(0, min((int)(x2+0.5), e));
        int y2_pix = max(0, min((int)(y2+0.5), f));
        draw_solid(colorsOfVertices[v0.colorId], colorsOfVertices[v1.colorId], colorsOfVertices[v2.colorId], x0_pix, y0_pix, x1_pix, y1_pix, x2_pix, y2_pix, z0, z1, z2);
    }

}

void Scene::drawWFOrtho(Camera& camera, Mesh* mesh){
    Matrix4 inv_tform;
    if(cullingEnabled){
        inv_tform = invTpose(mesh->tform_mat);
    }
    Matrix4 complete_mat = matxmat_mult(camera.projxcam_mat, mesh->tform_mat);
    double a = ((double)camera.horRes)/2.0;
    double b = ((double)camera.verRes)/2.0;
    int e = camera.horRes-1;
    int f = camera.verRes-1;
    double c = ((double)e)/2.0;
    double d = ((double)f)/2.0;
    for(int i = 0; i < mesh->numberOfTriangles; ++i){

        Vec4 v0{*vertices[mesh->triangles[i].vertexIds[0]-1], mesh->triangles[i].vertexIds[0]-1}; 

        v0 = matxvec_mult(mesh->tform_mat, v0);

        if(cullingEnabled){
            Vec3 normal_new = findNewNormal(inv_tform, mesh->tri_normals[i]);
            normal_new.normalize();
            Vec3 ctov{v0.x-camera.pos.x,v0.y-camera.pos.y,v0.z-camera.pos.z,-1};
            ctov.normalize();
            if(!(dot_prod(ctov, normal_new) < 0)){ 
                continue;
            }
        }

        Vec4 v1{*vertices[mesh->triangles[i].vertexIds[1]-1], mesh->triangles[i].vertexIds[1]-1};
        Vec4 v2{*vertices[mesh->triangles[i].vertexIds[2]-1], mesh->triangles[i].vertexIds[2]-1};

        v0 = matxvec_mult(camera.projxcam_mat, v0);
        v1 = matxvec_mult(complete_mat, v1);
        v2 = matxvec_mult(complete_mat, v2);

        Vec4 tbc_v0 = v0;
        Vec4 tbc_v1 = v1;
        Color c1;
        Color c2;
        if(clipLineOrtho(tbc_v0, tbc_v1, c1, c2)){
            double x0 = tbc_v0.x;
            double y0 = tbc_v0.y;
            double x1 = tbc_v1.x;
            double y1 = tbc_v1.y;

            x0 = a*x0 + c;
            y0 = b*y0 + d;
            x1 = a*x1 + c;
            y1 = b*y1 + d;

            float z0 = 0.5f*tbc_v0.z + 0.5f;
            float z1 = 0.5f*tbc_v1.z + 0.5f;

            int x0_pix = max(0, min((int)(x0+0.5), e));
            int y0_pix = max(0, min((int)(y0+0.5), f));
            int x1_pix = max(0, min((int)(x1+0.5), e));
            int y1_pix = max(0, min((int)(y1+0.5), f));
            drawMidpoint(x0_pix, y0_pix, x1_pix, y1_pix, z0, z1, &c1, &c2);
        }
        Vec4 v1_clone = v1;
        Vec4 v2_clone = v2;
        if(clipLineOrtho(v1_clone, v2_clone, c1, c2)){
            double x1 = v1_clone.x;
            double y1 = v1_clone.y;
            double x2 = v2_clone.x;
            double y2 = v2_clone.y;

            x1 = a*x1 + c;
            y1 = b*y1 + d;
            x2 = a*x2 + c;
            y2 = b*y2 + d;

            float z1 = 0.5f*v1_clone.z + 0.5f;
            float z2 = 0.5f*v2_clone.z + 0.5f;

            int x1_pix = max(0, min((int)(x1+0.5), e));
            int y1_pix = max(0, min((int)(y1+0.5), f));
            int x2_pix = max(0, min((int)(x2+0.5), e));
            int y2_pix = max(0, min((int)(y2+0.5), f));
            drawMidpoint(x1_pix, y1_pix, x2_pix, y2_pix, z1, z2, &c1, &c2);
        }
        Vec4 v0_cl = v0;
        Vec4 v2_cl = v2;
        if(clipLineOrtho(v2_cl, v0_cl, c2, c1)){
            double x0 = v0_cl.x;
            double y0 = v0_cl.y;
            double x2 = v2_cl.x;
            double y2 = v2_cl.y;

            x0 = a*x0 + c;
            y0 = b*y0 + d;
            x2 = a*x2 + c;
            y2 = b*y2 + d;

            float z0 = 0.5f*v0_cl.z + 0.5f;
            float z2 = 0.5f*v2_cl.z + 0.5f;

            int x0_pix = max(0, min((int)(x0+0.5), e));
            int y0_pix = max(0, min((int)(y0+0.5), f));
            int x2_pix = max(0, min((int)(x2+0.5), e));
            int y2_pix = max(0, min((int)(y2+0.5), f));
            drawMidpoint(x2_pix, y2_pix, x0_pix, y0_pix, z2, z0, &c2, &c1);
        }
    }


}


void Scene::render_image(Camera& camera)
{
    if(camera.projectionType){
        for(int i = 0; i < meshes.size();++i){
            if(meshes[i]->type){
                drawSolidPers(camera, meshes[i]);
            }else{
                drawWFPers(camera, meshes[i]);
            }
        }
    }else{
        for(int i = 0; i < meshes.size();++i){
            if(meshes[i]->type){
                drawSolidOrtho(camera, meshes[i]);
            }else{
                drawWFOrtho(camera, meshes[i]);
            }
        }
    }

}

/*
    Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
    const char *str;
    XMLDocument xmlDoc;
    XMLElement *pElement;

    xmlDoc.LoadFile(xmlPath);

    XMLNode *pRoot = xmlDoc.FirstChild();

    // read background color
    pElement = pRoot->FirstChildElement("BackgroundColor");
    str = pElement->GetText();
    sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

    // read culling
    pElement = pRoot->FirstChildElement("Culling");
    if (pElement != NULL) {
        str = pElement->GetText();
        
        if (strcmp(str, "enabled") == 0) {
            cullingEnabled = true;
        }
        else {
            cullingEnabled = false;
        }
    }

    // read cameras
    pElement = pRoot->FirstChildElement("Cameras");
    XMLElement *pCamera = pElement->FirstChildElement("Camera");
    XMLElement *camElement;
    while (pCamera != NULL)
    {
        Camera *cam = new Camera();

        pCamera->QueryIntAttribute("id", &cam->cameraId);

        // read projection type
        str = pCamera->Attribute("type");

        if (strcmp(str, "orthographic") == 0) {
            cam->projectionType = 0;
        }
        else {
            cam->projectionType = 1;
        }

        camElement = pCamera->FirstChildElement("Position");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

        camElement = pCamera->FirstChildElement("Gaze");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

        camElement = pCamera->FirstChildElement("Up");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

        cam->gaze.normalize();
        cam->u = cross_prod(cam->gaze, cam->v);
        cam->u.normalize();

        cam->w = -cam->gaze;
        cam->v = cross_prod(cam->u, cam->gaze);
        cam->v.normalize();

        camElement = pCamera->FirstChildElement("ImagePlane");
        str = camElement->GetText();
        sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
               &cam->left, &cam->right, &cam->bottom, &cam->top,
               &cam->near, &cam->far, &cam->horRes, &cam->verRes);

        camElement = pCamera->FirstChildElement("OutputName");
        str = camElement->GetText();
        cam->outputFileName = string(str);

        cameras.push_back(cam);

        pCamera = pCamera->NextSiblingElement("Camera");
    }

    // read vertices
    pElement = pRoot->FirstChildElement("Vertices");
    XMLElement *pVertex = pElement->FirstChildElement("Vertex");
    int vertexId = 1;

    while (pVertex != NULL)
    {
        Vec3 *vertex = new Vec3();
        Color *color = new Color();

        vertex->colorId = vertexId;

        str = pVertex->Attribute("position");
        sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

        str = pVertex->Attribute("color");
        sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

        vertices.push_back(vertex);
        colorsOfVertices.push_back(color);

        pVertex = pVertex->NextSiblingElement("Vertex");

        vertexId++;
    }

    // read translations
    pElement = pRoot->FirstChildElement("Translations");
    XMLElement *pTranslation = pElement->FirstChildElement("Translation");
    while (pTranslation != NULL)
    {
        Translation *translation = new Translation();

        pTranslation->QueryIntAttribute("id", &translation->translationId);

        str = pTranslation->Attribute("value");
        sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

        translations.push_back(translation);

        pTranslation = pTranslation->NextSiblingElement("Translation");
    }

    // read scalings
    pElement = pRoot->FirstChildElement("Scalings");
    XMLElement *pScaling = pElement->FirstChildElement("Scaling");
    while (pScaling != NULL)
    {
        Scaling *scaling = new Scaling();

        pScaling->QueryIntAttribute("id", &scaling->scalingId);
        str = pScaling->Attribute("value");
        sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

        scalings.push_back(scaling);

        pScaling = pScaling->NextSiblingElement("Scaling");
    }

    // read rotations
    pElement = pRoot->FirstChildElement("Rotations");
    XMLElement *pRotation = pElement->FirstChildElement("Rotation");
    while (pRotation != NULL)
    {
        Rotation *rotation = new Rotation();

        pRotation->QueryIntAttribute("id", &rotation->rotationId);
        str = pRotation->Attribute("value");
        sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

        rotations.push_back(rotation);

        pRotation = pRotation->NextSiblingElement("Rotation");
    }

    // read meshes
    pElement = pRoot->FirstChildElement("Meshes");

    XMLElement *pMesh = pElement->FirstChildElement("Mesh");
    XMLElement *meshElement;
    while (pMesh != NULL)
    {
        Mesh *mesh = new Mesh();

        pMesh->QueryIntAttribute("id", &mesh->meshId);

        // read projection type
        str = pMesh->Attribute("type");

        if (strcmp(str, "wireframe") == 0) {
            mesh->type = 0;
        }
        else {
            mesh->type = 1;
        }

        // read mesh transformations
        XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
        XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

        while (pTransformation != NULL)
        {
            char transformationType;
            int transformationId;

            str = pTransformation->GetText();
            sscanf(str, "%c %d", &transformationType, &transformationId);

            mesh->transformationTypes.push_back(transformationType);
            mesh->transformationIds.push_back(transformationId);

            pTransformation = pTransformation->NextSiblingElement("Transformation");
        }

        mesh->numberOfTransformations = mesh->transformationIds.size();

        // read mesh faces
        char *row;
        char *clone_str;
        int v1, v2, v3;
        XMLElement *pFaces = pMesh->FirstChildElement("Faces");
            str = pFaces->GetText();
        clone_str = strdup(str);

        row = strtok(clone_str, "\n");
        while (row != NULL)
        {
            int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);
            
            if (result != EOF) {
                mesh->triangles.push_back(Triangle(v1, v2, v3));
            }
            row = strtok(NULL, "\n");
        }
        mesh->numberOfTriangles = mesh->triangles.size();
        meshes.push_back(mesh);

        pMesh = pMesh->NextSiblingElement("Mesh");
    }
}


void Scene::initializeImage(Camera *camera)
{
    if (this->image.empty())
    {
        vector<Color> rowOfColors(camera->verRes, backgroundColor);
        vector<float> depths(camera->verRes, 1.0f);

        for (int i = 0; i < camera->horRes; i++)
        {

            this->image.push_back(rowOfColors);
            this->depth_buffer.push_back(depths);
        }
    }
    else
    {
        for (int i = 0; i < camera->horRes; i++)
        {
            for (int j = 0; j < camera->verRes; j++)
            {
                this->image[i][j].r = this->backgroundColor.r;
                this->image[i][j].g = this->backgroundColor.g;
                this->image[i][j].b = this->backgroundColor.b;
                this->depth_buffer[i][j] = 1.0f;
            }
        }
    }
}


int Scene::my_clamp(double value)
{
    if (value >= 255.0)
        return 255;
    if (value <= 0.0)
        return 0;
    return (int)(value);
}


void Scene::writeImageToPPMFile(Camera *camera)
{
    ofstream fout;

    fout.open(camera->outputFileName.c_str());

    fout << "P3" << endl;
    fout << "# " << camera->outputFileName << endl;
    fout << camera->horRes << " " << camera->verRes << endl;
    fout << "255" << endl;

    for (int j = camera->verRes - 1; j >= 0; j--)
    {
        for (int i = 0; i < camera->horRes; i++)
        {
                fout << my_clamp(this->image[i][j].r) << " "
                 << my_clamp(this->image[i][j].g) << " "
                 << my_clamp(this->image[i][j].b) << " ";
        }
        fout << endl;
    }
    fout.close();
}


void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
    string command;

    // call command on Ubuntu
    if (osType == 1)
    {
        command = "convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

    // call command on Windows
    else if (osType == 2)
    {
        command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
        system(command.c_str());
    }

    // default action - don't do conversion
    else
    {
    }
}
