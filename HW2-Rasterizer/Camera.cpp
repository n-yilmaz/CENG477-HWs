#include "Camera.h"
#include <string>
#include <iostream>
#include <iomanip>
#include "Helpers.h"

using namespace std;

Camera::Camera() {}

Camera::Camera(int cameraId,
               int projectionType,
               Vec3 pos, Vec3 gaze,
               Vec3 u, Vec3 v, Vec3 w,
               double left, double right, double bottom, double top,
               double near, double far,
               int horRes, int verRes,
               string outputFileName)
{

    this->cameraId = cameraId;
    this->projectionType = projectionType;
    this->pos = pos;
    this->gaze = gaze;
    this->u = u;
    this->v = v;
    this->w = w;
    this->left = left;
    this->right = right;
    this->bottom = bottom;
    this->top = top;
    this->near = near;
    this->far = far;
    this->horRes = horRes;
    this->verRes = verRes;
    this->outputFileName = outputFileName;
}

Camera::Camera(const Camera &other)
{
    this->cameraId = other.cameraId;
    this->projectionType = other.projectionType;
    this->pos = other.pos;
    this->gaze = other.gaze;
    this->u = other.u;
    this->v = other.v;
    this->w = other.w;
    this->left = other.left;
    this->right = other.right;
    this->bottom = other.bottom;
    this->top = other.top;
    this->near = other.near;
    this->far = other.far;
    this->horRes = other.horRes;
    this->verRes = other.verRes;
    this->outputFileName = other.outputFileName;
}

void Camera::initialize_camera(){
    
    double a = dot_prod(u, pos);
    double b = dot_prod(v, pos);
    double c = dot_prod(w, pos);
    
    Matrix4 cam_mat{};
    cam_mat.val[0][0] = u.x;
    cam_mat.val[0][1] = u.y;
    cam_mat.val[0][2] = u.z;
    cam_mat.val[0][3] = -a;
    cam_mat.val[1][0] = v.x;
    cam_mat.val[1][1] = v.y;
    cam_mat.val[1][2] = v.z;
    cam_mat.val[1][3] = -b;
    cam_mat.val[2][0] = w.x;
    cam_mat.val[2][1] = w.y;
    cam_mat.val[2][2] = w.z;
    cam_mat.val[2][3] = -c;
    cam_mat.val[3][3] = 1;


    Matrix4 proj_mat{};

    double r_l = right-left;
    double rpl = right+left;
    double tpb = top+bottom;
    double t_b = top-bottom;
    double f_n = far-near;
    double fpn = far+near;

    if(projectionType){
        double m1 = 2*near/r_l;
        double m2 = rpl/r_l;
        double m3 = 2*near/t_b;
        double m4 = tpb/t_b;
        double m5 = -fpn/f_n;
        double m6 = -2 * far * near / f_n;
        proj_mat.val[0][0] = m1;
        proj_mat.val[0][2] = m2;
        proj_mat.val[1][1] = m3;
        proj_mat.val[1][2] = m4;
        proj_mat.val[2][2] = m5;
        proj_mat.val[2][3] = m6;
        proj_mat.val[3][2] = -1;
    }else{
        double m1 = 2/r_l;
        double m2 = -rpl/r_l;
        double m3 = 2/t_b;
        double m4 = -tpb/t_b;
        double m5 = -2/f_n;
        double m6 = -fpn/f_n;
        proj_mat.val[0][0] = m1;
        proj_mat.val[0][3] = m2;
        proj_mat.val[1][1] = m3;
        proj_mat.val[1][3] = m4;
        proj_mat.val[2][2] = m5;
        proj_mat.val[2][3] = m6;
        proj_mat.val[3][3] = 1;
    }


    this->projxcam_mat = matxmat_mult(proj_mat, cam_mat);
}

    

ostream &operator<<(ostream &os, const Camera &c)
{
    const char *camType = c.projectionType ? "perspective" : "orthographic";

    os << fixed << setprecision(6) << "Camera " << c.cameraId << " (" << camType << ") => pos: " << c.pos << " gaze: " << c.gaze << endl
       << "\tu: " << c.u << " v: " << c.v << " w: " << c.w << endl
       << fixed << setprecision(3) << "\tleft: " << c.left << " right: " << c.right << " bottom: " << c.bottom << " top: " << c.top << endl
       << "\tnear: " << c.near << " far: " << c.far << " resolutions: " << c.horRes << "x" << c.verRes << " fileName: " << c.outputFileName;

    return os;
}
