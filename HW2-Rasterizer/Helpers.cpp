#include <iostream>
#include <cmath>
#include "Helpers.h"


using namespace std;


Vec3 cross_prod(Vec3& a, Vec3& b)
{

    double res_x = a.y * b.z - b.y * a.z;
    double res_y = b.x * a.z - a.x * b.z;
    double res_z = a.x * b.y - b.x * a.y;

    return {res_x, res_y, res_z, -1};
}


double dot_prod(Vec3& a, Vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/*
 * Prints elements in a vec3. Can be used for debugging purposes.
 */
void print_vec(Vec3& v)
{
    cout << "(" << v.x << "," << v.y << "," << v.z << ")" << endl;
}


/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
*/
Matrix4 id_matrix()
{
    Matrix4 result{};

    result.val[0][0] = result.val[1][1] = result.val[2][2] = result.val[3][3] = 1.0;

    return result;
}


Matrix4 matxmat_mult(Matrix4& m1, Matrix4& m2)
{
    Matrix4 result;
    double total;

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            total = 0;
            for (int k = 0; k < 4; k++)
            {
                total += m1.val[i][k] * m2.val[k][j];
            }

            result.val[i][j] = total;
        }
    }

    return result;
}

Matrix4 makeTMat(Translation* t){

    Matrix4 translation = id_matrix();
    translation.val[0][3] = t->tx;
    translation.val[1][3] = t->ty;
    translation.val[2][3] = t->tz;
    return translation;

}

Matrix4 makeSMat(Scaling* s){
    Matrix4 scaling{};
    scaling.val[0][0] = s->sx;
    scaling.val[1][1] = s->sy;
    scaling.val[2][2] = s->sz;
    scaling.val[3][3] = 1.0;
    return scaling;
}

Matrix4 makeRMat(Rotation* r){
    Matrix4 rot{};
    static const double deg2rad = PI / 180.0;
    double rot_cos = cos(r->angle * deg2rad);
    double one_minus = 1-rot_cos;
    double rot_sin = sin(r->angle * deg2rad);

    double len = std::sqrt(r->ux*r->ux + r->uy*r->uy + r->uz*r->uz);
    double my_ux = r->ux/len;
    double my_uy = r->uy/len;
    double my_uz = r->uz/len;

    double ux_uy = my_ux * my_uy;
    double uy_uz = my_uy * my_uz;
    double ux_uz = my_ux * my_uz;

    rot.val[0][0] = rot_cos + one_minus*my_ux*my_ux;
    rot.val[0][1] = ux_uy*one_minus - my_uz*rot_sin;
    rot.val[0][2] = ux_uz*one_minus + my_uy*rot_sin;
    rot.val[1][0] = ux_uy*one_minus + my_uz*rot_sin;
    rot.val[1][1] = rot_cos + one_minus*my_uy*my_uy;
    rot.val[1][2] = uy_uz*one_minus - my_ux*rot_sin;
    rot.val[2][0] = ux_uz*one_minus - my_uy*rot_sin;
    rot.val[2][1] = uy_uz*one_minus + my_ux*rot_sin;
    rot.val[2][2] = rot_cos + one_minus*my_uz*my_uz;
    rot.val[3][3] = 1.0;
    return rot;
}

Color interpolate_colors(Color* c0, Color* c1, Color* c2, double alpha, double beta, double gamma){
    double r_val = c0->r*alpha + c1->r*beta + c2->r*gamma;
    double g_val = c0->g*alpha + c1->g*beta + c2->g*gamma;
    double b_val = c0->b*alpha + c1->b*beta + c2->b*gamma;
    return {r_val, g_val, b_val};
}

float interpolate_depth(float z_0, float z_1, float z_2, double alpha, double beta, double gamma){

    float z0_inv = 1.f / z_0;
    float z1_inv = 1.f / z_1;
    float z2_inv = 1.f / z_2;

    return 1.f / (z0_inv * alpha + z1_inv * beta + z2_inv * gamma);
}

Matrix4 invTpose(Matrix4& tform){
    double a = tform.val[0][0] * (tform.val[1][1]*tform.val[2][2] - tform.val[1][2]*tform.val[2][1]);
    double b = tform.val[0][1] * (tform.val[1][0]*tform.val[2][2] - tform.val[1][2]*tform.val[2][0]);
    double c = tform.val[0][2] * (tform.val[1][0]*tform.val[2][1] - tform.val[1][1]*tform.val[2][0]);

    double det = a-b+c;

    det = 1.0/det;
    Matrix4 inv{};
    inv.val[0][0] = det*(tform.val[1][1]*tform.val[2][2] - tform.val[1][2]*tform.val[2][1]);
    inv.val[0][1] = det*(tform.val[1][2]*tform.val[2][0] - tform.val[1][0]*tform.val[2][2]);
    inv.val[0][2] = det*(tform.val[1][0]*tform.val[2][1] - tform.val[2][0]*tform.val[1][1]);
    inv.val[1][0] = det*(tform.val[0][2]*tform.val[2][1] - tform.val[0][1]*tform.val[2][2]);
    inv.val[1][1] = det*(tform.val[0][0]*tform.val[2][2] - tform.val[0][2]*tform.val[2][0]);
    inv.val[1][2] = det*(tform.val[2][0]*tform.val[0][1] - tform.val[0][0]*tform.val[2][1]);
    inv.val[2][0] = det*(tform.val[0][1]*tform.val[1][2] - tform.val[0][2]*tform.val[1][1]);
    inv.val[2][1] = det*(tform.val[1][0]*tform.val[0][2] - tform.val[0][0]*tform.val[1][2]);
    inv.val[2][2] = det*(tform.val[1][1]*tform.val[0][0] - tform.val[1][0]*tform.val[0][1]);

    return inv;
    
}

Vec3 findNewNormal(Matrix4& m, Vec3& v){
    double total;
    double vals[3];
    for(int i = 0; i < 3; ++i){
        total = m.val[i][0]*v.x + m.val[i][1]*v.y + m.val[i][2]*v.z;
        vals[i] = total;
    }
    return {vals[0], vals[1], vals[2], v.colorId};
} 
        

Vec4 matxvec_mult(Matrix4& m, Vec4& v)
{
    double values[4];
    double total;

    for (int i = 0; i < 4; i++)
    {
        total = m.val[i][0]*v.x+m.val[i][1]*v.y+m.val[i][2]*v.z+m.val[i][3]*v.w;
        values[i] = total;
    }

    return Vec4(values[0], values[1], values[2], values[3], v.colorId);
}
