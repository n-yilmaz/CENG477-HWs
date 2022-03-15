#ifndef __HELPERS_H__
#define __HELPERS_H__

#define EPSILON 0.000000001
#define PI 3.141592635897932

#include "Matrix4.h"
#include "Vec3.h"
#include "Vec4.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Color.h"

double dot_prod(Vec3& a, Vec3& b);

void print_vec(Vec3& a);

Vec3 cross_prod(Vec3& a, Vec3& b);
/*
 * Returns an identity matrix (values on the diagonal are 1, others are 0).
*/
Matrix4 id_matrix();

Matrix4 invTpose(Matrix4& tform);

Vec3 findNewNormal(Matrix4&, Vec3&);
/*
 * Multiply matrices m1 (Matrix4) and m2 (Matrix4) and return the result matrix r (Matrix4).
 */
Matrix4 matxmat_mult(Matrix4& m1, Matrix4& m2);

/*
 * Multiply matrix m (Matrix4) with vector v (vec4) and store the result in vector r (vec4).
 */
Vec4 matxvec_mult(Matrix4& m, Vec4& v);

Color interpolate_colors(Color*, Color*, Color*, double, double, double);

float interpolate_depth(float, float, float, double, double, double);

Matrix4 makeTMat(Translation* t);

Matrix4 makeSMat(Scaling* s);

Matrix4 makeRMat(Rotation* r);

#endif
