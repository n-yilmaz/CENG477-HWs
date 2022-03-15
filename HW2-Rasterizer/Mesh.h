#ifndef __MESH_H__
#define __MESH_H__

#include <vector>
#include "Triangle.h"
#include <iostream>
#include "Matrix4.h"
#include "Vec3.h"

using namespace std;

class Mesh
{

public:
    int meshId;
    int type; // 0 for wireframe, 1 for solid

    int numberOfTransformations;
    vector<int> transformationIds;
    vector<char> transformationTypes;

    Matrix4 tform_mat;

    int numberOfTriangles;
    vector<Triangle> triangles;
    vector<Vec3> tri_normals;

    Mesh();
    Mesh(int meshId, int type, int numberOfTransformations,
          vector<int> transformationIds,
          vector<char> transformationTypes,
          int numberOfTriangles,
          vector<Triangle> triangles);

    friend ostream &operator<<(ostream &os, const Mesh &m);
};

#endif
