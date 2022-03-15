#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"

using namespace std;

class Scene
{
private:

    Color backgroundColor;
    bool cullingEnabled;

    vector< vector<Color> > image;
    vector< vector<float> > depth_buffer;
    vector< Camera* > cameras;

    vector< Vec3* > vertices;
    vector< Color* > colorsOfVertices;

    vector< Scaling* > scalings;
    vector< Rotation* > rotations;
    vector< Translation* > translations;

    vector< Mesh* > meshes;

public:

    Scene(const char *xmlPath);

    void initializeImage(Camera* camera);

    void main_loop();

    void render_image(Camera& camera);

    void initialize_tri_normals(Mesh* mesh);

    int my_clamp(double value);


    void drawMidpoint(int, int, int, int, float, float, Color*, Color*);
    void mpLowHelper(int, int, int, int, float, float, Color*, Color*);
    void mpHighHelper(int, int, int, int, float, float, Color*, Color*);

    bool visible(double, double, double&, double&);


    bool clipLinePers(Vec4&, Vec4&, Color&, Color&);

    bool clipLineOrtho(Vec4&, Vec4&, Color&, Color&);

    void drawSolidPers(Camera& camera, Mesh* mesh);
    void drawWFPers(Camera&, Mesh*);
    void drawSolidOrtho(Camera&, Mesh*);
    void drawWFOrtho(Camera&, Mesh*);


    void draw_solid(Color*, Color*, Color*, int, int, int, int, int, int, float, float, float);
    void initialize_tform_mat(Mesh* mesh);
    void writeImageToPPMFile(Camera* camera);
    void convertPPMToPNG(string ppmFileName, int osType);
};

#endif
