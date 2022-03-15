#ifndef ECLIPSEMAP_H
#define ECLIPSEMAP_H

#include <vector>
#include <GL/glew.h>
#include <iostream>
#include "glm/glm/ext.hpp"
#include "Shader.h"
#include <vector>
#include "glm/glm/glm.hpp"
#include <GLFW/glfw3.h>
#include <jpeglib.h>
#include <GL/glew.h>

#define PI 3.14159265358979323
using namespace std;


class EclipseMap {
private:
    float heightFactor = 80;
    glm::vec3 lightPos = glm::vec3(0, 4000, 0);
    bool pKeyPressed = false;
    // DISPLAY SETTINGS
    enum displayFormatOptions {
        windowed = 1, fullScreen = 0
    };
    const char *windowName = "Ceng477 - HW3";
    int defaultScreenWidth = 1000;
    int defaultScreenHeight = 1000;
    int screenWidth = defaultScreenWidth;
    int screenHeight = defaultScreenHeight;
    int screenPosX = 0;
    int screenPosY = 0;

    int displayFormat = displayFormatOptions::windowed;
    // CAMERA SETTINGS

    float projectionAngle = 45;
    float aspectRatio = 1;
    float near = 0.1;
    float far = 10000;

    double mouse_x = screenWidth / 2.0;
    double mouse_y = screenHeight / 2.0;

    float startSpeed = 0;

    float speed = startSpeed;

    glm::vec3 cameraStartPosition = glm::vec3(0, 4000, 4000);
    glm::vec3 cameraStartDirection = glm::vec3(0, -1, -1);
    glm::vec3 cameraStartUp = glm::vec3(0, 0, 1);
    glm::vec3 cameraStartLeft;

    glm::vec3 cameraUp = cameraStartUp;
    glm::vec3 cameraPosition = cameraStartPosition;
    glm::vec3 cameraDirection = cameraStartDirection;
    glm::vec3 cameraLeft;

public:
    unsigned int textureColor;
    unsigned int textureGrey;
    unsigned int VAO;
    unsigned int VBO, EBO;
    float imageHeight;
    float imageWidth;
    float radius = 600;
    int horizontalSplitCount = 250;
    int verticalSplitCount = 125;

    unsigned int moonTextureColor;
    unsigned int moonVAO;
    unsigned int moonVBO, moonEBO;
    float moonImageHeight;
    float moonImageWidth;
    float moonRadius = 162;
    bool free_roam_mode = false;

    vector<float> worldVertices;
    vector<unsigned int> worldIndices;

    float own_axis_rot_amount = 180.f / 250.f;
    vector<float> moonVertices;
    vector<unsigned int> moonIndices;

    GLFWwindow *openWindow(const char *windowName, int width, int height);

    void Render(const char *coloredTexturePath, const char *greyTexturePath, const char *moonTexturePath);

    void handleKeyPress(GLFWwindow *window);

    void initColoredTexture(const char *filename, GLuint shader);

    void initGreyTexture(const char *filename, GLuint shader);

    void initMoonColoredTexture(const char *filename, GLuint shader);

    void initializeMoonVertices();

    void initializeWorldVertices();



};

#endif
