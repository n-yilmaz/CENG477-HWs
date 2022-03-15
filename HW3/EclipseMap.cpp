#include "EclipseMap.h"
#include <cmath>

using namespace std;


void EclipseMap::initializeMoonVertices(){

    float x, y, z, xy;
    float nx, ny, nz;
    float lengthInv = 1.0f / moonRadius;
    float s, t;

    float sectorStep = 2.f * PI / horizontalSplitCount;
    float stackStep = PI / verticalSplitCount;

    float sectorAngle, stackAngle;
    float piovertwo = PI/2.f;
    for(int i = 0; i <= verticalSplitCount; ++i){

        stackAngle = piovertwo - i*stackStep;
        xy = moonRadius * cosf(stackAngle);
        z = moonRadius * sinf(stackAngle);

        for(int j = 0; j <= horizontalSplitCount; ++j){

            sectorAngle = j*sectorStep;

            x = xy * cosf(sectorAngle);
            y = xy * sinf(sectorAngle);
            moonVertices.push_back(x);
            moonVertices.push_back(y);
            moonVertices.push_back(z);


            nx = x * lengthInv;
            ny = y * lengthInv;
            nz = z * lengthInv;

            moonVertices.push_back(nx);
            moonVertices.push_back(ny);
            moonVertices.push_back(nz);

            s = (float)j / horizontalSplitCount;
            t = (float)i / verticalSplitCount;
            moonVertices.push_back(s);
            moonVertices.push_back(t);



        }

    }


    int k1, k2;
    for(int i = 0; i < verticalSplitCount; ++i){

        k1 = i * (horizontalSplitCount+1);
        k2 = k1 + horizontalSplitCount + 1;

        for(int j = 0; j < horizontalSplitCount; ++j, ++k1, ++k2){

            if(i != 0){
                moonIndices.push_back(k1);
                moonIndices.push_back(k2);
                moonIndices.push_back(k1+1);
            }

            if(i != (verticalSplitCount-1)){

                moonIndices.push_back(k1+1);
                moonIndices.push_back(k2);
                moonIndices.push_back(k2+1);

            }

        }

    }

}

void EclipseMap::initializeWorldVertices(){

    float x, y, z, xy;
    float nx, ny, nz;
    float lengthInv = 1.0f / radius;
    float s, t;

    float sectorStep = 2.f * PI / horizontalSplitCount;
    float stackStep = PI / verticalSplitCount;

    float sectorAngle, stackAngle;
    float piovertwo = PI/2.f;

    for(int i = 0; i <= verticalSplitCount; ++i){

        stackAngle = piovertwo - i*stackStep;
        xy = radius * cosf(stackAngle);
        z = radius * sinf(stackAngle);

        for(int j = 0; j <= horizontalSplitCount; ++j){

            sectorAngle = j*sectorStep;

            x = xy * cosf(sectorAngle);
            y = xy * sinf(sectorAngle);
            worldVertices.push_back(x);
            worldVertices.push_back(y);
            worldVertices.push_back(z);


            nx = x * lengthInv;
            ny = y * lengthInv;
            nz = z * lengthInv;

            worldVertices.push_back(nx);
            worldVertices.push_back(ny);
            worldVertices.push_back(nz);

            s = (float)j;
            t = (float)i;
            worldVertices.push_back(s);
            worldVertices.push_back(t);



        }

    }



    int k1, k2;
    for(int i = 0; i < verticalSplitCount; ++i){

        k1 = i * (horizontalSplitCount+1);
        k2 = k1 + horizontalSplitCount + 1;

        for(int j = 0; j < horizontalSplitCount; ++j, ++k1, ++k2){

            if(i != 0){
                worldIndices.push_back(k1);
                worldIndices.push_back(k2);
                worldIndices.push_back(k1+1);
            }

            if(i != (verticalSplitCount-1)){

                worldIndices.push_back(k1+1);
                worldIndices.push_back(k2);
                worldIndices.push_back(k2+1);

            }

        }

    }

}

void fb_callback_function(GLFWwindow* window, int width, int height){

    glViewport(0, 0, width, height);

}


void EclipseMap::Render(const char *coloredTexturePath, const char *greyTexturePath, const char *moonTexturePath) {
    // Open window
    GLFWwindow *window = openWindow(windowName, screenWidth, screenHeight);
    glfwSetFramebufferSizeCallback(window, fb_callback_function);

    // Moon commands
    GLuint moonShaderID = initShaders("moonShader.vert", "moonShader.frag");

    initMoonColoredTexture(moonTexturePath, moonShaderID);

    initializeMoonVertices();

    glGenVertexArrays(1, &moonVAO);
    glGenBuffers(1, &moonVBO);
    glGenBuffers(1, &moonEBO);

    glBindVertexArray(moonVAO);
    glBindBuffer(GL_ARRAY_BUFFER, moonVBO);
    glBufferData(GL_ARRAY_BUFFER, (unsigned int) moonVertices.size() * sizeof(float), moonVertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, moonEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, (unsigned int) moonIndices.size() * sizeof(unsigned int), moonIndices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8*sizeof(float), (void*) 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8*sizeof(float), (void*) (3*sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8*sizeof(float), (void*) (6*sizeof(float)));
    glEnableVertexAttribArray(2);

    glBindVertexArray(0);
    

    // World commands
    GLuint worldShaderID = initShaders("worldShader.vert", "worldShader.frag");

    initColoredTexture(coloredTexturePath, worldShaderID);

    initializeWorldVertices();

    initGreyTexture(greyTexturePath, worldShaderID);


    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, (unsigned int) worldVertices.size() * sizeof(float), worldVertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, (unsigned int) worldIndices.size() * sizeof(unsigned int), worldIndices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8*sizeof(float), (void*) 0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8*sizeof(float), (void*) (3*sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8*sizeof(float), (void*) (6*sizeof(float)));
    glEnableVertexAttribArray(2);

    glBindVertexArray(0);
    
    // Program locations
    GLint moon_model, moon_cam_pos, moon_MVP_loc;
    GLint world_model, world_height, world_cam_pos, world_MVP_loc;

    cameraStartDirection = glm::normalize(cameraStartDirection);
    cameraStartLeft = glm::normalize(glm::cross(cameraStartDirection, cameraStartUp));
    cameraStartUp = glm::normalize(glm::cross(cameraStartLeft, cameraStartDirection));
    cameraDirection = cameraStartDirection;
    cameraLeft = cameraStartLeft;
    cameraUp = cameraStartUp;
    //
    glEnable(GL_DEPTH_TEST);

    glm::mat4 view_matrix = glm::lookAt(cameraPosition, cameraPosition + cameraDirection, cameraUp);

    glm::mat4 moon_model_matrix;
    
    glm::mat4 perspective_matrix = glm::perspective(projectionAngle, aspectRatio, near, far);

    glUseProgram(worldShaderID);
    glUniform3fv(glGetUniformLocation(worldShaderID, "lightPosition"), 1, glm::value_ptr(lightPos));
    world_model = glGetUniformLocation(worldShaderID, "ModelMatrix");
    world_height = glGetUniformLocation(worldShaderID, "heightFactor");
    world_cam_pos = glGetUniformLocation(worldShaderID, "cameraPosition");
    world_MVP_loc = glGetUniformLocation(worldShaderID, "MVP");

    glUseProgram(moonShaderID);
    glUniform3fv(glGetUniformLocation(moonShaderID, "lightPosition"), 1, glm::value_ptr(lightPos));
    moon_model = glGetUniformLocation(moonShaderID, "ModelMatrix");
    moon_cam_pos = glGetUniformLocation(moonShaderID, "cameraPosition");
    moon_MVP_loc = glGetUniformLocation(moonShaderID, "MVP");

    float self_rot_angle = 0.f;
    float moon_around_world_angle = 0.f;
    glm::mat4 self_rot;
    glm::mat4 moon_world_rot;
    glm::mat4 moon_translate = glm::translate(glm::vec3(0.f, 2600.f, 0.f));
    glm::mat4 moon_MVP;
    glm::mat4 world_MVP;

    do {
        

        glClearStencil(0);
        glClearDepth(1.0f);
        glClearColor(0, 0, 0, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        handleKeyPress(window);

        cameraPosition = free_roam_mode ? cameraPosition : cameraPosition + cameraDirection * speed;


        view_matrix = glm::lookAt(cameraPosition, cameraPosition + cameraDirection, cameraUp);

        self_rot_angle += own_axis_rot_amount;
        self_rot_angle = self_rot_angle >= 360.f ? self_rot_angle-360.f : self_rot_angle;

        moon_around_world_angle += 0.2f;
        moon_around_world_angle = moon_around_world_angle >= 360.f ? moon_around_world_angle-360.f : moon_around_world_angle;


        moon_model_matrix = glm::mat4(1.0f);
        self_rot = glm::rotate(glm::radians(self_rot_angle), glm::vec3(0.f, 0.f, 1.f));
        moon_world_rot = glm::rotate(glm::radians(moon_around_world_angle), glm::vec3(0.f, 0.f, -1.f));
        moon_model_matrix = moon_world_rot * moon_translate * self_rot * moon_model_matrix;

        moon_MVP = perspective_matrix * view_matrix * moon_model_matrix;

        glUseProgram(moonShaderID);
        glUniformMatrix4fv(moon_MVP_loc, 1, GL_FALSE, glm::value_ptr(moon_MVP));
        glUniform3fv(moon_cam_pos, 1, glm::value_ptr(cameraPosition));
        glUniformMatrix4fv(moon_model, 1, GL_FALSE, glm::value_ptr(moon_model_matrix));
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, moonTextureColor);
        glBindVertexArray(moonVAO);
        glDrawElements(GL_TRIANGLES, (unsigned int)moonIndices.size(), GL_UNSIGNED_INT, 0);
        
        
        /*************************/

        world_MVP = perspective_matrix * view_matrix * self_rot;
        glUseProgram(worldShaderID);
        glUniformMatrix4fv(world_MVP_loc, 1, GL_FALSE, glm::value_ptr(world_MVP));
        glUniform3fv(world_cam_pos, 1, glm::value_ptr(cameraPosition));
        glUniformMatrix4fv(world_model, 1, GL_FALSE, glm::value_ptr(self_rot));
        glUniform1f(world_height, heightFactor);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, textureColor);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, textureGrey);
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, (unsigned int)worldIndices.size(), GL_UNSIGNED_INT, 0);
        

        glfwSwapBuffers(window);
        glfwPollEvents();
    } while (!glfwWindowShouldClose(window));

    // Delete buffers
    glDeleteBuffers(1, &moonVAO);
    glDeleteBuffers(1, &moonVBO);
    glDeleteBuffers(1, &moonEBO);

    glDeleteBuffers(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
   
    glDeleteProgram(moonShaderID);
    glDeleteProgram(worldShaderID);

    glfwTerminate();
}

void EclipseMap::handleKeyPress(GLFWwindow *window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }else if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        if(!free_roam_mode){
        cameraDirection = glm::rotate(cameraDirection, 0.05f, cameraLeft);
        cameraUp = glm::cross(cameraLeft, cameraDirection);
        }else{
            cameraPosition += speed * cameraDirection;
        }
    }else if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        if(!free_roam_mode){        
        cameraDirection = glm::rotate(cameraDirection, -0.05f, cameraLeft);
        cameraUp = glm::cross(cameraLeft, cameraDirection);
        }else{
            cameraPosition -= speed * cameraDirection;
        }
    }else if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        if(!free_roam_mode){        
        cameraDirection = glm::rotate(cameraDirection, 0.05f, cameraUp);
        cameraLeft = glm::cross(cameraDirection, cameraUp);
        }else{
            cameraPosition -= speed * cameraLeft;
        }
    }else if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        if(!free_roam_mode){
        cameraDirection = glm::rotate(cameraDirection, -0.05f, cameraUp);
        cameraLeft = glm::cross(cameraDirection, cameraUp);
        }else{
            cameraPosition += speed * cameraLeft;
        }
    }else if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS) {
        speed += 0.1f; 
    }else if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS) {
        speed -= 0.1f; 
    }else if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) {
        cameraPosition = cameraStartPosition;
        cameraUp = cameraStartUp;
        cameraDirection = cameraStartDirection;
        cameraLeft = cameraStartLeft;
        speed = startSpeed;
    }else if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) {
        speed = 0.0f;
    }else if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
        heightFactor += 10;
    }else if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
        heightFactor -= 10;
    }else if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
        if(glfwGetKey(window, GLFW_KEY_P) == GLFW_RELEASE){
            if(displayFormat == displayFormatOptions::windowed){
                displayFormat = displayFormatOptions::fullScreen;
                glfwGetWindowPos(window, &screenPosX, &screenPosY);
                glfwGetWindowSize(window, &screenWidth, &screenHeight);
                GLFWmonitor* monitor = glfwGetPrimaryMonitor();
                const GLFWvidmode* mode = glfwGetVideoMode(monitor);
                glfwSetWindowMonitor(window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
            }else{
                displayFormat = displayFormatOptions::windowed;
                glfwSetWindowMonitor(window, nullptr, screenPosX, screenPosY, screenWidth, screenHeight, 0);
            }
        }
    }else if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) {
        if(glfwGetKey(window, GLFW_KEY_O) == GLFW_RELEASE){
            own_axis_rot_amount = own_axis_rot_amount > 0.f ? 0.f : 180.f/250.f;
        }
    }else if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) {
        if(glfwGetKey(window, GLFW_KEY_L) == GLFW_RELEASE) {
            free_roam_mode = !free_roam_mode;
        }
    }
}

GLFWwindow *EclipseMap::openWindow(const char *windowName, int width, int height) {
    if (!glfwInit()) {
        getchar();
        return 0;
    }

    const GLFWvidmode *mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow *window = glfwCreateWindow(width, height, windowName, NULL, NULL);
    glfwSetWindowMonitor(window, NULL, 1, 31, screenWidth, screenHeight, mode->refreshRate);

    if (window == NULL) {
        getchar();
        glfwTerminate();
        return 0;
    }

    glfwMakeContextCurrent(window);

    glewExperimental = true;
    if (glewInit() != GLEW_OK) {
        getchar();
        glfwTerminate();
        return 0;
    }

    glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
    glClearColor(0, 0, 0, 0);

    return window;
}


void EclipseMap::initColoredTexture(const char *filename, GLuint shader) {
    int width, height;
    glGenTextures(1, &textureColor);
    cout << shader << endl;
    glBindTexture(GL_TEXTURE_2D, textureColor);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;


    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
   

    imageWidth = width;
    imageHeight = height;

    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "TexColor"), 0);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}

void EclipseMap::initGreyTexture(const char *filename, GLuint shader) {

    glGenTextures(1, &textureGrey);
    glBindTexture(GL_TEXTURE_2D, textureGrey);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    int width, height;

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
  



    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "TexGrey"), 1);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}

void EclipseMap::initMoonColoredTexture(const char *filename, GLuint shader) {
    int width, height;
    glGenTextures(1, &moonTextureColor);
    cout << shader << endl;
    glBindTexture(GL_TEXTURE_2D, moonTextureColor);
    // set the texture wrapping parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                    GL_CLAMP_TO_EDGE);    // set texture wrapping to GL_REPEAT (default wrapping method)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // set texture filtering parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    unsigned char *raw_image = NULL;
    int bytes_per_pixel = 3;   /* or 1 for GRACYSCALE images */
    int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

    /* these are standard libjpeg structures for reading(decompression) */
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    /* libjpeg data structure for storing one row, that is, scanline of an image */
    JSAMPROW row_pointer[1];

    FILE *infile = fopen(filename, "rb");
    unsigned long location = 0;
    int i = 0, j = 0;

    if (!infile) {
        printf("Error opening jpeg file %s\n!", filename);
        return;
    }
    printf("Texture filename = %s\n", filename);

    /* here we set up the standard libjpeg error handler */
    cinfo.err = jpeg_std_error(&jerr);
    /* setup decompression process and source, then read JPEG header */
    jpeg_create_decompress(&cinfo);
    /* this makes the library read from infile */
    jpeg_stdio_src(&cinfo, infile);
    /* reading the image header which contains image information */
    jpeg_read_header(&cinfo, TRUE);
    /* Start decompression jpeg here */
    jpeg_start_decompress(&cinfo);

    /* allocate memory to hold the uncompressed image */
    raw_image = (unsigned char *) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    /* now actually read the jpeg into the raw buffer */
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);
    /* read one scan line at a time */
    while (cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++)
            raw_image[location++] = row_pointer[0][i];
    }

    height = cinfo.image_height;
    width = cinfo.image_width;


    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, raw_image);
   

    imageWidth = width;
    imageHeight = height;

    glGenerateMipmap(GL_TEXTURE_2D);

    glUseProgram(shader); // don't forget to activate/use the shader before setting uniforms!
    // either set it manually like so:

    glUniform1i(glGetUniformLocation(shader, "MoonTexColor"), 2);
    /* wrap up decompression, destroy objects, free pointers and close open files */
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    free(raw_image);
    fclose(infile);

}
