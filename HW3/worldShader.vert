#version 430

layout (location = 0) in vec3 VertexPosition;
layout (location = 1) in vec3 VertexNormal;
layout (location = 2) in vec2 VertexTex;

uniform vec3 lightPosition;
uniform vec3 cameraPosition;

uniform mat4 MVP;
uniform mat4 ModelMatrix;

uniform sampler2D TexColor;
uniform sampler2D TexGrey;

float stackCount = 125.f;
float sectorCount = 250.f;
float radius = 600.f;

float sectorStep = 0.0251327f;
float stackStep = 0.0251327f;
float piovertwo = 1.57079632679f;


uniform float heightFactor;


out Data
{
    vec3 Position;
    vec3 Normal;
    vec2 TexCoord;
} data;




vec3 findVertexPos(float sectorNum, float stackNum){
    float sectorAngle = sectorNum * sectorStep;
    float stackAngle = piovertwo - stackNum*stackStep;

    float x, y, xy, z;

    xy = radius * cos(stackAngle);
    x =  xy * cos(sectorAngle);
    y = xy * sin(sectorAngle);
    z = radius * sin(stackAngle);
    vec3 vertical_disp = normalize(vec3(x,y,z)) * texture(TexGrey, vec2(sectorNum/sectorCount, stackNum/stackCount)).r * heightFactor;
    return vertical_disp + vec3(x,y,z);
}

vec3 findNormal(vec3 real_vertex){

    if(VertexTex.t == 0 || VertexTex.t == stackCount){
        return VertexNormal;
    };

    float j = VertexTex.s;
    float i = VertexTex.t;
    if(i == 1.f){
        vec3 adj0 = vec3(0.f, 0.f, 600.f);
        adj0 += vec3(0.f, 0.f, 1.f) * texture(TexGrey, vec2(j/sectorCount, 0.f)).r * heightFactor;
        vec3 adj1, adj2, adj3, adj4;
        adj1 = findVertexPos(j+1, 1.f);
        adj2 = findVertexPos(j, 2.f);
        adj3 = findVertexPos(j-1, 2.f);
        adj4 = findVertexPos(j-1, 1.f);
        
        vec3 normal1 = cross(adj1-real_vertex, adj0-real_vertex);
        vec3 normal2 = cross(adj2-real_vertex, adj1-real_vertex);
        vec3 normal3 = cross(adj3-real_vertex, adj2-real_vertex);
        vec3 normal4 = cross(adj4-real_vertex, adj3-real_vertex);
        vec3 normal5 = cross(adj0-real_vertex, adj4-real_vertex);

        return normalize(normal1+normal2+normal3+normal4+normal5);

    }else if(i == 249.f){

        vec3 adj0 = vec3(0.f, 0.f, -600.f);
        adj0 += vec3(0.f,0.f,-1.f) * texture(TexGrey, vec2(j/sectorCount, 1.f)).r * heightFactor;
        vec3 adj1, adj2, adj3, adj4;
        adj1 = findVertexPos(j+1, 249.f);
        adj2 = findVertexPos(j+1, 248.f);
        adj3 = findVertexPos(j, 248.f);
        adj4 = findVertexPos(j-1, 249.f);
        
        vec3 normal1 = cross(adj1-real_vertex, adj2-real_vertex);
        vec3 normal2 = cross(adj2-real_vertex, adj3-real_vertex);
        vec3 normal3 = cross(adj3-real_vertex, adj4-real_vertex);
        vec3 normal4 = cross(adj4-real_vertex, adj0-real_vertex);
        vec3 normal5 = cross(adj0-real_vertex, adj1-real_vertex);

        return normalize(normal1+normal2+normal3+normal4+normal5);

    }else{

        vec3 adj0, adj1, adj2, adj3, adj4, adj5;
        adj0 = findVertexPos(j+1, i);
        adj1 = findVertexPos(j+1, i-1);
        adj2 = findVertexPos(j, i-1);
        adj3 = findVertexPos(j-1, i);
        adj4 = findVertexPos(j-1, i+1);
        adj5 = findVertexPos(j, i+1);

        vec3 normal0 = cross(adj0-real_vertex, adj1-real_vertex);
        vec3 normal1 = cross(adj1-real_vertex, adj2-real_vertex);
        vec3 normal2 = cross(adj2-real_vertex, adj3-real_vertex);
        vec3 normal3 = cross(adj3-real_vertex, adj4-real_vertex);
        vec3 normal4 = cross(adj4-real_vertex, adj5-real_vertex);
        vec3 normal5 = cross(adj5-real_vertex, adj0-real_vertex);

        return normalize(normal0+normal1+normal2+normal3+normal4+normal5);
    }




}

void main()
{
    vec3 with_height_pos;
    vec2 real_tex;
    if(VertexTex.s == 0 || VertexTex.s == sectorCount){
        real_tex = vec2(VertexTex.s / sectorCount, VertexTex.t / stackCount);
        float height_1 = texture(TexGrey, vec2(0.f, VertexTex.t / stackCount)).r;
        float height_2 = texture(TexGrey, vec2(1.f, VertexTex.t / stackCount)).r;
        vec3 heightDisp = VertexNormal * ((height_1 + height_2) / 2.f * heightFactor);
        with_height_pos = VertexPosition + heightDisp;
    }else{
        real_tex = vec2(VertexTex.s / sectorCount, VertexTex.t / stackCount);
        float height = texture(TexGrey, real_tex).r * heightFactor;
        vec3 heightDisplacement = VertexNormal * height;
        with_height_pos    = VertexPosition + heightDisplacement;
    }
    
    gl_Position = MVP * vec4(with_height_pos,1.f); 

    data.Position = with_height_pos;
    data.TexCoord = real_tex;
    vec3 in_world_normal = findNormal(with_height_pos);
    data.Normal   = (ModelMatrix * vec4(in_world_normal, 0.f)).xyz;

}
