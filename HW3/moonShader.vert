#version 430

layout (location = 0) in vec3 VertexPosition;
layout (location = 1) in vec3 VertexNormal;
layout (location = 2) in vec2 VertexTex;

uniform vec3 lightPosition;
uniform vec3 cameraPosition;

uniform mat4 MVP;
uniform mat4 ModelMatrix;


out Data
{
    vec3 Position;
    vec3 Normal;
    vec2 TexCoord;
} data;



void main()
{

    vec4 in_world_vertex = ModelMatrix * vec4(VertexPosition, 1.0f);
    gl_Position = MVP * vec4(VertexPosition, 1.0f);

    
    data.Position = in_world_vertex.xyz;
    data.Normal = (ModelMatrix * vec4(VertexNormal, 0.f)).xyz;
    data.TexCoord = VertexTex;

}
