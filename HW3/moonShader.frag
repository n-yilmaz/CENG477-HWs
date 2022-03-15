#version 430

in Data
{
    vec3 Position;
    vec3 Normal;
    vec2 TexCoord;
} data;


uniform sampler2D MoonTexColor;
uniform vec3 lightPosition;
uniform vec3 cameraPosition;

out vec4 FragColor;

vec4 ambientReflectenceCoefficient = vec4(0.5f, 0.5f, 0.5f, 1.0f);
vec4 ambientLightColor = vec4(0.6f, 0.6f, 0.6f, 1.0f);
vec4 specularLightColor = vec4(1.0f);
float SpecularExponent = 10;



void main()
{
    vec3 texColor = texture(MoonTexColor, data.TexCoord).rgb;

    vec3 normalized_normal = normalize(data.Normal);
    vec3 LightVector = normalize(lightPosition - data.Position);
    vec3 CameraVector = normalize(cameraPosition - data.Position);

    vec3 ambient = (ambientReflectenceCoefficient * ambientLightColor).xyz;

    float diff_cos = dot(normalized_normal, LightVector);
    diff_cos = diff_cos > 0.f ? diff_cos : 0.f;
    vec3 diffuse = texColor * diff_cos;



    vec3 reflected = reflect(-LightVector, normalized_normal);
    float spec_coeff = dot(reflected, CameraVector);
    spec_coeff = diff_cos > 0.f && spec_coeff > 0.f ? spec_coeff : 0.f;
    spec_coeff = pow(spec_coeff, SpecularExponent);
    vec3 spec = specularLightColor.xyz * spec_coeff;

    vec3 combined = ambient+diffuse+spec;

    FragColor = vec4(clamp(texColor*combined, 0.f, 1.f), 1.0f);

}
