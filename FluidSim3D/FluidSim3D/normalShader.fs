#version 150

uniform vec4 color;
uniform mat4 model;

uniform struct Light {
   vec3 position;
   vec3 intensities; //a.k.a the color of the light
} light;

uniform vec3 lightPos;
uniform vec3 lightIntensity;

in vec3 fragVert;
in vec3 fragNormal;

void main() {
	mat3 normalMatrix = transpose(inverse(mat3(model)));
	vec3 normal = normalize(normalMatrix * fragNormal);
	vec3 fragPos = vec3(model * vec4(fragVert, 1));

	//calculate the vector from this pixels surface to the light source
    vec3 surfaceToLight = lightPos - fragPos;

    //calculate the cosine of the angle of incidence
    float brightness = dot(normal, surfaceToLight) / (length(surfaceToLight) * length(normal));
    brightness = clamp(brightness, 0, 1);

    //calculate final color of the pixel, based on:
    // 1. The angle of incidence: brightness
    // 2. The color/intensities of the light: light.intensities
    // 3. The overall material color: color
    gl_FragColor = vec4(brightness * lightIntensity * color.xyz, color[3]);


	//gl_FragColor = vec4(0.5f + 0.5f * normal, 1.0f);
}