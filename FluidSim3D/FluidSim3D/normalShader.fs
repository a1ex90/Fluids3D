#version 150

uniform vec4 color;
uniform mat4 model;
uniform vec3 cameraPosition;

uniform float materialShininess;
uniform vec3 materialSpecularColor;

uniform struct Light {
	vec3 position;
	vec3 intensities; //a.k.a the color of the light
	float attenuation;
	float ambientCoefficient;
} light;

in vec3 fragVert;
in vec3 fragNormal;

void main() {
	vec3 normal = normalize(transpose(inverse(mat3(model))) * fragNormal);
	vec3 surfacePos = vec3(model * vec4(fragVert, 1));
	vec3 surfaceToLight = normalize(light.position - surfacePos);
	vec3 surfaceToCamera = normalize(cameraPosition - surfacePos);
   
	//AMBIENT COMPONENT
	vec3 ambient = light.ambientCoefficient * color.xyz * light.intensities;

    //DIFFUSE COMPONENT
    float diffuseCoefficient = max(0.0, dot(normal, surfaceToLight));
	vec3 diffuse = diffuseCoefficient * color.xyz * light.intensities;

	//SPECULAR COMPONENT
	float specularCoefficient = 0.0f;
    if(diffuseCoefficient > 0.0)
        specularCoefficient = pow(max(0.0, dot(surfaceToCamera, reflect(-surfaceToLight, normal))), materialShininess);
    vec3 specular = specularCoefficient * materialSpecularColor * light.intensities;

	//Attenuation
    float distanceToLight = length(light.position - surfacePos);
    float attenuation = 1.0 / (1.0 + light.attenuation * pow(distanceToLight, 2));

	vec3 linearColor = ambient + attenuation*(diffuse + specular);
    vec3 gamma = vec3(1.0/2.2);
    gl_FragColor = vec4(pow(linearColor, gamma), color[3]);
}