#version 150

uniform vec4 color;
uniform mat4 model;
uniform vec3 cameraPosition;

uniform vec3 lightPos;
uniform vec3 lightIntensity;

in vec3 fragVert;
in vec3 fragNormal;

void main() {
	vec3 normal = normalize(transpose(inverse(mat3(model))) * fragNormal);
	vec3 surfacePos = vec3(model * vec4(fragVert, 1));
	vec3 surfaceToLight = normalize(lightPos - surfacePos);
	vec3 surfaceToCamera = normalize(cameraPosition - surfacePos);
   
	//AMBIENT COMPONENT
	float lightAmbientCoefficient = 0.1f;
	vec3 ambient = lightAmbientCoefficient * color.xyz * lightIntensity;

    //DIFFUSE COMPONENT
    float diffuseCoefficient = max(0.0, dot(normal, surfaceToLight));
	vec3 diffuse = diffuseCoefficient * color.xyz * lightIntensity;

	//SPECULAR COMPONENT
	float materialShininess = 80.0f;
	vec3 materialSpecularColor = vec3(1.0f, 1.0f, 1.0f);
	float specularCoefficient = 0.0f;
    if(diffuseCoefficient > 0.0)
        specularCoefficient = pow(max(0.0, dot(surfaceToCamera, reflect(-surfaceToLight, normal))), materialShininess);
    vec3 specular = specularCoefficient * materialSpecularColor * lightIntensity;

	//Attenuation
	float lightAttenuation = 0.2f;
    float distanceToLight = length(lightPos - surfacePos);
    float attenuation = 1.0 / (1.0 + lightAttenuation * pow(distanceToLight, 2));

	vec3 linearColor = ambient + attenuation*(diffuse + specular);
    vec3 gamma = vec3(1.0/2.2);
    gl_FragColor = vec4(pow(linearColor, gamma), color[3]);
}