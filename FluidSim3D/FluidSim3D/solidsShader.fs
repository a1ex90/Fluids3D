#version 120

varying vec3 normal0;

void main() {
	float opacity = 0.1;

	vec3 color = vec3(0.392, 0.584, 0.929);
	gl_FragColor = vec4(color, opacity);

	//gl_FragColor = vec4(color, opacity) * clamp(dot(-vec3(0,0,1), normal0), 0.0, 1.0);
}