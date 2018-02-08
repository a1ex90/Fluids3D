#version 120

varying vec3 normal0;

uniform vec4 color;

void main() {
	//gl_FragColor = vec4(0.5, 0.0, 0.1, 1.0);
	gl_FragColor = color;
	//gl_FragColor = color * clamp(dot(-vec3(0,0,1), normal0), 0.0, 1.0);
}