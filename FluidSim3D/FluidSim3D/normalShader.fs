#version 120

varying vec3 normal0;

uniform vec4 color;

void main() {
	gl_FragColor = vec4(0.5f + 0.5f * normal0, 1.0f);
	//gl_FragColor = vec4(1.0f, 0.0f, 0.0f, 1.0f);
}