#version 120

uniform vec4 color;

void main() {
	//gl_FragColor = vec4(0.5, 0.0, 0.1, 1.0);
	gl_FragColor = color;
}