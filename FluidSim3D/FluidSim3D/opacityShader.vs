#version 120

attribute vec2 position;
attribute float opacity;

varying float opacity0;


void main() {
	opacity0 = opacity;
	gl_Position = vec4(position, 0.0, 1.0);
}