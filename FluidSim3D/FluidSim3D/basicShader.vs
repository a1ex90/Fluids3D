#version 120

attribute vec3 position;
attribute vec3 normal;

varying vec3 normal0;

uniform mat4 camera;
uniform mat4 model; 

void main() {
	gl_Position = camera * model * vec4(position, 1.0);
	normal0 = (model * vec4(normal, 0.0)).xyz;
}